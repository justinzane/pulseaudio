/***
    This file is part of PulseAudio.

    Copyright 2013 Justin Chudgar <justin@justinzane.com>

    PulseAudio is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published
    by the Free Software Foundation; either version 2.1 of the License,
    or (at your option) any later version.

    PulseAudio is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
    General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with PulseAudio; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
    USA.
***/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>

#include <pulse/gccmacro.h>
#include <pulse/xmalloc.h>
#include <pulse/def.h>

#include <pulsecore/i18n.h>
#include <pulsecore/namereg.h>
#include <pulsecore/sink.h>
#include <pulsecore/module.h>
#include <pulsecore/core-util.h>
#include <pulsecore/modargs.h>
#include <pulsecore/log.h>
#include <pulsecore/rtpoll.h>
#include <pulsecore/sample-util.h>
#include <pulsecore/ltdl-helper.h>

#include "module-lfe-lp-symdef.h"

PA_MODULE_AUTHOR("Justin Chudgar");
PA_MODULE_DESCRIPTION(_("LFE LP Filter"));
PA_MODULE_VERSION(PACKAGE_VERSION);
PA_MODULE_LOAD_ONCE(false);
PA_MODULE_USAGE(
        _("sink_name=<name for the sink> "
          "sink_properties=<properties for the sink> "
          "lpfreq=low pass cutoff freq 50-500 Hz"
          "master=<name of sink to filter> "
          "rate=<sample rate> "
          "channels=<number of channels> "
          "channel_map=<channel map> "
          "use_volume_sharing=<yes or no> "
          "force_flat_volume=<yes or no> "
        ));

#define MEMBLOCKQ_MAXLENGTH (16*1024*1024)

struct userdata {
    pa_module *module;
    bool autoloaded;
    pa_sink *sink;
    pa_sink_input *sink_input;
    pa_memblockq *memblockq;
    bool auto_desc;
    pa_sample_spec sample_spec;
    double lpfreq;
    unsigned int lfe_index;
    unsigned int center_index;
    double hp_factor, lp_factor;
    double prev_sample[PA_CHANNELS_MAX];
    double filtered[PA_CHANNELS_MAX];
    // TODO: debugging, remove me
    FILE *lfe_lp_log;
};

static const char* const valid_modargs[] = {
    "sink_name",
    "sink_properties",
    "lpfreq",
    "master",
    "rate",
    "channels",
    "channel_map",
    "use_volume_sharing",
    "force_flat_volume",
    NULL
};

/* Called from I/O thread context */
static int sink_process_msg_cb(pa_msgobject *o, int code, void *data, int64_t offset, pa_memchunk *chunk) {
    struct userdata *u = PA_SINK(o)->userdata;

    switch (code) {
        case PA_SINK_MESSAGE_GET_LATENCY:
            /* The sink is _put() before the sink input is, so let's make sure
             * we don't access it in that time. Also, the sink input is first
             * shut down, the sink second. */
            if (!PA_SINK_IS_LINKED(u->sink->thread_info.state) ||
                !PA_SINK_INPUT_IS_LINKED(u->sink_input->thread_info.state)) {
                *((pa_usec_t*) data) = 0;
                return 0;
            }
            /* Get the latency of the master sink and add the latency internal
               to our sink input on top */
            *((pa_usec_t*) data) =
                pa_sink_get_latency_within_thread(u->sink_input->sink) +
                pa_bytes_to_usec(pa_memblockq_get_length(u->sink_input->thread_info.render_memblockq),
                                 &u->sink_input->sink->sample_spec);
            return 0;
    }
    return pa_sink_process_msg(o, code, data, offset, chunk);
}

/* Called from main context */
static int sink_set_state_cb(pa_sink *s, pa_sink_state_t state) {
    struct userdata *u;

    pa_sink_assert_ref(s);
    pa_assert_se(u = s->userdata);

    if (!PA_SINK_IS_LINKED(state) ||
        !PA_SINK_INPUT_IS_LINKED(pa_sink_input_get_state(u->sink_input)))
        return 0;

    pa_sink_input_cork(u->sink_input, state == PA_SINK_SUSPENDED);
    return 0;
}

/* Called from I/O thread context */
static void sink_request_rewind_cb(pa_sink *s) {
    struct userdata *u;

    pa_sink_assert_ref(s);
    pa_assert_se(u = s->userdata);

    if (!PA_SINK_IS_LINKED(u->sink->thread_info.state) ||
        !PA_SINK_INPUT_IS_LINKED(u->sink_input->thread_info.state))
        return;

    /* Just hand this one over to the master sink */
    pa_sink_input_request_rewind(u->sink_input,
                                 s->thread_info.rewind_nbytes +
                                 pa_memblockq_get_length(u->memblockq), true, false, false);
}

/* Called from I/O thread context */
static void sink_update_requested_latency_cb(pa_sink *s) {
    struct userdata *u;

    pa_sink_assert_ref(s);
    pa_assert_se(u = s->userdata);

    if (!PA_SINK_IS_LINKED(u->sink->thread_info.state) ||
        !PA_SINK_INPUT_IS_LINKED(u->sink_input->thread_info.state))
        return;

    /* Just hand this one over to the master sink */
    pa_sink_input_set_requested_latency_within_thread(
            u->sink_input,
            pa_sink_get_requested_latency_within_thread(s));
}

/* Called from main context */
static void sink_set_volume_cb(pa_sink *s) {
    struct userdata *u;

    pa_sink_assert_ref(s);
    pa_assert_se(u = s->userdata);

    if (!PA_SINK_IS_LINKED(pa_sink_get_state(s)) ||
        !PA_SINK_INPUT_IS_LINKED(pa_sink_input_get_state(u->sink_input)))
        return;

    pa_sink_input_set_volume(u->sink_input, &s->real_volume, s->save_volume, true);
}

/* Called from main context */
static void sink_set_mute_cb(pa_sink *s) {
    struct userdata *u;

    pa_sink_assert_ref(s);
    pa_assert_se(u = s->userdata);

    if (!PA_SINK_IS_LINKED(pa_sink_get_state(s)) ||
        !PA_SINK_INPUT_IS_LINKED(pa_sink_input_get_state(u->sink_input)))
        return;

    pa_sink_input_set_mute(u->sink_input, s->muted, s->save_muted);
}

/**
 * \fn sink_input_pop_cb
 * \return int: What is the significance of the return value?
 * \param sink_input
 * \param nbytes
 * \param chunk
 * Called from I/O thread context
 * Low pass filters the LFE channel and high pass filters the remaining channels.
 */
static int sink_input_pop_cb(pa_sink_input *sink_input,
                             size_t nbytes,
                             pa_memchunk *chunk) {
    struct userdata *u;
    float *src, *dst, *cur_sample, *cur_frame, *dst_frame, *dst_sample;
    size_t framesize;
    unsigned num_frames, frame_index, channel_index;
    pa_memchunk tchunk;
    pa_usec_t current_latency PA_GCC_UNUSED;

    pa_sink_input_assert_ref(sink_input);
    pa_assert(chunk);
    pa_assert_se(u = sink_input->userdata);

    /* Hmm, process any rewind request that might be queued up */
    pa_sink_process_rewind(u->sink, 0);

    /* (1) IF YOU NEED A FIXED BLOCK SIZE USE
     * pa_memblockq_peek_fixed_size() HERE INSTEAD. NOTE THAT FILTERS
     * WHICH CAN DEAL WITH DYNAMIC BLOCK SIZES ARE HIGHLY
     * PREFERRED. */
    while (pa_memblockq_peek(u->memblockq, &tchunk) < 0) {
        pa_memchunk nchunk;
        pa_sink_render(u->sink, nbytes, &nchunk);
        pa_memblockq_push(u->memblockq, &nchunk);
        pa_memblock_unref(nchunk.memblock);
    }

    /* (2) IF YOU NEED A FIXED BLOCK SIZE, THIS NEXT LINE IS NOT NECESSARY */
    tchunk.length = PA_MIN(nbytes, tchunk.length);
    pa_assert(tchunk.length > 0);
    //pa_log_debug ("JZ: %s[%d] tchunk.length=%d\n", __FILE__, __LINE__, tchunk.length);

    framesize = pa_frame_size(&sink_input->sample_spec);
    //pa_log_debug ("JZ: %s[%d] framesize=%d", __FILE__, __LINE__, framesize);
    num_frames = (unsigned) (tchunk.length / framesize);
    //pa_log_debug ("JZ: %s[%d] num_frames=%d", __FILE__, __LINE__, num_frames);

    pa_assert(num_frames > 0);

    chunk->index = 0;
    chunk->length = num_frames*framesize;
    chunk->memblock = pa_memblock_new(sink_input->sink->core->mempool, chunk->length);

    pa_memblockq_drop(u->memblockq, chunk->length);

    src = pa_memblock_acquire_chunk(&tchunk);
    dst = pa_memblock_acquire(chunk->memblock);

    /* (3) PUT YOUR CODE HERE TO DO SOMETHING WITH THE DATA */

    for (frame_index = 0; frame_index < num_frames; frame_index++) {
        cur_frame = src + frame_index * u->sample_spec.channels;
        dst_frame = dst + frame_index * u->sample_spec.channels;
        for (channel_index = 0; channel_index < u->sample_spec.channels; channel_index++) {
            cur_sample = cur_frame + channel_index;
            dst_sample = dst_frame + channel_index;
            if (channel_index == u->lfe_index) {                    // low pass
                u->filtered[channel_index] += (u->lp_factor *
                                               ((double)*cur_sample -
                                                u->filtered[channel_index]));
            /*} else if (channel_index == u->center_index) {            // high pass
                u->filtered[channel_index] = (u->hp_factor *
                                              (u->filtered[channel_index] +
                                               (double)*cur_sample -
                                               u->prev_sample[channel_index]));
                u->prev_sample[channel_index] = (double)*cur_sample; */
            } else {
                u->filtered[channel_index] = (double)*cur_sample;
                u->prev_sample[channel_index] = (double)*cur_sample;
            }
            *dst_sample = (float)u->filtered[channel_index];
            /* TODO: debugging, remove me
            fprintf (u->lfe_lp_log,
                     "%u,%+01.12f,%+01.12f\n",
                     channel_index,
                     u->filtered[channel_index],
                     (double)*cur_sample);
            */
        }
    }

    pa_memblock_release(tchunk.memblock);
    pa_memblock_release(chunk->memblock);

    pa_memblock_unref(tchunk.memblock);

    /* (4) IF YOU NEED THE LATENCY FOR SOMETHING ACQUIRE IT LIKE THIS: */
    /* Get the latency of master and add the latency internal to our sink input*/
    current_latency = (pa_sink_get_latency_within_thread(sink_input->sink) +
                       pa_bytes_to_usec(pa_memblockq_get_length(sink_input->thread_info.render_memblockq),
                                        &sink_input->sink->sample_spec)
                      );

    return 0;
}

/* Called from I/O thread context */
static void sink_input_process_rewind_cb(pa_sink_input *i, size_t nbytes) {
    struct userdata *u;
    size_t amount = 0;

    pa_sink_input_assert_ref(i);
    pa_assert_se(u = i->userdata);

    if (u->sink->thread_info.rewind_nbytes > 0) {
        size_t max_rewrite;

        max_rewrite = nbytes + pa_memblockq_get_length(u->memblockq);
        amount = PA_MIN(u->sink->thread_info.rewind_nbytes, max_rewrite);
        u->sink->thread_info.rewind_nbytes = 0;

        if (amount > 0) {
            pa_memblockq_seek(u->memblockq, - (int64_t) amount, PA_SEEK_RELATIVE, true);
            /* (5) PUT YOUR CODE HERE TO RESET YOUR FILTER  */
            // TODO: Determine how to reset.
        }
    }

    pa_sink_process_rewind(u->sink, amount);
    pa_memblockq_rewind(u->memblockq, nbytes);
}

/* Called from I/O thread context */
static void sink_input_update_max_rewind_cb(pa_sink_input *i, size_t nbytes) {
    struct userdata *u;
    pa_sink_input_assert_ref(i);
    pa_assert_se(u = i->userdata);
    /* FIXME: Too small max_rewind:
     * https://bugs.freedesktop.org/show_bug.cgi?id=53709 */
    pa_memblockq_set_maxrewind(u->memblockq, nbytes);
    pa_sink_set_max_rewind_within_thread(u->sink, nbytes);
}

/* Called from I/O thread context */
static void sink_input_update_max_request_cb(pa_sink_input *i, size_t nbytes) {
    struct userdata *u;
    pa_sink_input_assert_ref(i);
    pa_assert_se(u = i->userdata);
    /* (6) IF YOU NEED A FIXED BLOCK SIZE ROUND nbytes UP TO MULTIPLES
     * OF IT HERE. THE PA_ROUND_UP MACRO IS USEFUL FOR THAT. */
    pa_sink_set_max_request_within_thread(u->sink, nbytes);
}

/* Called from I/O thread context */
static void sink_input_update_sink_latency_range_cb(pa_sink_input *i) {
    struct userdata *u;
    pa_sink_input_assert_ref(i);
    pa_assert_se(u = i->userdata);
    pa_sink_set_latency_range_within_thread(u->sink,
                                            i->sink->thread_info.min_latency,
                                            i->sink->thread_info.max_latency);
}

/* Called from I/O thread context */
static void sink_input_update_sink_fixed_latency_cb(pa_sink_input *i) {
    struct userdata *u;
    pa_sink_input_assert_ref(i);
    pa_assert_se(u = i->userdata);
    /* (7) IF YOU NEED A FIXED BLOCK SIZE ADD THE LATENCY FOR ONE BLOCK MINUS
     * ONE SAMPLE HERE. pa_usec_to_bytes_round_up() IS USEFUL FOR THAT. */
    pa_sink_set_fixed_latency_within_thread(u->sink, i->sink->thread_info.fixed_latency);
}

/* Called from I/O thread context */
static void sink_input_detach_cb(pa_sink_input *i) {
    struct userdata *u;
    pa_sink_input_assert_ref(i);
    pa_assert_se(u = i->userdata);
    pa_sink_detach_within_thread(u->sink);
    pa_sink_set_rtpoll(u->sink, NULL);
}

/* Called from I/O thread context */
static void sink_input_attach_cb(pa_sink_input *i) {
    struct userdata *u;
    pa_sink_input_assert_ref(i);
    pa_assert_se(u = i->userdata);

    pa_sink_set_rtpoll(u->sink, i->sink->thread_info.rtpoll);
    pa_sink_set_latency_range_within_thread(u->sink,
                                            i->sink->thread_info.min_latency,
                                            i->sink->thread_info.max_latency);
    /* (8.1) IF YOU NEED A FIXED BLOCK SIZE ADD THE LATENCY FOR ONE BLOCK MINUS
     * ONE SAMPLE HERE. SEE (7) */
    pa_sink_set_fixed_latency_within_thread(u->sink, i->sink->thread_info.fixed_latency);
    /* (8.2) IF YOU NEED A FIXED BLOCK SIZE ROUND pa_sink_input_get_max_request(i)
     * UP TO MULTIPLES OF IT HERE. SEE (6) */
    pa_sink_set_max_request_within_thread(u->sink, pa_sink_input_get_max_request(i));
    /* FIXME: Too small max_rewind:
     * https://bugs.freedesktop.org/show_bug.cgi?id=53709 */
    pa_sink_set_max_rewind_within_thread(u->sink, pa_sink_input_get_max_rewind(i));

    pa_sink_attach_within_thread(u->sink);
}

/* Called from main context */
static void sink_input_kill_cb(pa_sink_input *i) {
    struct userdata *u;

    pa_sink_input_assert_ref(i);
    pa_assert_se(u = i->userdata);

    /* The order here matters! We first kill the sink input, followed by the sink.
     * That means the sink callbacks must be protected against an unconnected sink input! */
    pa_sink_input_unlink(u->sink_input);
    pa_sink_unlink(u->sink);

    pa_sink_input_unref(u->sink_input);
    u->sink_input = NULL;

    pa_sink_unref(u->sink);
    u->sink = NULL;

    pa_module_unload_request(u->module, true);
}

/* Called from IO thread context */
static void sink_input_state_change_cb(pa_sink_input *i, pa_sink_input_state_t state) {
    struct userdata *u;
    pa_sink_input_assert_ref(i);
    pa_assert_se(u = i->userdata);
    /* If we are added for the first time, ask for a rewinding so that we are heard right-away. */
    if (PA_SINK_INPUT_IS_LINKED(state) &&
        i->thread_info.state == PA_SINK_INPUT_INIT) {
        pa_log_debug("Requesting rewind due to state change.");
        pa_sink_input_request_rewind(i, 0, false, true, true);
    }
}

/* Called from main context */
static bool sink_input_may_move_to_cb(pa_sink_input *i, pa_sink *dest) {
    struct userdata *u;
    pa_sink_input_assert_ref(i);
    pa_assert_se(u = i->userdata);
    return u->sink != dest;
}

/* Called from main context */
static void sink_input_moving_cb(pa_sink_input *i, pa_sink *dest) {
    struct userdata *u;

    pa_sink_input_assert_ref(i);
    pa_assert_se(u = i->userdata);

    if (dest) {
        pa_sink_set_asyncmsgq(u->sink, dest->asyncmsgq);
        pa_sink_update_flags(u->sink, PA_SINK_LATENCY|PA_SINK_DYNAMIC_LATENCY, dest->flags);
    } else
        pa_sink_set_asyncmsgq(u->sink, NULL);

    if (u->auto_desc && dest) {
        const char *z;
        pa_proplist *pl;

        pl = pa_proplist_new();
        z = pa_proplist_gets(dest->proplist, PA_PROP_DEVICE_DESCRIPTION);
        pa_proplist_setf(pl, PA_PROP_DEVICE_DESCRIPTION, "lfe-lp %s on %s",
                         pa_proplist_gets(u->sink->proplist, "device.vsink.name"), z ? z : dest->name);

        pa_sink_update_proplist(u->sink, PA_UPDATE_REPLACE, pl);
        pa_proplist_free(pl);
    }
}

/* Called from main context */
static void sink_input_volume_changed_cb(pa_sink_input *i) {
    struct userdata *u;
    pa_sink_input_assert_ref(i);
    pa_assert_se(u = i->userdata);
    pa_sink_volume_changed(u->sink, &i->volume);
}

/* Called from main context */
static void sink_input_mute_changed_cb(pa_sink_input *i) {
    struct userdata *u;
    pa_sink_input_assert_ref(i);
    pa_assert_se(u = i->userdata);
    pa_sink_mute_changed(u->sink, i->muted);
}

int pa__init(pa_module *m) {
    struct userdata *u;
    pa_channel_map map;
    pa_modargs *ma;
    pa_sink *master=NULL;
    pa_sink_input_new_data sink_input_data;
    pa_sink_new_data sink_data;
    bool use_volume_sharing = true;
    bool force_flat_volume = false;
    float arg_lpfreq = 100.0;
    double  rc_lp, rc_hp, dt;
    pa_memchunk silence;

    pa_assert(m);

    if (!(ma = pa_modargs_new(m->argument, valid_modargs))) {
        pa_log("Failed to parse module arguments.");
        goto fail;
    }

    if (!(master = pa_namereg_get(m->core, pa_modargs_get_value(ma, "master", NULL), PA_NAMEREG_SINK))) {
        pa_log("Master sink not found");
        goto fail;
    }
    pa_assert(master);

    u = pa_xnew0(struct userdata, 1);

    /* TODO: debugging, remove me
    u->lfe_lp_log = fopen("/tmp/lfe_lp_log","w"); */

    // get, validate and assign lowpass cutoff freq
    arg_lpfreq = atof(pa_modargs_get_value(ma, "lpfreq", "100.0"));
    if ((arg_lpfreq < 20.0) || (arg_lpfreq > 500.0)) {
        pa_log ("JZ: %s[%d] lpfreq must be between 20.0 and 500.0.", __FILE__, __LINE__);
        goto fail;
    }
    u->lpfreq = arg_lpfreq;
    pa_log_debug ("JZ: %s[%d] lpfreq=%f\n", __FILE__, __LINE__, u->lpfreq );

    u->sample_spec = master->sample_spec;
    u->sample_spec.format = PA_SAMPLE_FLOAT32;

    // setup filter factors
    rc_lp = 1.0 / (u->lpfreq * 2.0 * 3.1415926535897932384626433);
    rc_hp = 1.0 / (u->lpfreq * 2.0 * 3.1415926535897932384626433);
    dt = (1.0 / (double)u->sample_spec.rate);
    u->lp_factor = dt / (dt + rc_lp);
    u->hp_factor = rc_hp / (dt + rc_hp);
    pa_log ("JZ: %s[%d]\n\tlp_factor=%f\n\thp_factor=%f\n",
                  __FILE__, __LINE__, u->lp_factor, u->hp_factor );

    // init filter data
    for (unsigned int i=0; i<PA_CHANNELS_MAX; i++) {
        u->filtered[i] = (double)0.0;
        u->prev_sample[i] = (double)0.0;
    }

    map = master->channel_map;
    if (pa_modargs_get_sample_spec_and_channel_map(ma, &u->sample_spec, &map, PA_CHANNEL_MAP_DEFAULT) < 0) {
        pa_log("Invalid sample format specification or channel map");
        goto fail;
    }

    // make sure we have and lfe channel
    u->lfe_index = PA_CHANNELS_MAX + 1;
    u->center_index = PA_CHANNELS_MAX + 2;
    for (int i=0; i<u->sample_spec.channels; i++) {
        if (map.map[i] == PA_CHANNEL_POSITION_LFE) {
            u->lfe_index = i;
            pa_log ("JZ: %s[%d] lfe_index=%u", __FILE__, __LINE__, u->lfe_index);
        }
        if (map.map[i] == PA_CHANNEL_POSITION_CENTER) {
            u->center_index = i;
            pa_log ("JZ: %s[%d] center_index=%u", __FILE__, __LINE__, u->center_index);
        }
    }
    if (u->lfe_index > PA_CHANNELS_MAX) {
        pa_log("JZ: %s[%d] no lfe channel found", __FILE__, __LINE__ );
        // TODO: debugging, remove me and reenable fail
        //u->lfe_index = 0;
        //goto fail;
    }

    if (pa_modargs_get_value_boolean(ma, "use_volume_sharing", &use_volume_sharing) < 0) {
        pa_log("use_volume_sharing= expects a boolean argument");
        goto fail;
    }

    if (pa_modargs_get_value_boolean(ma, "force_flat_volume", &force_flat_volume) < 0) {
        pa_log("force_flat_volume= expects a boolean argument");
        goto fail;
    }

    if (use_volume_sharing && force_flat_volume) {
        pa_log("Flat volume can't be forced when using volume sharing.");
        goto fail;
    }

    u->module = m;
    m->userdata = u;

    /* Create sink */
    pa_sink_new_data_init(&sink_data);
    sink_data.driver = __FILE__;
    sink_data.module = m;
    if (!(sink_data.name =
            pa_xstrdup(pa_modargs_get_value(ma, "sink_name", NULL)))) {
        sink_data.name = pa_sprintf_malloc("%s.lfe_lp", master->name);
    }
    pa_sink_new_data_set_sample_spec(&sink_data, &u->sample_spec);
    pa_sink_new_data_set_channel_map(&sink_data, &map);
    pa_proplist_sets(sink_data.proplist,
                     PA_PROP_DEVICE_MASTER_DEVICE,
                     master->name);
    pa_proplist_sets(sink_data.proplist,
                     PA_PROP_DEVICE_CLASS,
                     "filter");
    pa_proplist_sets(sink_data.proplist,
                     "device.lfe_lp.name",
                     sink_data.name);

    if (pa_modargs_get_proplist(ma,
                                "sink_properties",
                                sink_data.proplist,
                                PA_UPDATE_REPLACE) < 0) {
        pa_log("Invalid properties");
        pa_sink_new_data_done(&sink_data);
        goto fail;
    }

    if ((u->auto_desc = !pa_proplist_contains(sink_data.proplist,
                                              PA_PROP_DEVICE_DESCRIPTION))) {
        const char *z;
        z = pa_proplist_gets(master->proplist, PA_PROP_DEVICE_DESCRIPTION);
        pa_proplist_setf(sink_data.proplist,
                         PA_PROP_DEVICE_DESCRIPTION,
                         "lfe-lp %s on %s",
                         sink_data.name,
                         z ? z : master->name);
    }

    // FIXME: should be PA_SINK_SHARE_VOLUME_WITH_MASTER not 0x1000000U
    u->sink = pa_sink_new(m->core, &sink_data, (master->flags & (PA_SINK_LATENCY|PA_SINK_DYNAMIC_LATENCY))
                          | (use_volume_sharing ? 0x1000000U : 0));
    pa_sink_new_data_done(&sink_data);

    if (!u->sink) {
        pa_log("Failed to create sink.");
        goto fail;
    }

    u->sink->parent.process_msg = sink_process_msg_cb;
    u->sink->set_state = sink_set_state_cb;
    u->sink->update_requested_latency = sink_update_requested_latency_cb;
    u->sink->request_rewind = sink_request_rewind_cb;
    pa_sink_set_set_mute_callback(u->sink, sink_set_mute_cb);
    if (!use_volume_sharing) {
        pa_sink_set_set_volume_callback(u->sink, sink_set_volume_cb);
        pa_sink_enable_decibel_volume(u->sink, true);
    }
    /* Normally this flag would be enabled automatically be we can force it. */
    if (force_flat_volume)
        u->sink->flags |= PA_SINK_FLAT_VOLUME;
    u->sink->userdata = u;

    pa_sink_set_asyncmsgq(u->sink, master->asyncmsgq);

    /* Create sink input */
    pa_sink_input_new_data_init(&sink_input_data);
    sink_input_data.driver = __FILE__;
    sink_input_data.module = m;
    pa_sink_input_new_data_set_sink(&sink_input_data, master, false);
    sink_input_data.origin_sink = u->sink;
    pa_proplist_setf(sink_input_data.proplist,
                     PA_PROP_MEDIA_NAME,
                     "lfe-lp Stream from %s",
                     pa_proplist_gets(u->sink->proplist,
                                      PA_PROP_DEVICE_DESCRIPTION)
                    );
    pa_proplist_sets(sink_input_data.proplist, PA_PROP_MEDIA_ROLE, "filter");
    pa_sink_input_new_data_set_sample_spec(&sink_input_data, &u->sample_spec);
    pa_sink_input_new_data_set_channel_map(&sink_input_data, &map);

    pa_sink_input_new(&u->sink_input, m->core, &sink_input_data);
    pa_sink_input_new_data_done(&sink_input_data);

    if (!u->sink_input)
        goto fail;

    u->sink_input->pop = sink_input_pop_cb;
    u->sink_input->process_rewind = sink_input_process_rewind_cb;
    u->sink_input->update_max_rewind = sink_input_update_max_rewind_cb;
    u->sink_input->update_max_request = sink_input_update_max_request_cb;
    u->sink_input->update_sink_latency_range = sink_input_update_sink_latency_range_cb;
    u->sink_input->update_sink_fixed_latency = sink_input_update_sink_fixed_latency_cb;
    u->sink_input->kill = sink_input_kill_cb;
    u->sink_input->attach = sink_input_attach_cb;
    u->sink_input->detach = sink_input_detach_cb;
    u->sink_input->state_change = sink_input_state_change_cb;
    u->sink_input->may_move_to = sink_input_may_move_to_cb;
    u->sink_input->moving = sink_input_moving_cb;
    u->sink_input->volume_changed = use_volume_sharing ? NULL : sink_input_volume_changed_cb;
    u->sink_input->mute_changed = sink_input_mute_changed_cb;
    u->sink_input->userdata = u;

    u->sink->input_to_master = u->sink_input;

    pa_sink_input_get_silence(u->sink_input, &silence);
    u->memblockq = pa_memblockq_new("module-lfe-lp memblockq",     //name
                                    0,                             //start index
                                    MEMBLOCKQ_MAXLENGTH,         //max length
                                    0,                             //target length
                                    &u->sample_spec,             //sample_spec
                                    1,                             //prebuf
                                    1,                             //minreq
                                    0,                             //maxrewind
                                    &silence);                    //silence
    pa_memblock_unref(silence.memblock);

    /* (9) INITIALIZE ANYTHING ELSE YOU NEED HERE */
    pa_sink_put(u->sink);
    pa_sink_input_put(u->sink_input);
    pa_modargs_free(ma);
    pa_log_debug ("JZ: %s[%d] finished lfe-lp.pa__init().", __FILE__, __LINE__ );
    return 0;

fail:
    if (ma)
        pa_modargs_free(ma);
    pa__done(m);
    return -1;
}

int pa__get_n_used(pa_module *m) {
    struct userdata *u;
    pa_assert(m);
    pa_assert_se(u = m->userdata);
    return pa_sink_linked_by(u->sink);
}

void pa__done(pa_module*m) {
    struct userdata *u;
    pa_log_debug ("JZ: %s[%d] ", __FILE__, __LINE__ );
    pa_assert(m);

    if (!(u = m->userdata))
        return;

    /* TODO: debugging, remove me
    fflush(u->lfe_lp_log);
    fclose(u->lfe_lp_log); */

    /* See comments in sink_input_kill_cb() above regarding
     * destruction order! */

    if (u->sink_input)
        pa_sink_input_unlink(u->sink_input);

    if (u->sink)
        pa_sink_unlink(u->sink);

    if (u->sink_input)
        pa_sink_input_unref(u->sink_input);

    if (u->sink)
        pa_sink_unref(u->sink);

    if (u->memblockq)
        pa_memblockq_free(u->memblockq);

    pa_xfree(u);
}
