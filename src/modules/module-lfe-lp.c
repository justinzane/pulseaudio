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
#include <math.h>

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
#define MEMBLOCKQ_MAXLENGTH (16*1024*1024)
#ifndef PI
#define PI 3.1415926535897932384626433
#endif
#define MAX_POLES 10
#define MIN_POLES 1
#define MIN_CUTOFF_FREQ 20.0
#define MAX_CUTOFF_FREQ 500.0

typedef struct monopole_factors {
    double a0;
} monopole_factors;

typedef struct monopole_data {
    double h0;
    double h1;
} monopole_data;

typedef struct biquad_factors {
    double a0;
    double a1;
    double a2;
    double b0;
    double b1;
    double b2;
} biquad_factors;

typedef struct biquad_data{
    double y0;
    double y1;
    double y2;
    double w0;
    double w1;
    double w2;
} biquad_data;

PA_MODULE_AUTHOR("Justin Chudgar");
PA_MODULE_DESCRIPTION(_("LFE LP Filter"));
#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.0.0-bogus"
#endif
PA_MODULE_VERSION(PACKAGE_VERSION);
PA_MODULE_LOAD_ONCE(FALSE);
PA_MODULE_USAGE( _("sink_name=<name for the sink> "
        "sink_properties=<properties for the sink> "
        "master=<name of sink to filter> "
        "lpfreq=low pass cutoff freq 50-500 Hz"
        "lppoles=<number of poles in filter>"
        "use_volume_sharing=<yes or no> "
        "force_flat_volume=<yes or no> " ));

struct userdata {
    pa_module *module;
    pa_bool_t autoloaded;
    pa_sink *sink;
    pa_sink_input *sink_input;
    pa_memblockq *memblockq;
    pa_bool_t auto_desc;
    pa_sample_spec sample_spec;
    double lpfreq;
    unsigned int poles;
    unsigned int lfeidx;
    struct biquad_factors *bqfs;
    struct biquad_data *bqdt;
    struct monopole_factors *mpfs;
    struct monopole_data *mpdt;
    //debugging
    FILE *logfile;
};

/**
 * \brief calculate biquad butterworth second order lowpass
 * \param [in]      bqfs    polynomial factors calculated before
 * \param [in, out] data    history data
 * \param [in]      src     input sample
 */
static void biquad_stage(struct biquad_factors *bqfs,
                           struct biquad_data *bqdt,
                           double src) {
    //#y0= (b0 * x0 + b1 * x1 + b2 * x2) − (a1 * y1 + a2 * y2);
    (*bqdt).y0 =     src * (*bqfs).b0 + (*bqdt).w1 * (*bqfs).b1 +
              (*bqdt).w2 * (*bqfs).b2 - ((*bqdt).y1 * (*bqfs).a1 +
              (*bqdt).y2 * (*bqfs).a2);
    (*bqdt).w2 = (*bqdt).w1;
    (*bqdt).w1 = (*bqdt).w0;
    (*bqdt).w0 = src;
    (*bqdt).y2 = (*bqdt).y1;
    (*bqdt).y1 = (*bqdt).y0;
}

/**
 * \brief calculates the first order RC lowpass
 * \param [in]      mpfs  the factor as previously calculated
 * \param [in/out]  mpdt  the history data
 * \param [in]      src   the input sample
 */
static void monopole_stage(pa_sink *sink,
                           double src) {
    struct userdata *u;
    pa_sink_assert_ref(sink);
    pa_assert_se(u = sink->userdata);

    u->mpdt->h0 = u->mpdt->h1 + u->mpfs->a0 * (src - u->mpdt->h1);
    u->mpdt->h1 = u->mpdt->h0;
}

/**
 * \brief do the low pass filtering
 * \param [in]  sink    this sink
 * \param [in]  src     the input sample
 * \param [out] dst     the filtered sample
 */
static void do_filter(pa_sink *sink,
                      float *src,
                      float *dst) {
    struct userdata         *u;
    static unsigned int     poles;
    double           sample;
    pa_sink_assert_ref(sink);
    pa_assert_se(u = sink->userdata);

    poles = u->poles;
    sample = (double)*src;

    if ((poles % 2) == 1) {
        poles -= 1;
        monopole_stage(sink, sample);
        sample = u->mpdt->h0;
    }
    while (poles >= 2) {
        poles -= 2;
        biquad_stage(u->bqfs, u->bqdt, sample);
        sample = u->bqdt->y0;
    }
    *dst = (float)sample;
}
/**
 * \brief function to calculate filter factors
 */
static void calc_filter_factors(pa_sink *sink) {
    struct userdata *u;
    double dt,rc,w0,Q,alpha;
    pa_sink_assert_ref(sink);
    pa_assert_se(u = sink->userdata);

    // for monopole
    dt = 1.0 / (double)u->sample_spec.rate;
    rc = 1.0 / (u->lpfreq * 2.0 * PI);
    u->mpfs->a0 = dt / (dt + rc);

    // for biquad
    w0 = 2.0 * PI * u->lpfreq / (double)u->sample_spec.rate;
    Q = sqrt(2.0)/2.0;
    alpha = sin(w0) / (2.0 * Q);

    u->bqfs->a0 = 1.0 + alpha;
    // note that we normalize by a0 from here
    u->bqfs->a1 = (-2.0 * cos(w0)) / u->bqfs->a0;
    u->bqfs->a2 = (1.0 - alpha) / u->bqfs->a0;
    u->bqfs->b0 = ((1.0 - cos(w0))/2.0) / u->bqfs->a0;
    u->bqfs->b1 = ((1.0 - cos(w0))) / u->bqfs->a0;
    u->bqfs->b2 = ((1.0 - cos(w0))/2.0) / u->bqfs->a0;
    u->bqfs->a0 = 1.0;

    pa_log_debug ("JZ: %s[%d]\n"
                  "\tdt=%0.8f, rc=%0.8f\n",
                  __FILE__, __LINE__, dt, rc);
    pa_log_debug ("\tw0=%0.8f, Q=%0.8f, alpha=%0.8f", w0, Q, alpha);
    pa_log_debug ("\tmonopole_factors\n"
                  "\ta0=[%0.8f]\n", u->mpfs->a0);
    pa_log_debug ("\tbiquad_factors\n"
                  "\t[b0, b1, b2]=[%0.8f, %0.8f, %0.8f]\n"
                  "\t[a0, a1, a2]=[%0.8f, %0.8f, %0.8f]\n",
                  u->bqfs->b0, u->bqfs->b1, u->bqfs->b2,
                  u->bqfs->a0, u->bqfs->a1, u->bqfs->a2
                 );

    return;
}

/**
 * \brief set filter history data to 0.0
 */
static void init_filter(pa_sink *sink) {
    struct userdata *u;
    pa_sink_assert_ref(sink);
    pa_assert_se(u = sink->userdata);

    u->bqdt->w0 = 0.0;
    u->bqdt->w1 = 0.0;
    u->bqdt->w2 = 0.0;
    u->bqdt->y0 = 0.0;
    u->bqdt->y1 = 0.0;
    u->bqdt->y2 = 0.0;

    u->mpdt->h1 = 0.0;
}

static const char* const valid_modargs[] = { "sink_name", "sink_properties",
                                             "master", "lpfreq", "lppoles",
                                             "use_volume_sharing",
                                             "force_flat_volume", NULL };

/* Called from I/O thread context */
static int sink_process_msg_cb(pa_msgobject *o,
                               int code,
                               void *data,
                               int64_t offset,
                               pa_memchunk *chunk) {
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
                    pa_bytes_to_usec(
                            pa_memblockq_get_length(
                                    u->sink_input->thread_info.render_memblockq),
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

    if (!PA_SINK_IS_LINKED(u->sink->thread_info.state) || !PA_SINK_INPUT_IS_LINKED(
            u->sink_input->thread_info.state))
        return;

    /* Just hand this one over to the master sink */
    pa_sink_input_request_rewind(
            u->sink_input,
            s->thread_info.rewind_nbytes + pa_memblockq_get_length(
                    u->memblockq),
            TRUE, FALSE, FALSE );
}

/* Called from I/O thread context */
static void sink_update_requested_latency_cb(pa_sink *s) {
    struct userdata *u;

    pa_sink_assert_ref(s);
    pa_assert_se(u = s->userdata);

    if (!PA_SINK_IS_LINKED(u->sink->thread_info.state) || !PA_SINK_INPUT_IS_LINKED(
            u->sink_input->thread_info.state))
        return;

    /* Just hand this one over to the master sink */
    pa_sink_input_set_requested_latency_within_thread(
            u->sink_input, pa_sink_get_requested_latency_within_thread(s));
}

/* Called from main context */
static void sink_set_volume_cb(pa_sink *s) {
    struct userdata *u;

    pa_sink_assert_ref(s);
    pa_assert_se(u = s->userdata);

    if (!PA_SINK_IS_LINKED(pa_sink_get_state(s)) || !PA_SINK_INPUT_IS_LINKED(
            pa_sink_input_get_state(u->sink_input)))
        return;

    pa_sink_input_set_volume(u->sink_input, &s->real_volume, s->save_volume,
                             TRUE );
}

/* Called from main context */
static void sink_set_mute_cb(pa_sink *s) {
    struct userdata *u;

    pa_sink_assert_ref(s);
    pa_assert_se(u = s->userdata);

    if (!PA_SINK_IS_LINKED(pa_sink_get_state(s)) || !PA_SINK_INPUT_IS_LINKED(
            pa_sink_input_get_state(u->sink_input)))
        return;

    pa_sink_input_set_mute(u->sink_input, s->muted, s->save_muted);
}

/**
 * \brief Low pass filters the audio data from sink_input
 * \return always returns (int)0
 * \note Called from I/O thread context
 */
static int sink_input_pop_cb(pa_sink_input *sink_input,
                             size_t nbytes,
                             pa_memchunk *chunk) {
    struct userdata *u;
    float *src, *dst, *cur_sample, *cur_frame, *dst_frame, *dst_sample;
    size_t framesize;
    unsigned num_frames, frm_idx, chan_idx;
    pa_memchunk tchunk;
    pa_usec_t current_latency PA_GCC_UNUSED;

    pa_sink_input_assert_ref(sink_input);
    pa_assert(chunk);
    pa_assert_se(u = sink_input->userdata);

    /* Hmm, process any rewind request that might be queued up */
    pa_sink_process_rewind(u->sink, 0);

    while (pa_memblockq_peek(u->memblockq, &tchunk) < 0) {
        pa_memchunk nchunk;
        pa_sink_render(u->sink, nbytes, &nchunk);
        pa_memblockq_push(u->memblockq, &nchunk);
        pa_memblock_unref(nchunk.memblock);
    }

    tchunk.length = PA_MIN(nbytes, tchunk.length);
    pa_assert(tchunk.length > 0);

    framesize = pa_frame_size(&sink_input->sample_spec);
    num_frames = (unsigned) (tchunk.length / framesize);
    pa_assert(num_frames > 0);

    chunk->index = 0;
    chunk->length = num_frames * framesize;
    chunk->memblock = pa_memblock_new(sink_input->sink->core->mempool,
                                      chunk->length);

    pa_memblockq_drop(u->memblockq, chunk->length);

    src = pa_memblock_acquire_chunk(&tchunk);
    dst = pa_memblock_acquire(chunk->memblock);


    /* (3) PUT YOUR CODE HERE TO DO SOMETHING WITH THE DATA */
    for (frm_idx = 0; frm_idx < num_frames; frm_idx++) {
        cur_frame = src + frm_idx * u->sample_spec.channels;
        dst_frame = dst + frm_idx * u->sample_spec.channels;
        for (chan_idx = 0; chan_idx < u->sample_spec.channels; chan_idx++) {
            cur_sample = cur_frame + chan_idx;
            dst_sample = dst_frame + chan_idx;
            if (chan_idx == u->lfeidx) {
                do_filter(u->sink, cur_sample, dst_sample);
            } else {
                *dst_sample = *cur_sample;
            }
        }
    }

    pa_memblock_release(tchunk.memblock);
    pa_memblock_release(chunk->memblock);

    pa_memblock_unref(tchunk.memblock);

    /* (4) IF YOU NEED THE LATENCY FOR SOMETHING ACQUIRE IT LIKE THIS: */
    /* Get the latency of master and add the latency internal to our sink input*/
    current_latency = (pa_sink_get_latency_within_thread(sink_input->sink)
            + pa_bytes_to_usec(
                    pa_memblockq_get_length(
                            sink_input->thread_info.render_memblockq),
                    &sink_input->sink->sample_spec));

    return (0);
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
            pa_memblockq_seek(u->memblockq, -(int64_t) amount, PA_SEEK_RELATIVE,
                              TRUE );
            /* (5) PUT YOUR CODE HERE TO RESET YOUR FILTER  */
            // init filter data
            calc_filter_factors(u->sink);
            init_filter(u->sink);
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
    pa_sink_set_fixed_latency_within_thread(u->sink,
                                            i->sink->thread_info.fixed_latency);
}

/* Called from I/O thread context */
static void sink_input_detach_cb(pa_sink_input *i) {
    struct userdata *u;
    pa_sink_input_assert_ref(i);
    pa_assert_se(u = i->userdata);
    pa_sink_detach_within_thread(u->sink);
    pa_sink_set_rtpoll(u->sink, NULL );
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
    pa_sink_set_fixed_latency_within_thread(u->sink,
                                            i->sink->thread_info.fixed_latency);
    pa_sink_set_max_request_within_thread(u->sink,
                                          pa_sink_input_get_max_request(i));
    /* FIXME: Too small max_rewind:
     * https://bugs.freedesktop.org/show_bug.cgi?id=53709 */
    pa_sink_set_max_rewind_within_thread(u->sink,
                                         pa_sink_input_get_max_rewind(i));

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

    pa_module_unload_request(u->module, TRUE );
}

/* Called from IO thread context */
static void sink_input_state_change_cb(pa_sink_input *i,
                                       pa_sink_input_state_t state) {
    struct userdata *u;
    pa_sink_input_assert_ref(i);
    pa_assert_se(u = i->userdata);
    /* If we are added for the first time, ask for a rewinding so that we are
     * heard right-away. */
    if (PA_SINK_INPUT_IS_LINKED(state) && i->thread_info.state
            == PA_SINK_INPUT_INIT) {
        pa_log_debug("Requesting rewind due to state change.");
        pa_sink_input_request_rewind(i, 0, FALSE, TRUE, TRUE );
    }
}

/* Called from main context */
static pa_bool_t sink_input_may_move_to_cb(pa_sink_input *i, pa_sink *dest) {
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
        pa_sink_update_flags(u->sink, PA_SINK_LATENCY | PA_SINK_DYNAMIC_LATENCY,
                             dest->flags);
    } else
        pa_sink_set_asyncmsgq(u->sink, NULL );

    if (u->auto_desc && dest) {
        const char *z;
        pa_proplist *pl;

        pl = pa_proplist_new();
        z = pa_proplist_gets(dest->proplist, PA_PROP_DEVICE_DESCRIPTION);
        pa_proplist_setf(
                pl, PA_PROP_DEVICE_DESCRIPTION, "lfe-lp %s on %s",
                pa_proplist_gets(u->sink->proplist, "device.vsink.name"),
                z ? z : dest->name);

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
    pa_sink *master = NULL;
    pa_sink_input_new_data sink_input_data;
    pa_sink_new_data sink_data;
    pa_bool_t use_volume_sharing = TRUE;
    pa_bool_t force_flat_volume = FALSE;
    pa_memchunk silence;
    struct biquad_factors bqfs =   {0.0,0.0,0.0,0.0,0.0,0.0};
    struct biquad_data bqdt =      {0.0,0.0,0.0,0.0,0.0,0.0};
    struct monopole_factors mpfs = {0.0};
    struct monopole_data mpdt =    {0.0};

    pa_assert(m);

    if (!(ma = pa_modargs_new(m->argument, valid_modargs))) {
        pa_log("Failed to parse module arguments.");
        goto fail;
    }

    if (!(master = pa_namereg_get(m->core,
                                  pa_modargs_get_value(ma, "master", NULL ),
                                  PA_NAMEREG_SINK))) {
        pa_log("Master sink not found");
        goto fail;
    }
    pa_assert(master);

    u = pa_xnew0(struct userdata, 1);
    u->bqfs = &bqfs;
    u->bqdt = &bqdt;
    u->mpfs = &mpfs;
    u->mpdt = &mpdt;

    // TODO: debugging, delete
    u->logfile = fopen("/home/justin/tmp/logfile.pcm","w");
    setvbuf(u->logfile, NULL, _IOFBF, BUFSIZ);

    // get, validate and assign lowpass cutoff freq
    u->lpfreq = atof(pa_modargs_get_value(ma, "lpfreq", "100.0"));
    if (u->lpfreq < MIN_CUTOFF_FREQ) {
        pa_log ("JZ: %s[%d] lpfreq must be between %f and %f.",
                __FILE__, __LINE__, MIN_CUTOFF_FREQ, MAX_CUTOFF_FREQ);
        u->lpfreq = MIN_CUTOFF_FREQ;
    }
    if (u->lpfreq > MAX_CUTOFF_FREQ) {
        pa_log ("JZ: %s[%d] lpfreq must be between %f and %f.",
                __FILE__, __LINE__, MIN_CUTOFF_FREQ, MAX_CUTOFF_FREQ);
        u->lpfreq = MAX_CUTOFF_FREQ;
    }
    pa_log_debug("JZ: %s[%d] lpfreq=%f\n", __FILE__, __LINE__, u->lpfreq);

    // get, validate and assign lowpass poles
    u->poles = atoi(pa_modargs_get_value(ma, "lppoles", "1"));
    if (u->poles < MIN_POLES) {
        pa_log ("JZ: %s[%d] poles must be from 1-10.", __FILE__, __LINE__);
        u->poles = 1;
    }
    if (u->poles > MAX_POLES) {
        pa_log ("JZ: %s[%d] poles must be from 1-10.", __FILE__, __LINE__);
        u->poles = 10;
    }
    pa_log_debug("JZ: %s[%d] poles=%d\n", __FILE__, __LINE__, u->poles);

    u->sample_spec = master->sample_spec;
    u->sample_spec.format = PA_SAMPLE_FLOAT32;

    map = master->channel_map;
    if (pa_modargs_get_sample_spec_and_channel_map(ma, &u->sample_spec, &map,
                                                   PA_CHANNEL_MAP_DEFAULT)
        < 0) {
        pa_log("Invalid sample format specification or channel map");
        goto fail;
    }

    // make sure we have an lfe channel
    u->lfeidx = PA_CHANNELS_MAX + 1;
    for (int i = 0; i < u->sample_spec.channels; i++) {
        if (map.map[i] == PA_CHANNEL_POSITION_LFE)
            u->lfeidx = i;
    }
    if (u->lfeidx > PA_CHANNELS_MAX) {
        pa_log("JZ: %s[%d] no lfe channel found", __FILE__, __LINE__ );
        // TODO: debugging, remove me and reenable fail
        u->lfeidx = 0;
        //goto fail;
    }
    pa_log_debug("JZ: %s[%d] lfe_index=%u", __FILE__, __LINE__, u->lfeidx);

    if (pa_modargs_get_value_boolean(ma, "use_volume_sharing",
                                     &use_volume_sharing)
        < 0) {
        pa_log("use_volume_sharing= expects a boolean argument");
        goto fail;
    }

    if (pa_modargs_get_value_boolean(ma, "force_flat_volume",
                                     &force_flat_volume)
        < 0) {
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
    if (!(sink_data.name = pa_xstrdup(
            pa_modargs_get_value(ma, "sink_name", NULL )))) {
        sink_data.name = pa_sprintf_malloc("%s.lfe_lp", master->name);
    }
    pa_sink_new_data_set_sample_spec(&sink_data, &u->sample_spec);
    pa_sink_new_data_set_channel_map(&sink_data, &map);
    pa_proplist_sets(sink_data.proplist, PA_PROP_DEVICE_MASTER_DEVICE,
                     master->name);
    pa_proplist_sets(sink_data.proplist, PA_PROP_DEVICE_CLASS, "filter");
    pa_proplist_sets(sink_data.proplist, "device.name", sink_data.name);

    if (pa_modargs_get_proplist(ma, "sink_properties", sink_data.proplist,
                                PA_UPDATE_REPLACE)
        < 0) {
        pa_log("Invalid properties");
        pa_sink_new_data_done(&sink_data);
        goto fail;
    }

    if ((u->auto_desc = !pa_proplist_contains(sink_data.proplist,
                                              PA_PROP_DEVICE_DESCRIPTION))) {
        const char *z;
        z = pa_proplist_gets(master->proplist, PA_PROP_DEVICE_DESCRIPTION);
        pa_proplist_setf(sink_data.proplist, PA_PROP_DEVICE_DESCRIPTION,
                         "lfe_lp %s on %s", sink_data.name,
                         z ? z : master->name);
    }

    // FIXME: should be PA_SINK_SHARE_VOLUME_WITH_MASTER not 0x1000000U
    u->sink = pa_sink_new(
            m->core,
            &sink_data,
            (master->flags & (PA_SINK_LATENCY | PA_SINK_DYNAMIC_LATENCY)) | (
                    use_volume_sharing ? 0x1000000U : 0));
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
        pa_sink_enable_decibel_volume(u->sink, TRUE );
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
    pa_sink_input_new_data_set_sink(&sink_input_data, master, FALSE );
    sink_input_data.origin_sink = u->sink;
    pa_proplist_setf(
            sink_input_data.proplist, PA_PROP_MEDIA_NAME,
            "lfe-lp Stream from %s",
            pa_proplist_gets(u->sink->proplist, PA_PROP_DEVICE_DESCRIPTION));
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
    u->sink_input->update_sink_latency_range =
            sink_input_update_sink_latency_range_cb;
    u->sink_input->update_sink_fixed_latency =
            sink_input_update_sink_fixed_latency_cb;
    u->sink_input->kill = sink_input_kill_cb;
    u->sink_input->attach = sink_input_attach_cb;
    u->sink_input->detach = sink_input_detach_cb;
    u->sink_input->state_change = sink_input_state_change_cb;
    u->sink_input->may_move_to = sink_input_may_move_to_cb;
    u->sink_input->moving = sink_input_moving_cb;
    u->sink_input->volume_changed =
            use_volume_sharing ? NULL : sink_input_volume_changed_cb;
    u->sink_input->mute_changed = sink_input_mute_changed_cb;
    u->sink_input->userdata = u;

    u->sink->input_to_master = u->sink_input;

    pa_sink_input_get_silence(u->sink_input, &silence);
    u->memblockq = pa_memblockq_new("module-lfe-lp memblockq", 	//name
            0, 							//start index
            MEMBLOCKQ_MAXLENGTH, 		//max length
            0, 							//target length
            &u->sample_spec, 			//sample_spec
            1, 							//prebuf
            1, 							//minreq
            0, 							//maxrewind
            &silence);					//silence
    pa_memblock_unref(silence.memblock);

    /* (9) INITIALIZE ANYTHING ELSE YOU NEED HERE */
    // init filter data
    calc_filter_factors(u->sink);
    init_filter(u->sink);

    pa_sink_put(u->sink);
    pa_sink_input_put(u->sink_input);
    pa_modargs_free(ma);
    pa_log_debug("JZ: %s[%d] finished lfe-lp.pa__init().", __FILE__, __LINE__);
    return 0;

    fail: if (ma)
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
    pa_log_debug("JZ: %s[%d] ", __FILE__, __LINE__);
    pa_assert(m);

    if (!(u = m->userdata))
        return;

    // TODO: delete, debugging
    fflush(u->logfile);
    fclose(u->logfile);

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
/*
Command line: /www/usr/fisher/helpers/mkfilter -Bu -Lp -o 2 -a 2.0833333333e-03 0.0000000000e+00
raw alpha1    =   0.0020833333
warped alpha1 =   0.0020833631
gain at dc    :   mag = 2.356080688e+04   phase =   0.0000000000 pi
gain at centre:   mag = 1.666000631e+04   phase =  -0.5000000000 pi
S-plane poles:
     -0.0092561383 + j   0.0092561383
     -0.0092561383 + j  -0.0092561383
Z-plane zeros:
     -1.0000000000 + j   0.0000000000   2 times
Z-plane poles:
      0.9907442546 + j   0.0091708587
      0.9907442546 + j  -0.0091708587
Recurrence relation:
y[n] = (  1 * x[n- 2])
     + (  2 * x[n- 1])
     + (  1 * x[n- 0])
     + ( -0.9816582826 * y[n- 2])
     + (  1.9814885091 * y[n- 1])

Ansi ``C'' Code
#define NZEROS 2
#define NPOLES 2
#define GAIN   2.356080688e+04
static float u->src_factors[NZEROS+1], u->dst_factors[NPOLES+1];
static void filterloop()
  { for (;;)
      { u->src_factors[0] = u->src_factors[1]; u->src_factors[1] = u->src_factors[2];
        u->src_factors[2] = next input value / GAIN;
        u->dst_factors[0] = u->dst_factors[1]; u->dst_factors[1] = u->dst_factors[2];
        u->dst_factors[2] =   (u->src_factors[0] + u->src_factors[2]) + 2 * u->src_factors[1]  + ( -0.9816582826 * u->dst_factors[0]) + (  1.9814885091 * u->dst_factors[1]);
        next output value = u->dst_factors[2];
      }
  }
 */
