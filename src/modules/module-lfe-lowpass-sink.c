/**
 * \file        module-lfe-lp.c
 * \date        Apr 2, 2013
 * \author      Justin Chudgar, justin@justinzane.com
 * \copyright   Justin Chudgar
 * \license     GPLv3
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
 */

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

#include <pulsecore/biquad-filter.h>
#include "module-lfe-lowpass-sink-symdef.h"

PA_MODULE_AUTHOR("Justin Chudgar");
PA_MODULE_DESCRIPTION(_("LFE Lowpass Filter Sink"));
#ifndef PACKAGE_VERSION     /* Stop Eclipse from complaining */
#define PACKAGE_VERSION ("0.0.0-bogus")
#endif
PA_MODULE_VERSION(PACKAGE_VERSION);
PA_MODULE_LOAD_ONCE(FALSE);
PA_MODULE_USAGE(_("sink_name=<name for the sink> "
                  "sink_properties=<properties for the sink> "
                  "master=<name of sink to filter> "
                  "lpfreq=low pass cutoff freq 20-500 Hz"));

static const char* const valid_modargs[] = {"sink_name",
                                            "sink_properties",
                                            "master",
                                            "lpfreq",
                                            NULL };

/* persistent user data structure */
struct userdata {
    pa_module *module;
    pa_bool_t autoloaded;
    pa_sink *sink;
    pa_sink_input *sink_input;
    pa_memblockq *memblockq;
    pa_bool_t auto_desc;
    pa_sample_spec sample_spec;
    size_t sz_smp, sz_frm, sz_bqf;
    double lpfreq;                               /* corner/cutoff frequency, user defined */
//    biquad_factors *s1lpfs, *s1hpfs, *s1apfs;    /* lowpass, highpass and allpass coefficients */
//    biquad_factors *s2lpfs, *s2hpfs, *s2apfs;    /* lowpass, highpass and allpass coefficients */
//    biquad_data *s1lpdt, *s1hpdt, *s1apdt;       /* history data for the various filters stage 1*/
//    biquad_history *s1histbuf;                   /* rewind buffer for biquad_data stage 1*/
//    biquad_data *s2lpdt, *s2hpdt, *s2apdt;       /* history data for the various filters stage 2*/
//    biquad_history *s2histbuf;                   /* rewind buffer for biquad_data stage 2*/
//    biquad_types filter_map[PA_CHANNELS_MAX];    /* map of channels index to biquad_types */
    biquad_filter_map_4 *filter_map;
};

static int sink_process_msg_cb(pa_msgobject *msgobject, int code, void *data, int64_t offset, pa_memchunk *chunk) {
    struct userdata *u = PA_SINK(msgobject)->userdata;

    switch (code) {
        case PA_SINK_MESSAGE_GET_LATENCY:
            /* The sink is _put() before the sink input is, so let's make sure
             * we don't access it in that time. Also, the sink input is first
             * shut down, the sink second. */
            if (!PA_SINK_IS_LINKED(u->sink->thread_info.state) || !PA_SINK_INPUT_IS_LINKED(u->sink_input->thread_info.state)) {
                * ((pa_usec_t*)data) = 0;
                return 0;
            }
            /* Get the latency of the master sink and add the latency internal
             to our sink input on top */
            * ((pa_usec_t*)data) = pa_sink_get_latency_within_thread(u->sink_input->sink)
                    + pa_bytes_to_usec(pa_memblockq_get_length(u->sink_input->thread_info.render_memblockq),
                                       &u->sink_input->sink->sample_spec);
            return 0;
    }
    return pa_sink_process_msg(msgobject, code, data, offset, chunk);
}

static int sink_set_state_cb(pa_sink *sink, pa_sink_state_t state) {
    struct userdata *u;
    pa_sink_assert_ref(sink);
    pa_assert_se(u = sink->userdata);

    if (!PA_SINK_IS_LINKED(state) || !PA_SINK_INPUT_IS_LINKED(pa_sink_input_get_state(u->sink_input)))
        return 0;

    pa_sink_input_cork(u->sink_input, state == PA_SINK_SUSPENDED);
    return 0;
}

static void sink_request_rewind_cb(pa_sink *sink) {
    struct userdata *u;
    pa_sink_assert_ref(sink);
    pa_assert_se(u = sink->userdata);

    if (!PA_SINK_IS_LINKED(u->sink->thread_info.state) || !PA_SINK_INPUT_IS_LINKED(u->sink_input->thread_info.state)) {
        pa_log_warn("[%d]%s sink or sink-input not linked, cannot rewind", __LINE__, __func__);
        return;
    }

    /* Just hand this one over to the master sink */
    pa_sink_input_request_rewind(u->sink_input, (sink->thread_info.rewind_nbytes + pa_memblockq_get_length(u->memblockq)), TRUE,
                                 FALSE, FALSE );
}

static void sink_update_requested_latency_cb(pa_sink *sink) {
    struct userdata *u;
    pa_sink_assert_ref(sink);
    pa_assert_se(u = sink->userdata);

    if (!PA_SINK_IS_LINKED(u->sink->thread_info.state) || !PA_SINK_INPUT_IS_LINKED(u->sink_input->thread_info.state))
        return;

    /* Just hand this one over to the master sink */
    pa_sink_input_set_requested_latency_within_thread(u->sink_input, pa_sink_get_requested_latency_within_thread(sink));
}

/* filters the audio data from sink_input
 * \return always returns (int)0
 * \note Called from I/O thread context */
__attribute__((optimize(3))) static int sink_input_pop_cb(pa_sink_input *sink_input,
                             size_t nbytes,
                             pa_memchunk *dst_chunk) {
    struct userdata *u;
    pa_memchunk *src, *dst;
//  float *cur_sample, *cur_frame, *dst_frame, *dst_sample;
//    float hp = 0.0f;
//    float lp = 0.0f;
//    float ap = 0.0f;
    size_t framesize, num_frames, chan_idx; //, frm_idx
//    biquad_data_element bqdtel;
    pa_memchunk tchunk;
    pa_memchunk deint_chunk, reint_chunk;
    biquad_map_item_4 *current_map_item;

    pa_sink_input_assert_ref(sink_input);
    pa_assert(dst_chunk);
    pa_assert_se(u = sink_input->userdata);

    /* Removed: https://bugs.freedesktop.org/show_bug.cgi?id=53915
     pa_sink_process_rewind(u->sink, 0); */

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

    dst_chunk->index = 0;
    dst_chunk->length = num_frames * framesize;
    dst_chunk->memblock = pa_memblock_new(sink_input->sink->core->mempool, dst_chunk->length);

    pa_memblockq_drop(u->memblockq, dst_chunk->length);

    src = pa_memblock_acquire_chunk(&tchunk);
    dst = pa_memblock_acquire(dst_chunk->memblock);
    biquad_deinterleave_chunk(src,
                              &deint_chunk,
                              &u->sample_spec,
                              (num_frames * framesize));
    for (chan_idx = 0; chan_idx < u->sample_spec.channels; chan_idx++) {
        current_map_item = &u->filter_map->map[chan_idx];
        pa_biquad_chunk(current_map_item->bqdt1,
                        current_map_item->bqfs1,
                        (&deint_chunk + chan_idx*num_frames),
                        (&reint_chunk + chan_idx*num_frames),
                        num_frames);
    }
    biquad_reinterleave_chunk(&reint_chunk,
                              dst,
                              &u->sample_spec,
                              (num_frames * framesize));

    /* (3) PUT YOUR CODE HERE TO DO SOMETHING WITH THE DATA
    for (frm_idx = 0; frm_idx < num_frames; frm_idx++) {
        cur_frame = src + frm_idx * u->sample_spec.channels;
        dst_frame = dst + frm_idx * u->sample_spec.channels;
        for (chan_idx = 0; chan_idx < u->sample_spec.channels; chan_idx++) {
            cur_sample = cur_frame + chan_idx;
            dst_sample = dst_frame + chan_idx;
            if (u->filter_map[chan_idx] == LOWPASS) {
                // stage 1
                lp = pa_biquad(& (u->s1lpdt[chan_idx]), * (u->s1lpfs), cur_sample);
                bqdtel.w0 = u->s1lpdt[chan_idx].w0;
                bqdtel.y0 = u->s1lpdt[chan_idx].y0;
                pa_store_history(u->s1histbuf, &bqdtel);
                // stage 2
                *dst_sample = pa_biquad(& (u->s2lpdt[chan_idx]), * (u->s2lpfs), &lp);
                bqdtel.w0 = u->s2lpdt[chan_idx].w0;
                bqdtel.y0 = u->s2lpdt[chan_idx].y0;
                pa_store_history(u->s2histbuf, &bqdtel);
            } else if (u->filter_map[chan_idx] == HIGHPASS) {
                // stage 1
                hp = pa_biquad(& (u->s1hpdt[chan_idx]), * (u->s1hpfs), cur_sample);
                bqdtel.w0 = u->s1hpdt[chan_idx].w0;
                bqdtel.y0 = u->s1hpdt[chan_idx].y0;
                pa_store_history(u->s1histbuf, &bqdtel);
                // stage 2
                *dst_sample = pa_biquad(& (u->s2hpdt[chan_idx]), * (u->s2hpfs), &hp);
                bqdtel.w0 = u->s2hpdt[chan_idx].w0;
                bqdtel.y0 = u->s2hpdt[chan_idx].y0;
                pa_store_history(u->s2histbuf, &bqdtel);
            } else if (u->filter_map[chan_idx] == ALLPASS) {
                // stage 1
                ap = pa_biquad(& (u->s1apdt[chan_idx]), * (u->s1apfs), cur_sample);
                bqdtel.w0 = u->s1apdt[chan_idx].w0;
                bqdtel.y0 = u->s1apdt[chan_idx].y0;
                pa_store_history(u->s1histbuf, &bqdtel);
                // stage 2
                *dst_sample = pa_biquad(& (u->s2apdt[chan_idx]), * (u->s2apfs), &ap);
                bqdtel.w0 = u->s2apdt[chan_idx].w0;
                bqdtel.y0 = u->s2apdt[chan_idx].y0;
                pa_store_history(u->s2histbuf, &bqdtel);
            } else {
                pa_log_error("Should never get here, even in Jersey.");
            }
        }
    } */

    pa_memblock_release(tchunk.memblock);
    pa_memblock_release(dst_chunk->memblock);

    pa_memblock_unref(tchunk.memblock);

    return 0;
}

/* Callback function used called from the IO thread context. Causes the sink to
 * rewind the biquad history buffer. */
static void sink_input_process_rewind_cb(pa_sink_input *sink_input, size_t rewind_bytes) {
    struct userdata *u;
    size_t amount = 0;
    size_t rewind_frames = 0;
    size_t rewind_samples = 0;
    size_t max_rewrite = 0;
    size_t frame_size = 0;
    unsigned int i = 0;

    pa_sink_input_assert_ref(sink_input);
    pa_assert_se(u = sink_input->userdata);

    //TODO: call biquad_filter's rewind
/*
    if (u->sink->thread_info.rewind_nbytes > 0) {
        max_rewrite = rewind_bytes + pa_memblockq_get_length(u->memblockq);
        amount = PA_MIN(u->sink->thread_info.rewind_nbytes, max_rewrite);
        u->sink->thread_info.rewind_nbytes = 0;

        if (amount > 0) {
            pa_memblockq_seek(u->memblockq, -(int64_t)amount, PA_SEEK_RELATIVE, TRUE );
             (5) PUT YOUR CODE HERE TO REWIND YOUR FILTER
            rewind_samples = amount / pa_sample_size_of_format(u->sample_spec.format);
            rewind_frames = rewind_samples / u->sample_spec.channels;
            // add the 2 frame backlog for the biquad
            rewind_frames += 2;
            rewind_samples = rewind_frames * u->sample_spec.channels;
            // check if we have to wrap the ring buffer backwards.
            if (u->s1histbuf->idx > rewind_samples) {
                //we're cool, no wrap required
                u->s1histbuf->idx -= rewind_samples;
            } else {
                //doh! gotta wrap.
                u->s1histbuf->idx = u->s1histbuf->length - (rewind_samples - u->s1histbuf->idx);
            }
            frame_size = u->sample_spec.channels * sizeof(biquad_data_element);
            for (i = 0; i < u->sample_spec.channels; i++) {
                switch (u->filter_map[i]) {
                    case LOWPASS: {
                        // stage 1
                        u->s1lpdt[i].w2 = u->s1histbuf->buffer[i].w0;
                        u->s1lpdt[i].w1 = u->s1histbuf->buffer[i + frame_size].w0;
                        u->s1lpdt[i].w0 = u->s1histbuf->buffer[i + 2 * frame_size].w0;
                        u->s1lpdt[i].y2 = u->s1histbuf->buffer[i].y0;
                        u->s1lpdt[i].y1 = u->s1histbuf->buffer[i + frame_size].y0;
                        u->s1lpdt[i].y0 = u->s1histbuf->buffer[i + 2 * frame_size].y0;
                        // stage 2
                        u->s2lpdt[i].w2 = u->s1histbuf->buffer[i].w0;
                        u->s2lpdt[i].w1 = u->s1histbuf->buffer[i + frame_size].w0;
                        u->s2lpdt[i].w0 = u->s1histbuf->buffer[i + 2 * frame_size].w0;
                        u->s2lpdt[i].y2 = u->s1histbuf->buffer[i].y0;
                        u->s2lpdt[i].y1 = u->s1histbuf->buffer[i + frame_size].y0;
                        u->s2lpdt[i].y0 = u->s1histbuf->buffer[i + 2 * frame_size].y0;
                        break;
                    }
                    case HIGHPASS: {
                        // stage 1
                        u->s1hpdt[i].w2 = u->s1histbuf->buffer[i].w0;
                        u->s1hpdt[i].w1 = u->s1histbuf->buffer[i + frame_size].w0;
                        u->s1hpdt[i].w0 = u->s1histbuf->buffer[i + 2 * frame_size].w0;
                        u->s1hpdt[i].y2 = u->s1histbuf->buffer[i].y0;
                        u->s1hpdt[i].y1 = u->s1histbuf->buffer[i + frame_size].y0;
                        u->s1hpdt[i].y0 = u->s1histbuf->buffer[i + 2 * frame_size].y0;
                        // stage 2
                        u->s2hpdt[i].w2 = u->s1histbuf->buffer[i].w0;
                        u->s2hpdt[i].w1 = u->s1histbuf->buffer[i + frame_size].w0;
                        u->s2hpdt[i].w0 = u->s1histbuf->buffer[i + 2 * frame_size].w0;
                        u->s2hpdt[i].y2 = u->s1histbuf->buffer[i].y0;
                        u->s2hpdt[i].y1 = u->s1histbuf->buffer[i + frame_size].y0;
                        u->s2hpdt[i].y0 = u->s1histbuf->buffer[i + 2 * frame_size].y0;
                        break;
                    }
                    case ALLPASS: {
                        // stage 1
                        u->s1apdt[i].w2 = u->s1histbuf->buffer[i].w0;
                        u->s1apdt[i].w1 = u->s1histbuf->buffer[i + frame_size].w0;
                        u->s1apdt[i].w0 = u->s1histbuf->buffer[i + 2 * frame_size].w0;
                        u->s1apdt[i].y2 = u->s1histbuf->buffer[i].y0;
                        u->s1apdt[i].y1 = u->s1histbuf->buffer[i + frame_size].y0;
                        u->s1apdt[i].y0 = u->s1histbuf->buffer[i + 2 * frame_size].y0;
                        // stage 2
                        u->s2apdt[i].w2 = u->s1histbuf->buffer[i].w0;
                        u->s2apdt[i].w1 = u->s1histbuf->buffer[i + frame_size].w0;
                        u->s2apdt[i].w0 = u->s1histbuf->buffer[i + 2 * frame_size].w0;
                        u->s2apdt[i].y2 = u->s1histbuf->buffer[i].y0;
                        u->s2apdt[i].y1 = u->s1histbuf->buffer[i + frame_size].y0;
                        u->s2apdt[i].y0 = u->s1histbuf->buffer[i + 2 * frame_size].y0;
                        break;
                    }
                }
            }
            // move the buffer forward 2 frames so we can write right
            u->s1histbuf->idx += frame_size;
            u->s2histbuf->idx += frame_size;
        }
    }
*/

    pa_sink_process_rewind(u->sink, amount);
    pa_memblockq_rewind(u->memblockq, rewind_bytes);
}

/* FIXME: Too small max_rewind: https://bugs.freedesktop.org/show_bug.cgi?id=53709 */
static void sink_input_update_max_rewind_cb(pa_sink_input *sink_input, size_t max_rewind) {
    struct userdata *u;
    size_t max_rewind_frames = 0;
    pa_sink_input_assert_ref(sink_input);
    pa_assert_se(u = sink_input->userdata);

    if (max_rewind == u->sink->thread_info.max_rewind) {
        pa_log_warn("[%d]%s\n\t called without changing size.\n\tmax_rewind = %lu samples\n",
                    __LINE__, __func__, max_rewind / u->sz_smp);
        pa_memblockq_set_maxrewind(u->memblockq, max_rewind);
        pa_sink_set_max_rewind_within_thread(u->sink, max_rewind);
        return;
    }

    max_rewind_frames = (max_rewind / u->sz_frm);

    if ( max_rewind_frames < MIN_MAX_REWIND_FRAMES) {
        // This is an attempt to prevent excessive realloc shrinks and grows
        pa_log_warn("[%d]%s\n\t called with too small size.\n\tmax_rewind = %lu samples\n",
                    __LINE__, __func__, max_rewind / u->sz_smp);
        pa_memblockq_set_maxrewind(u->memblockq, MIN_MAX_REWIND_FRAMES * u->sz_frm);
        pa_sink_set_max_rewind_within_thread(u->sink, MIN_MAX_REWIND_FRAMES * u->sz_frm);
        return;
    }

    biquad_resize_rewind_buffer(max_rewind_frames, u->filter_map);

    pa_memblockq_set_maxrewind(u->memblockq, max_rewind);
    pa_sink_set_max_rewind_within_thread(u->sink, max_rewind);

}

static void sink_input_update_max_request_cb(pa_sink_input *sink_input, size_t max_request) {
    struct userdata *u;
    pa_sink_input_assert_ref(sink_input);
    pa_assert_se(u = sink_input->userdata);
    pa_sink_set_max_request_within_thread(u->sink, max_request);
}

static void sink_input_update_sink_latency_range_cb(pa_sink_input *sink_input) {
    struct userdata *u;
    pa_sink_input_assert_ref(sink_input);
    pa_assert_se(u = sink_input->userdata);
    pa_sink_set_latency_range_within_thread(u->sink, sink_input->sink->thread_info.min_latency,
                                            sink_input->sink->thread_info.max_latency);
}

static void sink_input_update_sink_fixed_latency_cb(pa_sink_input *i) {
    struct userdata *u;
    pa_sink_input_assert_ref(i);
    pa_assert_se(u = i->userdata);
    pa_sink_set_fixed_latency_within_thread(u->sink, i->sink->thread_info.fixed_latency);
}

static void sink_input_detach_cb(pa_sink_input *sink_input) {
    struct userdata *u;
    pa_sink_input_assert_ref(sink_input);
    pa_assert_se(u = sink_input->userdata);
    pa_sink_detach_within_thread(u->sink);
    pa_sink_set_rtpoll(u->sink, NULL );
}

/* Callback used from IO thread to attach an input to this sink. Importantly,
 * this is where the rewind buffer gets allocated, not in pa__init. */
static void sink_input_attach_cb(pa_sink_input *input) {
    struct userdata *u;
    biquad_map_item_4 *cmi;
    size_t i;
    pa_sink_input_assert_ref(input);
    pa_assert_se(u = input->userdata);

    pa_sink_set_rtpoll(u->sink, input->sink->thread_info.rtpoll);
    pa_sink_set_latency_range_within_thread(u->sink, input->sink->thread_info.min_latency,
                                            input->sink->thread_info.max_latency);
    pa_sink_set_fixed_latency_within_thread(u->sink, input->sink->thread_info.fixed_latency);
    pa_sink_set_max_request_within_thread(u->sink, pa_sink_input_get_max_request(input));

    /* Alloc the rewind buffers */
    //TODO: Move allocation of rewind buffer into biquad-filter.c/h
    for (i = 0; i < u->sample_spec.channels; i++) {
        cmi = &u->filter_map->map[i];
        if ( (u->sink->thread_info.max_rewind / u->sz_frm) < MIN_MAX_REWIND_FRAMES) {
            cmi->bqhs1->length = MIN_MAX_REWIND_FRAMES;
            cmi->bqhs2->length = MIN_MAX_REWIND_FRAMES;
            u->sink->thread_info.max_rewind = MIN_MAX_REWIND_FRAMES * u->sz_frm;
        } else {
            cmi->bqhs1->length = u->sink->thread_info.max_rewind / u->sz_frm;
            cmi->bqhs2->length = u->sink->thread_info.max_rewind / u->sz_frm;
        }
        cmi->bqhs1->buffer = calloc(cmi->bqhs1->length, sizeof(biquad_data_element));
        cmi->bqhs2->buffer = calloc(cmi->bqhs2->length, sizeof(biquad_data_element));
    }

    /* FIXME: Too small max_rewind:
     * https://bugs.freedesktop.org/show_bug.cgi?id=53709 */
    pa_sink_set_max_rewind_within_thread(u->sink, pa_sink_input_get_max_rewind(input));

    pa_sink_attach_within_thread(u->sink);
}

static void sink_input_kill_cb(pa_sink_input *sink_input) {
    struct userdata *u;
    pa_sink_input_assert_ref(sink_input);
    pa_assert_se(u = sink_input->userdata);

    pa_sink_input_unlink(u->sink_input);
    pa_sink_unlink(u->sink);

    pa_sink_input_unref(u->sink_input);
    u->sink_input = NULL;

    pa_sink_unref(u->sink);
    u->sink = NULL;

    pa_module_unload_request(u->module, TRUE );
}

static void sink_input_state_change_cb(pa_sink_input *i, pa_sink_input_state_t state) {
    struct userdata *u;
    pa_sink_input_assert_ref(i);
    pa_assert_se(u = i->userdata);
    /* If we are added for the first time, ask for a rewinding so that we are
     * heard right-away. */
    if (PA_SINK_INPUT_IS_LINKED(state) && i->thread_info.state == PA_SINK_INPUT_INIT) {
        pa_log_debug("Requesting rewind due to state change.");
        pa_sink_input_request_rewind(i, 0, FALSE, TRUE, TRUE );
    }
}

static pa_bool_t sink_input_may_move_to_cb(pa_sink_input *i, pa_sink *dest) {
    struct userdata *u;
    pa_sink_input_assert_ref(i);
    pa_assert_se(u = i->userdata);
    return (u->sink != dest);
}

static void sink_input_moving_cb(pa_sink_input *i, pa_sink *dest) {
    struct userdata *u;

    pa_sink_input_assert_ref(i);
    pa_assert_se(u = i->userdata);

    if (dest) {
        pa_sink_set_asyncmsgq(u->sink, dest->asyncmsgq);
        pa_sink_update_flags(u->sink, PA_SINK_LATENCY | PA_SINK_DYNAMIC_LATENCY, dest->flags);
    } else
        pa_sink_set_asyncmsgq(u->sink, NULL );

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

int pa__init(pa_module *module) {
    struct userdata *u;
    pa_channel_map map;
    biquad_map_item_4 *current_map_item;
    pa_modargs *ma;
    pa_sink *master = NULL;
    pa_sink_input_new_data sink_input_data;
    pa_sink_new_data sink_data;
    pa_memchunk silence;

    pa_assert(module);

    if (! (ma = pa_modargs_new(module->argument, valid_modargs))) {
        pa_log("Failed to parse module arguments.");
        goto fail;
    }

    if (! (master = pa_namereg_get(module->core, pa_modargs_get_value(ma, "master", NULL ), PA_NAMEREG_SINK))) {
        pa_log("Master sink not found");
        goto fail;
    }
    pa_assert(master);

    u = pa_xnew0(struct userdata, 1);

    // get, validate and assign lowpass cutoff freq
    u->lpfreq = atof(pa_modargs_get_value(ma, "lpfreq", "100.0"));
    if (u->lpfreq < MIN_CUTOFF_FREQ) {
        pa_log_error("[%d]%s lpfreq must be between %f and %f.",
                     __LINE__, __func__, MIN_CUTOFF_FREQ, MAX_CUTOFF_FREQ);
        u->lpfreq = MIN_CUTOFF_FREQ;
    }
    if (u->lpfreq > MAX_CUTOFF_FREQ) {
        pa_log_error("[%d]%s lpfreq must be between %f and %f.",
                     __LINE__, __func__, MIN_CUTOFF_FREQ, MAX_CUTOFF_FREQ);
        u->lpfreq = MAX_CUTOFF_FREQ;
    }
    pa_log_info("[%d]%s lpfreq=%f\n", __LINE__, __func__, u->lpfreq);

    u->sample_spec = master->sample_spec;
    u->sample_spec.format = PA_SAMPLE_FLOAT32;

    map = master->channel_map;
/*
    if (pa_modargs_get_sample_spec_and_channel_map(ma, &u->sample_spec, &map, PA_CHANNEL_MAP_DEFAULT) < 0) {
        pa_log_error("Invalid sample format specification or channel map");
        goto fail;
    }
*/

    u->module = module;
    module->userdata = u;

    /* Create sink */
    pa_sink_new_data_init(&sink_data);
    sink_data.driver = __FILE__;
    sink_data.module = module;
    if (! (sink_data.name = pa_xstrdup(pa_modargs_get_value(ma, "sink_name", NULL )))) {
        sink_data.name = pa_sprintf_malloc("%s.lfe_lowpass", master->name);
    }
    pa_sink_new_data_set_sample_spec(&sink_data, &u->sample_spec);
    pa_sink_new_data_set_channel_map(&sink_data, &map);
    pa_proplist_sets(sink_data.proplist, PA_PROP_DEVICE_MASTER_DEVICE, master->name);
    pa_proplist_sets(sink_data.proplist, PA_PROP_DEVICE_CLASS, "filter");
    pa_proplist_sets(sink_data.proplist, "device.name", sink_data.name);

    if (pa_modargs_get_proplist(ma, "sink_properties", sink_data.proplist, PA_UPDATE_REPLACE) < 0) {
        pa_log("Invalid properties");
        pa_sink_new_data_done(&sink_data);
        goto fail;
    }

    if ( (u->auto_desc = !pa_proplist_contains(sink_data.proplist, PA_PROP_DEVICE_DESCRIPTION))) {
        const char *z;
        z = pa_proplist_gets(master->proplist, PA_PROP_DEVICE_DESCRIPTION);
        pa_proplist_setf(sink_data.proplist, PA_PROP_DEVICE_DESCRIPTION, "%s on %s", sink_data.name,
                         z ? z : master->name);
    }

    u->sink = pa_sink_new(
            module->core,
            &sink_data,
            (master->flags & (PA_SINK_LATENCY | PA_SINK_DYNAMIC_LATENCY)) | 0);
    pa_sink_new_data_done(&sink_data);

    if (!u->sink) {
        pa_log("Failed to create sink.");
        goto fail;
    }

    u->sink->parent.process_msg = sink_process_msg_cb;
    u->sink->set_state = sink_set_state_cb;
    u->sink->update_requested_latency = sink_update_requested_latency_cb;
    u->sink->request_rewind = sink_request_rewind_cb;
    u->sink->userdata = u;

    pa_sink_set_asyncmsgq(u->sink, master->asyncmsgq);

    /* Create sink input */
    pa_sink_input_new_data_init(&sink_input_data);
    sink_input_data.driver = __FILE__;
    sink_input_data.module = module;
    pa_sink_input_new_data_set_sink(&sink_input_data, master, FALSE );
    sink_input_data.origin_sink = u->sink;
    pa_proplist_setf(sink_input_data.proplist, PA_PROP_MEDIA_NAME, "lfe-lp Stream from %s",
                     pa_proplist_gets(u->sink->proplist, PA_PROP_DEVICE_DESCRIPTION));
    pa_proplist_sets(sink_input_data.proplist, PA_PROP_MEDIA_ROLE, "filter");
    pa_sink_input_new_data_set_sample_spec(&sink_input_data, &u->sample_spec);
    pa_sink_input_new_data_set_channel_map(&sink_input_data, &map);

    pa_sink_input_new(&u->sink_input, module->core, &sink_input_data);
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
    u->sink_input->volume_changed = NULL;
    u->sink_input->mute_changed = NULL;
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
    // setup shortcuts to common type sizes
    u->sz_smp = pa_sample_size(& (u->sample_spec));
    u->sz_frm = u->sz_smp * u->sample_spec.channels;

    /* If the sink is unlinked, max_rewind doesn't have any meaning anyway. You
     * get notifications about connection state changes via the attach() and
     * detach() callbacks. You can check in attach() what the current max_rewind
     * is and allocate the rewind buffer accordingly.
     * Before the initial attach() call, and between detach() and attach() calls
     * when the sink input is moved, you shouldn't have any need to access the
     * rewind buffer.
     * In addition to updating the rewind buffer in update_max_rewind(), you
     * need to initialize/update it also in attach(), but I didn't think about
     * this before now. <tanuk> */

    /* setup filter_map */
    u->filter_map = malloc(sizeof(biquad_filter_map_4));
    u->filter_map->num_chans = u->sample_spec.channels;
    //u->filter_map->map = malloc(sizeof(biquad_map_item_4) * u->sample_spec.channels);
    for (int i = 0; i < u->sample_spec.channels; i++) {
        current_map_item = &u->filter_map->map[i];
        switch (map.map[i]) {
            case PA_CHANNEL_POSITION_CENTER:
            case PA_CHANNEL_POSITION_REAR_CENTER:
            case PA_CHANNEL_POSITION_FRONT_LEFT_OF_CENTER:
            case PA_CHANNEL_POSITION_FRONT_RIGHT_OF_CENTER:
            case PA_CHANNEL_POSITION_TOP_CENTER:
            case PA_CHANNEL_POSITION_TOP_FRONT_CENTER:
            case PA_CHANNEL_POSITION_TOP_REAR_CENTER:
                current_map_item->type = HIGHPASS;
                current_map_item->bqfs1 = malloc(sizeof(biquad_factors));
                pa_calc_factors(current_map_item->bqfs1, u->sink->sample_spec.rate, u->lpfreq, LOWPASS, 1, 2);
                current_map_item->bqfs2 = malloc(sizeof(biquad_factors));
                pa_calc_factors(current_map_item->bqfs2, u->sink->sample_spec.rate, u->lpfreq, LOWPASS, 2, 2);
                current_map_item->bqdt1 = malloc(sizeof(biquad_data));
                pa_init_bqdt(current_map_item->bqdt1, u->sample_spec.channels);
                current_map_item->bqdt2 = malloc(sizeof(biquad_data));
                pa_init_bqdt(current_map_item->bqdt2, u->sample_spec.channels);
                current_map_item->bqhs1 = malloc(sizeof(biquad_history));
                current_map_item->bqhs1->idx = 0;
                current_map_item->bqhs1->start = 0;
                current_map_item->bqhs1->length = 0;
                current_map_item->bqhs1->buffer = NULL;
                current_map_item->bqhs2 = malloc(sizeof(biquad_history));
                current_map_item->bqhs2->idx = 0;
                current_map_item->bqhs2->start = 0;
                current_map_item->bqhs2->length = 0;
                current_map_item->bqhs2->buffer = NULL;
                break;
            case PA_CHANNEL_POSITION_LFE:
                current_map_item->type = LOWPASS;
                current_map_item->bqfs1 = malloc(sizeof(biquad_factors));
                pa_calc_factors(current_map_item->bqfs1, u->sink->sample_spec.rate, u->lpfreq, HIGHPASS, 1, 2);
                current_map_item->bqfs2 = malloc(sizeof(biquad_factors));
                pa_calc_factors(current_map_item->bqfs2, u->sink->sample_spec.rate, u->lpfreq, HIGHPASS, 2, 2);
                current_map_item->bqdt1 = malloc(sizeof(biquad_data));
                pa_init_bqdt(current_map_item->bqdt1, u->sample_spec.channels);
                current_map_item->bqdt2 = malloc(sizeof(biquad_data));
                pa_init_bqdt(current_map_item->bqdt2, u->sample_spec.channels);
                current_map_item->bqhs1 = malloc(sizeof(biquad_history));
                current_map_item->bqhs1->idx = 0;
                current_map_item->bqhs1->start = 0;
                current_map_item->bqhs1->length = 0;
                current_map_item->bqhs1->buffer = NULL;
                current_map_item->bqhs2 = malloc(sizeof(biquad_history));
                current_map_item->bqhs2->idx = 0;
                current_map_item->bqhs2->start = 0;
                current_map_item->bqhs2->length = 0;
                current_map_item->bqhs2->buffer = NULL;
                break;
            default:
                current_map_item->type = ALLPASS;
                current_map_item->bqfs1 = malloc(sizeof(biquad_factors));
                pa_calc_factors(current_map_item->bqfs1, u->sink->sample_spec.rate, u->lpfreq, ALLPASS, 1, 2);
                current_map_item->bqfs2 = malloc(sizeof(biquad_factors));
                pa_calc_factors(current_map_item->bqfs1, u->sink->sample_spec.rate, u->lpfreq, ALLPASS, 2, 2);
                current_map_item->bqdt1 = malloc(sizeof(biquad_data));
                pa_init_bqdt(current_map_item->bqdt1, u->sample_spec.channels);
                current_map_item->bqdt2 = malloc(sizeof(biquad_data));
                pa_init_bqdt(current_map_item->bqdt2, u->sample_spec.channels);
                current_map_item->bqhs1 = malloc(sizeof(biquad_history));
                current_map_item->bqhs1->idx = 0;
                current_map_item->bqhs1->start = 0;
                current_map_item->bqhs1->length = 0;
                current_map_item->bqhs1->buffer = NULL;
                current_map_item->bqhs2 = malloc(sizeof(biquad_history));
                current_map_item->bqhs2->idx = 0;
                current_map_item->bqhs2->start = 0;
                current_map_item->bqhs2->length = 0;
                current_map_item->bqhs2->buffer = NULL;
        }
    }

    pa_sink_put(u->sink);
    pa_sink_input_put(u->sink_input);
    pa_modargs_free(ma);
    pa_log_debug("Finished lfe-lp.pa__init().\n");
    return 0;

    fail: if (ma)
        pa_modargs_free(ma);
    pa__done(module);
    return -1;
}

int pa__get_n_used(pa_module *m) {
    struct userdata *u;
    pa_assert(m);
    pa_assert_se(u = m->userdata);
    return pa_sink_linked_by(u->sink);
}

/* \note    See comments in sink_input_kill_cb() above about destruction order! */
void pa__done(pa_module*m) {
    struct userdata *u;
    size_t i;
    biquad_map_item_4 *cmi;
    pa_assert(m);

    if (! (u = m->userdata))
        return;

    for (i = 0; i < u->sample_spec.channels; i++) {
        cmi = &u->filter_map->map[i];
        if (cmi->bqdt1)
            free(cmi->bqdt1);
        if (cmi->bqdt2)
            free(cmi->bqdt1);
        if (cmi->bqfs1)
            free(cmi->bqdt1);
        if (cmi->bqfs2)
            free(cmi->bqdt1);
        if (cmi->bqhs1->buffer)
            free(cmi->bqhs1->buffer);
        if (cmi->bqhs2->buffer)
            free(cmi->bqhs2->buffer);
        if (cmi->bqhs1)
            free(cmi->bqhs1);
        if (cmi->bqhs2)
            free(cmi->bqhs2);
    }

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
    pa_log_debug("All done. Bye from module-lfe-lowpass-sink.");
}
