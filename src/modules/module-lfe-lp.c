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
#include <time.h>
#include <float.h>
#include <xmmintrin.h>

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

// TODO: Why am I defined here? Shouldn't I come from pulse/pulsecore?
#define MEMBLOCKQ_MAXLENGTH (16*1024*1024)
#ifndef PI
    #ifndef M_PI
        #define PI 3.1415926535897932384626433
    #else
        #define PI M_PI
    #endif
#endif
#define MIN_CUTOFF_FREQ 20.0
#define MAX_CUTOFF_FREQ 500.0
#define MIN_REWIND_FRAMES 1024

/**
 * \struct biquad_factors
 * \brief  holds biquad filter coefficients/factors for a specific filter type
 * \url    http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt
 */
typedef struct biquad_factors {
    double a0;
    double a1;
    double a2;
    double b0;
    double b1;
    double b2;
} biquad_factors;

/**
 * \struct biquad_data
 * \brief  holds one iteration of biquad filter history data
 * \url    http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt
 */
typedef struct biquad_data{
    double y0;
    double y1;
    double y2;
    double w0;
    double w1;
    double w2;
} biquad_data;

/**
 * \struct biquad_data_element
 * \brief  holds one sample of biquad filter history data -- 1/3 of total.
 *         used for the rewind buffer
 */
typedef struct biquad_data_element {
    double y0;
    double w0;
} biquad_data_element;

/**
 * \struct biquad_history
 * \brief  holds the rewind history of filter data
 * \var idx     the current position in the buffer
 * \var start   the initial position in the buffer
 * \var length  the length, in sizeof(biquad_data) of the buffer
 * \var buffer  the buffer
 */
typedef struct biquad_history {
    size_t idx;
    size_t start;
    size_t length;
    biquad_data_element *buffer;
} biquad_history;


PA_MODULE_AUTHOR( _("Justin Chudgar"));
PA_MODULE_DESCRIPTION( _("LFE LP Filter"));
#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION ( _("0.0.0-bogus"))
#endif
PA_MODULE_VERSION(PACKAGE_VERSION);
PA_MODULE_LOAD_ONCE(FALSE);
PA_MODULE_USAGE( _("sink_name=<name for the sink> "
        "sink_properties=<properties for the sink> "
        "master=<name of sink to filter> "
        "lpfreq=low pass cutoff freq 50-500 Hz"
        "use_volume_sharing=<yes or no> "
        "force_flat_volume=<yes or no> " ));

static const char* const valid_modargs[] = { "sink_name", "sink_properties",
                                             "master", "lpfreq",
                                             "use_volume_sharing",
                                             "force_flat_volume", NULL };

/** \var persistent user data structure */
struct userdata {
    pa_module *module;
    pa_bool_t autoloaded;
    pa_sink *sink;
    pa_sink_input *sink_input;
    pa_memblockq *memblockq;
    pa_bool_t auto_desc;
    pa_sample_spec sample_spec;
    /** \var corner/cutoff frequency, user defined */
    double lpfreq;
    /** \var lowpass, highpass and allpass coefficients, respectively */
    struct biquad_factors *lpfs, *hpfs, *apfs;
    /** \var history data for the various filters */
    struct biquad_data *lpdt, *hpdt, *apdt;
    /** \var rewind buffer for biquad_data */
    struct biquad_history *rewind_buf;
    /** \var map of channels index to 'l','h' or 'a' to indicate filter type */
    char filter_map[PA_CHANNELS_MAX];
};

/**
 * \brief do the filtering
 * \param [in]  bqdt    the biquad history data
 * \param [in]  bqfs    the filter factors
 * \param [in]  src     the input sample
 * \return              the filtered sample
 */
static float filter_biquad (struct biquad_data *bqdt, struct biquad_factors bqfs,
                   float *src) {
    //#y0= (b0 * x0 + b1 * x1 + b2 * x2) − (a1 * y1 + a2 * y2);
    (*bqdt).w0 = (double)*src;
    (*bqdt).y0 = (*bqdt).w0 * bqfs.b0 +  (*bqdt).w1 * bqfs.b1 + (*bqdt).w2 * bqfs.b2
                                - ((*bqdt).y1 * bqfs.a1 + (*bqdt).y2 * bqfs.a2);
    (*bqdt).w2 = (*bqdt).w1;
    (*bqdt).w1 = (*bqdt).w0;
    (*bqdt).y2 = (*bqdt).y1;
    (*bqdt).y1 = (*bqdt).y0;

    return((float)(*bqdt).y0);
}

/**
 * \brief   function to calculate filter factors
 * \url     http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt
 * \details
 * LPF:        H(s) = 1 / (s^2 + s/Q + 1)
 *               a0 =   1 + alpha
 *               a1 =  -2*cos(w0)
 *               a2 =   1 - alpha
 *               b0 =  (1 - cos(w0))/2
 *               b1 =  (1 - cos(w0))
 *               b2 =  (1 - cos(w0))/2
 * HPF:        H(s) = s^2 / (s^2 + s/Q + 1)
 *               a0 =   1 + alpha
 *               a1 =  -2*cos(w0)
 *               a2 =   1 - alpha
 *               b0 =  (1 + cos(w0))/2
 *               b1 = -(1 + cos(w0))
 *               b2 =  (1 + cos(w0))/2
 * APF:        H(s) = (s^2 - s/Q + 1) / (s^2 + s/Q + 1)
 *               a0 =   1 + alpha
 *               a1 =  -2*cos(w0)
 *               a2 =   1 - alpha
 *               b0 =   1 - alpha
 *               b1 =  -2*cos(w0)
 *               b2 =   1 + alpha
 * All coefficients are normalized by dividing by the respective a0.
 */
static void filter_calc_factors(pa_sink *sink) {
    struct userdata *u;
    double w0,Q,alpha;
    pa_sink_assert_ref(sink);
    pa_assert_se(u = sink->userdata);

    w0 = 2.0 * PI * u->lpfreq / (double)u->sample_spec.rate;
    Q = sqrt(2.0)/2.0;
    alpha = sin(w0) / (2.0 * Q);

    u->lpfs->a0 = (1.0 + alpha);
    u->lpfs->a1 = (-2.0 * cos(w0))        / u->lpfs->a0;
    u->lpfs->a2 = (1.0 - alpha)           / u->lpfs->a0;
    u->lpfs->b0 = ((1.0 - cos(w0)) / 2.0) / u->lpfs->a0;
    u->lpfs->b1 = ((1.0 - cos(w0))      ) / u->lpfs->a0;
    u->lpfs->b2 = ((1.0 - cos(w0)) / 2.0) / u->lpfs->a0;
    u->lpfs->a0 = 1.0;

    u->hpfs->a0 = (1.0 + alpha);
    u->hpfs->a1 = (-2.0 * cos(w0))           / u->hpfs->a0;
    u->hpfs->a2 = (1.0 - alpha)              / u->hpfs->a0;
    u->hpfs->b0 = (     (1.0 + cos(w0))/2.0) / u->hpfs->a0;
    u->hpfs->b1 = (-1.0*(1.0 + cos(w0))    ) / u->hpfs->a0;
    u->hpfs->b2 = (     (1.0 + cos(w0))/2.0) / u->hpfs->a0;
    u->hpfs->a0 = 1.0;

    u->apfs->a0 = (1.0 + alpha);
    u->apfs->a1 = (-2.0 * cos(w0))  / u->apfs->a0;
    u->apfs->a2 = (1.0 - alpha)     / u->apfs->a0;
    u->apfs->b0 = (1.0 - alpha)     / u->apfs->a0;
    u->apfs->b1 = (-2.0 * cos(w0))  / u->apfs->a0;
    u->apfs->b2 = (1.0 + alpha)     / u->apfs->a0;
    u->apfs->a0 = 1.0;

    pa_log_info ("JZ: %s[%d]\n", __FILE__, __LINE__);
    pa_log_info ("\tw0=%0.8f, Q=%0.8f, alpha=%0.8f", w0, Q, alpha);
    pa_log_info ("\tlowpass_factors\n"
             "\t[b0, b1, b2]=[%0.8f, %0.8f, %0.8f]\n"
             "\t[a0, a1, a2]=[%0.8f, %0.8f, %0.8f]\n",
             u->lpfs->b0, u->lpfs->b1, u->lpfs->b2,
             u->lpfs->a0, u->lpfs->a1, u->lpfs->a2
            );
    pa_log_info ("\thighpass_factors\n"
             "\t[b0, b1, b2]=[%0.8f, %0.8f, %0.8f]\n"
             "\t[a0, a1, a2]=[%0.8f, %0.8f, %0.8f]\n",
             u->hpfs->b0, u->hpfs->b1, u->hpfs->b2,
             u->hpfs->a0, u->hpfs->a1, u->hpfs->a2
            );
    pa_log_info ("\tallpass_factors\n"
             "\t[b0, b1, b2]=[%0.8f, %0.8f, %0.8f]\n"
             "\t[a0, a1, a2]=[%0.8f, %0.8f, %0.8f]\n",
             u->apfs->b0, u->apfs->b1, u->apfs->b2,
             u->apfs->a0, u->apfs->a1, u->apfs->a2
            );
    return;
}

/**
 * \brief initialize filter history data to 0.0
 * \param [in/out]  sink    this sink
 */
static void filter_init_bqdt (pa_sink *sink) {
    struct userdata *u;
    unsigned int i;
    struct biquad_data bqdt;
    pa_sink_assert_ref(sink);
    pa_assert_se(u = sink->userdata);

    for (i=0; i<u->sample_spec.channels; i++) {
        bqdt.w0 = 0.0;
        bqdt.w1 = 0.0;
        bqdt.w2 = 0.0;
        bqdt.y0 = 0.0;
        bqdt.y1 = 0.0;
        bqdt.y2 = 0.0;
        u->lpdt[i] = bqdt;
        u->hpdt[i] = bqdt;
        u->apdt[i] = bqdt;
    }
}

/**
 * \brief store a biquad_data_element struct in the history buffer so that we
 *        can rewind without audio inconsistencies.
 * \param [in/out]  sink    pointer to this sink
 * \param [in]      bqdtel  pointer to the biquad_data_element to be stored
 */
static void filter_store_history (pa_sink *sink,
                                  biquad_data_element *bqdtel) {
    struct userdata *u;
    pa_sink_assert_ref(sink);
    pa_assert_se(u = sink->userdata);

    u->rewind_buf->buffer[u->rewind_buf->idx] = *bqdtel;
    u->rewind_buf->idx += 1;
    if (u->rewind_buf->idx >= u->rewind_buf->length)
        u->rewind_buf->idx = 0;
}

/**
 * \fn sink_process_msg_cb
 * \brief callback function used from the IO thread context
 * @param [in]  msgobject   ???
 * @param [in]  code        ???
 * @param [in]  data        ???
 * @param [in]  offset      ???
 * @param [in]  chunk       ???
 * @return  ???
 * \note Currently only handles PA_SINK_MESSAGE_GET_LATENCY.
 */
static int sink_process_msg_cb(pa_msgobject *msgobject,
                               int code,
                               void *data,
                               int64_t offset,
                               pa_memchunk *chunk) {
    struct userdata *u = PA_SINK(msgobject)->userdata;

    switch (code) {
        case PA_SINK_MESSAGE_GET_LATENCY:
            /* The sink is _put() before the sink input is, so let's make sure
             * we don't access it in that time. Also, the sink input is first
             * shut down, the sink second. */
            if (!PA_SINK_IS_LINKED(u->sink->thread_info.state) ||
                !PA_SINK_INPUT_IS_LINKED(u->sink_input->thread_info.state)) {
                *((pa_usec_t*) data) = 0;
                return(0);
            }
            /* Get the latency of the master sink and add the latency internal
             to our sink input on top */
            *((pa_usec_t*) data) =
                    pa_sink_get_latency_within_thread(u->sink_input->sink) +
                    pa_bytes_to_usec(
                            pa_memblockq_get_length(
                                    u->sink_input->thread_info.render_memblockq),
                            &u->sink_input->sink->sample_spec);
            return(0);
    }
    return(pa_sink_process_msg(msgobject, code, data, offset, chunk));
}

/**
 * \brief callback used from main thread context to set this sink's state
 * @param sink
 * @param state
 * @return always returns 0
 */
//TODO: document me better
static int sink_set_state_cb(pa_sink *sink, pa_sink_state_t state) {
    struct userdata *u;
    pa_sink_assert_ref(sink);
    pa_assert_se(u = sink->userdata);

    if (!PA_SINK_IS_LINKED(state) ||
        !PA_SINK_INPUT_IS_LINKED(pa_sink_input_get_state(u->sink_input)))
        return(0);

    pa_sink_input_cork(u->sink_input, state == PA_SINK_SUSPENDED);
    return(0);
}

/**
 * \brief   callback used from IO thread context to request a rewind
 * \param sink
 */
//TODO: document me better
static void sink_request_rewind_cb(pa_sink *sink) {
    struct userdata *u;
    pa_sink_assert_ref(sink);
    pa_assert_se(u = sink->userdata);

    if (!PA_SINK_IS_LINKED(u->sink->thread_info.state) ||
        !PA_SINK_INPUT_IS_LINKED(u->sink_input->thread_info.state)) {
        pa_log_debug("JZ: %s[%d] sink or sink-input not linked, cannot rewind",
                     __FILE__, __LINE__);
        return;
    }

    /* Just hand this one over to the master sink */
    pa_sink_input_request_rewind(u->sink_input,
                                 (sink->thread_info.rewind_nbytes +
                                  pa_memblockq_get_length(u->memblockq)),
                                 TRUE,
                                 FALSE,
                                 FALSE );
}

/**
 * \brief   called from the IO thread context. ???
 * \param sink  this sink
 */
//TODO: document me better
static void sink_update_requested_latency_cb(pa_sink *sink) {
    struct userdata *u;
    pa_sink_assert_ref(sink);
    pa_assert_se(u = sink->userdata);

    if (!PA_SINK_IS_LINKED(u->sink->thread_info.state) ||
        !PA_SINK_INPUT_IS_LINKED(u->sink_input->thread_info.state))
        return;

    /* Just hand this one over to the master sink */
    pa_sink_input_set_requested_latency_within_thread(
            u->sink_input,
            pa_sink_get_requested_latency_within_thread(sink));
}

/**
 * \brief   callback used from the main thread context to set the sink's volume
 * \param   sink  this sink
 * \note    This seems like a silly thing. A filter should not influence the
 *          general volume of the stream.
 */
//TODO: Remove me safely.
static void sink_set_volume_cb(pa_sink *sink) {
    struct userdata *u;

    pa_sink_assert_ref(sink);
    pa_assert_se(u = sink->userdata);

    if (!PA_SINK_IS_LINKED(pa_sink_get_state(sink)) ||
        !PA_SINK_INPUT_IS_LINKED(pa_sink_input_get_state(u->sink_input)))
        return;

    pa_sink_input_set_volume(u->sink_input,
                             &sink->real_volume,
                             sink->save_volume,
                             TRUE );
}

/**
 * \brief   callback used from the main thread context to mute/unmute the sink
 * \param   sink  this sink
 * \note    This seems like a silly thing. A filter should not influence the
 *          general volume of the stream.
 */
//TODO: Remove me safely or keep me after discussion with reviewers.
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
 * \brief filters the audio data from sink_input
 * \param [in]      sink_input  from whence cometh the noise
 * \param [in]      nbytes      ???
 * \param [in/out]  audio data
 * \return always returns (int)0
 * \note Called from I/O thread context
 */
static int sink_input_pop_cb(pa_sink_input *sink_input,
                             size_t nbytes,
                             pa_memchunk *chunk) {
    struct userdata *u;
    float *src, *dst, *cur_sample, *cur_frame, *dst_frame, *dst_sample;
    float hp = 0.0f;
    float lp = 0.0f;
    float ap = 0.0f;
    size_t framesize;
    unsigned num_frames, frm_idx, chan_idx;
    biquad_data_element bqdtel;
    pa_memchunk tchunk;

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
            if (u->filter_map[chan_idx] == 'l') {
                lp = filter_biquad(&(u->lpdt[chan_idx]), *(u->lpfs), cur_sample);
                *dst_sample = filter_biquad(&(u->lpdt[chan_idx]), *(u->lpfs), &lp);
                bqdtel.w0 = u->lpdt[chan_idx].w0;
                bqdtel.y0 = u->lpdt[chan_idx].y0;
                filter_store_history(u->sink, &bqdtel);
            } else
            if (u->filter_map[chan_idx] == 'h') {
                hp = filter_biquad(&(u->hpdt[chan_idx]), *(u->hpfs), cur_sample);
                *dst_sample = filter_biquad(&(u->hpdt[chan_idx]), *(u->hpfs), &hp);
                bqdtel.w0 = u->hpdt[chan_idx].w0;
                bqdtel.y0 = u->hpdt[chan_idx].y0;
                filter_store_history(u->sink, &bqdtel);
            } else
            if (u->filter_map[chan_idx] == 'a') {
                ap = filter_biquad(&(u->apdt[chan_idx]), *(u->apfs), cur_sample);
                *dst_sample = filter_biquad(&(u->apdt[chan_idx]), *(u->apfs), &ap);
                bqdtel.w0 = u->apdt[chan_idx].w0;
                bqdtel.y0 = u->apdt[chan_idx].y0;
                filter_store_history(u->sink, &bqdtel);
            } else {
                pa_log_error("JZ %s[%d] Should never get here, even in Jersey.",
                             __FILE__, __LINE__);
            }
        }
    }

    pa_memblock_release(tchunk.memblock);
    pa_memblock_release(chunk->memblock);

    pa_memblock_unref(tchunk.memblock);

    return (0);
}

/**
 * \brief Callback function used called from the IO thread context. Causes the
 *        sink to rewind the biquad history buffer.
 * \param [in] sink_input       pointer to the sink input
 * \param [in] rewind_bytes     number of bytes to rewind
 */
static void sink_input_process_rewind_cb(pa_sink_input *sink_input,
                                         size_t rewind_bytes) {
    struct userdata *u;
    size_t amount = 0;
    size_t rewind_frames = 0;
    size_t rewind_samples = 0;
    size_t max_rewrite = 0;
    size_t frame_size = 0;
    unsigned int i = 0;

    pa_sink_input_assert_ref(sink_input);
    pa_assert_se(u = sink_input->userdata);

    if (u->sink->thread_info.rewind_nbytes > 0) {
        max_rewrite = rewind_bytes + pa_memblockq_get_length(u->memblockq);
        amount = PA_MIN(u->sink->thread_info.rewind_nbytes, max_rewrite);
        u->sink->thread_info.rewind_nbytes = 0;

        if (amount > 0) {
            pa_memblockq_seek(u->memblockq, -(int64_t) amount, PA_SEEK_RELATIVE,
                              TRUE );
            /* (5) PUT YOUR CODE HERE TO REWIND YOUR FILTER  */
            rewind_samples = amount / pa_sample_size_of_format(u->sample_spec.format);
            rewind_frames = rewind_samples / u->sample_spec.channels;
            // add the 2 frame backlog for the biquad
            rewind_frames += 2;
            rewind_samples = rewind_frames * u->sample_spec.channels;
            // check if we have to wrap the ring buffer backwards.
            if (u->rewind_buf->idx > rewind_samples) {
                //we're cool, no wrap required
                u->rewind_buf->idx -= rewind_samples;
            } else {
                //doh! gotta wrap.
                u->rewind_buf->idx = u->rewind_buf->length -
                                     (rewind_samples - u->rewind_buf->idx);
            }
            frame_size = u->sample_spec.channels * sizeof(biquad_data_element);
            for (i=0; i < u->sample_spec.channels; i++) {
                switch (u->filter_map[i]) {
                    case 'l': {
                        u->lpdt[i].w2 = u->rewind_buf->buffer[i].w0;
                        u->lpdt[i].w1 = u->rewind_buf->buffer[i+frame_size].w0;
                        u->lpdt[i].w0 = u->rewind_buf->buffer[i+2*frame_size].w0;
                        u->lpdt[i].y2 = u->rewind_buf->buffer[i].y0;
                        u->lpdt[i].y1 = u->rewind_buf->buffer[i+frame_size].y0;
                        u->lpdt[i].y0 = u->rewind_buf->buffer[i+2*frame_size].y0;
                        break;
                    }
                    case 'h': {
                        u->hpdt[i].w2 = u->rewind_buf->buffer[i].w0;
                        u->hpdt[i].w1 = u->rewind_buf->buffer[i+frame_size].w0;
                        u->hpdt[i].w0 = u->rewind_buf->buffer[i+2*frame_size].w0;
                        u->hpdt[i].y2 = u->rewind_buf->buffer[i].y0;
                        u->hpdt[i].y1 = u->rewind_buf->buffer[i+frame_size].y0;
                        u->hpdt[i].y0 = u->rewind_buf->buffer[i+2*frame_size].y0;
                        break;
                    }
                    case 'a': {
                        u->apdt[i].w2 = u->rewind_buf->buffer[i].w0;
                        u->apdt[i].w1 = u->rewind_buf->buffer[i+frame_size].w0;
                        u->apdt[i].w0 = u->rewind_buf->buffer[i+2*frame_size].w0;
                        u->apdt[i].y2 = u->rewind_buf->buffer[i].y0;
                        u->apdt[i].y1 = u->rewind_buf->buffer[i+frame_size].y0;
                        u->apdt[i].y0 = u->rewind_buf->buffer[i+2*frame_size].y0;
                        break;
                    }
                }
            }
            // move the buffer forward 2 frames so we can write right
            u->rewind_buf->idx += frame_size;
        }
    }

    pa_sink_process_rewind(u->sink, amount);
    pa_memblockq_rewind(u->memblockq, rewind_bytes);
}

/**
 * \brief callback function used from IO thread context to update this sink's
 *        input's max_rewind
 * \param [in/out] sink_input   pointer to this sink's input
 * \param [in]     max_rewind   new max_rewind size in bytes
 */
/* FIXME: Too small max_rewind: https://bugs.freedesktop.org/show_bug.cgi?id=53709 */
static void sink_input_update_max_rewind_cb(pa_sink_input *sink_input,
                                            size_t max_rewind) {
    struct userdata *u;
    size_t shrinkage = 0;
    size_t growth = 0;
    size_t i = 0;
    pa_sink_input_assert_ref(sink_input);
    pa_assert_se(u = sink_input->userdata);

    if ((max_rewind / pa_frame_size(&(u->sample_spec))) < MIN_REWIND_FRAMES) {
        if (max_rewind == u->sink->thread_info.max_rewind)
            pa_log("JZ: %s[%d] sink_input_update_max_rewind_cb called "
                   "without changing size.\n\tmax_rewind=%u",
                   __FILE__, __LINE__, max_rewind);
        // This is an attempt to prevent excessive realloc shrinks and grows
        return;
    }

    // grow u->rewind_buf to the new size
    if (max_rewind > u->sink->thread_info.max_rewind) {
        // convert to num samples
        growth = (max_rewind /
                  pa_sample_size_of_format(u->sample_spec.format)) -
                 u->rewind_buf->length;
        pa_log_error("JZ: %s[%d]\n\tmax_rewind=%u, growth=%u, idx=%u, length=%u",
                     __FILE__, __LINE__, max_rewind, growth,
                     u->rewind_buf->idx, u->rewind_buf->length);
        u->rewind_buf = realloc(u->rewind_buf,
                                (sizeof(biquad_data_element) *
                                (u->rewind_buf->length + growth)));
        for (i=u->rewind_buf->length; i<(u->rewind_buf->length+growth); i++) {
            u->rewind_buf->buffer[i].w0 = 0.0;
            u->rewind_buf->buffer[i].y0 = 0.0;
        }
        u->rewind_buf->length += growth;
    }

    // shrink to new size
    if (max_rewind < u->sink->thread_info.max_rewind) {
        // convert to num samples
        shrinkage = u->sink->thread_info.max_rewind -
                    (max_rewind /
                     pa_sample_size_of_format(u->sample_spec.format));
        // loop through and move everybody back by shrinkage
        pa_log_error("JZ: %s[%d]\n\tmax_rewind=%u, shrinkage=%u, idx=%u, length=%u",
                     __FILE__, __LINE__, max_rewind, shrinkage,
                     u->rewind_buf->idx, u->rewind_buf->length);
        for (i=0; i<shrinkage; i++) {
            pa_log_error("\ti=%d",i);
            u->rewind_buf->buffer[i + (u->rewind_buf->length - shrinkage)] =
                    u->rewind_buf->buffer[i];
        }
        for (i=0; i<(u->rewind_buf->length - shrinkage); i++) {
            u->rewind_buf->buffer[i] = u->rewind_buf->buffer[i+shrinkage];
        }
        u->rewind_buf->length -= shrinkage;
        u->rewind_buf = realloc(u->rewind_buf,
                                (sizeof(biquad_data_element) *
                                u->rewind_buf->length));
        if (u->rewind_buf->idx > shrinkage)
            u->rewind_buf->idx -= shrinkage;
        else
            u->rewind_buf->idx = u->rewind_buf->length -
                                 (shrinkage - u->rewind_buf->idx);
    }

    pa_memblockq_set_maxrewind(u->memblockq, max_rewind);
    pa_sink_set_max_rewind_within_thread(u->sink, max_rewind);
}

/**
 * \fn sink_input_update_max_request_cb
 * \brief Callback function used from IO thread context to update this sink's
 *        input's max_request.
 * \param [in/out]  sink_input  this sink's input
 * \param [in]      max_request the new max request size
 */
static void sink_input_update_max_request_cb(pa_sink_input *sink_input,
                                             size_t max_request) {
    struct userdata *u;
    pa_sink_input_assert_ref(sink_input);
    pa_assert_se(u = sink_input->userdata);
    pa_sink_set_max_request_within_thread(u->sink, max_request);
}

/**
 * \brief   callback used from IO thread context to update this sink's latency
 *          with the latency from this sink's input ???
 * \param   sink_input  pointer to this sink's input
 */
//TODO: document me better
static void sink_input_update_sink_latency_range_cb(pa_sink_input *sink_input) {
    struct userdata *u;
    pa_sink_input_assert_ref(sink_input);
    pa_assert_se(u = sink_input->userdata);
    pa_sink_set_latency_range_within_thread(u->sink,
                                            sink_input->sink->thread_info.min_latency,
                                            sink_input->sink->thread_info.max_latency);
}

/**
 * \brief   callback used from IO thread context to update this sink's latency
 *          with the latency from this sink's input ???
 * \param   sink_input  pointer to this sink's input
 */
//TODO: document me better
static void sink_input_update_sink_fixed_latency_cb(pa_sink_input *i) {
    struct userdata *u;
    pa_sink_input_assert_ref(i);
    pa_assert_se(u = i->userdata);
    pa_sink_set_fixed_latency_within_thread(u->sink,
                                            i->sink->thread_info.fixed_latency);
}

/* Called from I/O thread context */
/**
 * \brief   Callback used from IO thread context to detach an input from this
 *          sink.
 * \param   sink_input
 */
//TODO: document me better
static void sink_input_detach_cb(pa_sink_input *sink_input) {
    struct userdata *u;
    pa_sink_input_assert_ref(sink_input);
    pa_assert_se(u = sink_input->userdata);
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

/**
 * \fn sink_input_kill_cb
 * \brief   Callback used by main thread context to remove this sink's input.
 * \details The order here matters! We first kill the sink input, followed by
 *          the sink. That means the sink callbacks must be protected against
 *          an unconnected sink input!
 * \param [in]  sink_input  pointer to the victim, this sink's input
 */
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

    //TODO: document why we unload
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
    return(u->sink != dest);
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

/**
 * \brief   Setup the sink.
 * \param   module  pointer to this module
 * \return  always returns 0
 */
int pa__init(pa_module *module) {
    struct userdata *u;
    pa_channel_map map;
    pa_modargs *ma;
    pa_sink *master = NULL;
    pa_sink_input_new_data sink_input_data;
    pa_sink_new_data sink_data;
    pa_bool_t use_volume_sharing = TRUE;
    pa_bool_t force_flat_volume = FALSE;
    pa_memchunk silence;

    pa_assert(module);

    if (!(ma = pa_modargs_new(module->argument, valid_modargs))) {
        pa_log("Failed to parse module arguments.");
        goto fail;
    }

    if (!(master = pa_namereg_get(module->core,
                                  pa_modargs_get_value(ma, "master", NULL ),
                                  PA_NAMEREG_SINK))) {
        pa_log("Master sink not found");
        goto fail;
    }
    pa_assert(master);

    u = pa_xnew0(struct userdata, 1);

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
    pa_log_info("JZ: %s[%d] lpfreq=%f\n", __FILE__, __LINE__, u->lpfreq);

    u->sample_spec = master->sample_spec;
    u->sample_spec.format = PA_SAMPLE_FLOAT32;

    map = master->channel_map;
    if (pa_modargs_get_sample_spec_and_channel_map(ma,
                                                   &u->sample_spec,
                                                   &map,
                                                   PA_CHANNEL_MAP_DEFAULT) < 0) {
        pa_log("Invalid sample format specification or channel map");
        goto fail;
    }

    /* setup filter_map
     * 'l' for lowpass, 'h' for highpass, 'a' allpass
     * lfe --> lowpass
     * .*center --> high
     * .*left/right/aux -> allpass
     */
    pa_log_info("JZ: %s[%d] filter_map=\n", __FILE__, __LINE__);
    for (int i = 0; i < u->sample_spec.channels; i++) {
        switch (map.map[i]) {
            case PA_CHANNEL_POSITION_CENTER:
            case PA_CHANNEL_POSITION_REAR_CENTER:
            case PA_CHANNEL_POSITION_FRONT_LEFT_OF_CENTER:
            case PA_CHANNEL_POSITION_FRONT_RIGHT_OF_CENTER:
            case PA_CHANNEL_POSITION_TOP_CENTER:
            case PA_CHANNEL_POSITION_TOP_FRONT_CENTER:
            case PA_CHANNEL_POSITION_TOP_REAR_CENTER:
                u->filter_map[i] = 'h';
                break;
            case PA_CHANNEL_POSITION_LFE:
                u->filter_map[i] = 'l';
                break;
            default:
                u->filter_map[i] = 'a';
        }
        pa_log_info("\t%d %c\n", i, u->filter_map[i]);
    }

    if (pa_modargs_get_value_boolean(ma, "use_volume_sharing",
                                     &use_volume_sharing) < 0) {
        pa_log("use_volume_sharing= expects a boolean argument");
        goto fail;
    }

    if (pa_modargs_get_value_boolean(ma, "force_flat_volume",
                                     &force_flat_volume) < 0) {
        pa_log("force_flat_volume= expects a boolean argument");
        goto fail;
    }

    if (use_volume_sharing && force_flat_volume) {
        pa_log("Flat volume can't be forced when using volume sharing.");
        goto fail;
    }

    u->module = module;
    module->userdata = u;

    /* Create sink */
    pa_sink_new_data_init(&sink_data);
    sink_data.driver = __FILE__;
    sink_data.module = module;
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

// here to make Eclipse CDT stop whining
#ifndef PA_SINK_SHARE_VOLUME_WITH_MASTER
#define PA_SINK_SHARE_VOLUME_WITH_MASTER 0x1000000U
#endif
    u->sink = pa_sink_new(module->core,
                          &sink_data,
                          (master->flags &
                           (PA_SINK_LATENCY | PA_SINK_DYNAMIC_LATENCY)) |
                          (use_volume_sharing ? PA_SINK_SHARE_VOLUME_WITH_MASTER : 0));
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
    sink_input_data.module = module;
    pa_sink_input_new_data_set_sink(&sink_input_data, master, FALSE );
    sink_input_data.origin_sink = u->sink;
    pa_proplist_setf(
            sink_input_data.proplist, PA_PROP_MEDIA_NAME,
            "lfe-lp Stream from %s",
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
    // setup filter factors and filter history buffer
    u->lpfs = malloc(sizeof(biquad_factors));
    u->hpfs = malloc(sizeof(biquad_factors));
    u->apfs = malloc(sizeof(biquad_factors));
    u->lpdt = malloc(u->sample_spec.channels * sizeof(biquad_data));
    u->hpdt = malloc(u->sample_spec.channels * sizeof(biquad_data));
    u->apdt = malloc(u->sample_spec.channels * sizeof(biquad_data));
    u->rewind_buf = malloc(4 * sizeof(size_t));
    if ((u->sink->thread_info.max_rewind /
        (pa_sample_size_of_format(u->sample_spec.format) /
         u->sample_spec.channels)) < MIN_REWIND_FRAMES) {
        u->rewind_buf->length = MIN_REWIND_FRAMES * u->sample_spec.channels;
    } else {
        u->rewind_buf->length = u->sink->thread_info.max_rewind /
                                pa_sample_size_of_format(u->sample_spec.format);
    }
    u->rewind_buf->idx = 0;
    u->rewind_buf->start = 0;
    u->rewind_buf->buffer = calloc(u->rewind_buf->length,
                                   sizeof(biquad_data_element));
    // init filter data
    filter_calc_factors(u->sink);
    filter_init_bqdt(u->sink);

    pa_sink_put(u->sink);
    pa_sink_input_put(u->sink_input);
    pa_modargs_free(ma);
    pa_log_debug("JZ: %s[%d] finished lfe-lp.pa__init().", __FILE__, __LINE__);
    return(0);

    fail: if (ma)
        pa_modargs_free(ma);
    pa__done(module);
    return(-1);
}

/**
 * \brief   Return the number of objects linked with this sink.
 * \param   module  pointer to this module
 * \return  number of objects linked with this sink
 */
int pa__get_n_used(pa_module *m) {
    struct userdata *u;
    pa_assert(m);
    pa_assert_se(u = m->userdata);
    return(pa_sink_linked_by(u->sink));
}

/**
 * \fn pa__done
 * \brief   We're done, let's cleanup.
 * \param m
 * \note    See comments in sink_input_kill_cb() above about destruction order!
 */
void pa__done(pa_module*m) {
    struct userdata *u;
    pa_assert(m);

    if (!(u = m->userdata))
        return;

    if (u->lpfs)
        free(u->lpfs);
    if (u->hpfs)
        free(u->hpfs);
    if (u->apfs)
        free(u->apfs);
    if (u->lpdt)
        free(u->lpdt);
    if (u->hpdt)
        free(u->hpdt);
    if (u->apdt)
        free(u->apdt);
    if ((*(u->rewind_buf)).buffer)
        free((*(u->rewind_buf)).buffer);
    if (u->rewind_buf)
        free(u->rewind_buf);

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
    pa_log_debug("JZ: %s[%d] All done. Bye.", __FILE__, __LINE__);
}
