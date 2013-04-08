/**
 * \file		biquad-filter.c
 * \date		Apr 2, 2013
 * \author      Justin Chudgar, justin@justinzane.com
 * \copyright   Justin Chudgar
 * \license		GPLv3
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
#include <pulsecore/biquad-filter.h>

__attribute__((hot)) void pa_biquad_chunk(struct biquad_data *bqdt,
                                          struct biquad_factors bqfs,
                                          float *src,
                                          float *dst,
                                          size_t num_samples) {
    size_t i;
    for (i = 0; i < num_samples; i++) {
        dst[i] = pa_biquad(bqdt, bqfs, &(src[i]));
    }
}

__attribute__((hot)) float pa_biquad(struct biquad_data *bqdt,
                                     struct biquad_factors bqfs,
                                     float *src) {
    //#y0= (b0 * x0 + b1 * x1 + b2 * x2) − (a1 * y1 + a2 * y2);
    (*bqdt).w0 = (double)*src;
    (*bqdt).y0 = ((*bqdt).w0 * bqfs.b0 + (*bqdt).w1 * bqfs.b1 + (*bqdt).w2 * bqfs.b2) -
                 ((*bqdt).y1 * bqfs.a1 +(*bqdt).y2 * bqfs.a2);
    (*bqdt).w2 = (*bqdt).w1;
    (*bqdt).w1 = (*bqdt).w0;
    (*bqdt).y2 = (*bqdt).y1;
    (*bqdt).y1 = (*bqdt).y0;
    return ((float) (*bqdt).y0);
}

void pa_calc_factors(biquad_factors *bqfs, double sample_rate, double cutoff_freq,
                                char type, unsigned int stage, unsigned int num_stages) {
    double w0, alpha;

    /* \see http://www.linkwitzlab.com/filters.htm */
    #define _SQRT_2_2       0.70710678118654757273731092936941422522068023681640625
    #define _SQRT_11875_2   0.54486236794258424698256249030237086117267608642578125
    #define _SQRT_71875_2   1.3404756618454509720095302327536046504974365234375
    const double _LINKWITZ_RILEY_Q_1[4] = {          0.5,     _SQRT_2_2,           0.5, _SQRT_11875_2};
    const double _LINKWITZ_RILEY_Q_2[4] = {    _SQRT_2_2,     _SQRT_2_2,           1.0, _SQRT_71875_2};
    const double _LINKWITZ_RILEY_Q_3[4] = {          0.5,           1.0,           1.0, _SQRT_11875_2};
    const double _LINKWITZ_RILEY_Q_4[4] = {_SQRT_11875_2, _SQRT_71875_2, _SQRT_11875_2, _SQRT_71875_2};

    /* TODO: validate input args */
    pa_assert(bqfs);
    pa_assert(stage <= num_stages);
    pa_assert(num_stages < 4);
    pa_assert(stage > 0);

    w0 = 2.0 * M_PI * cutoff_freq / sample_rate;
    switch (num_stages) {
        case 1:
            alpha = sin(w0) / (2.0 * _LINKWITZ_RILEY_Q_1[stage-1]);
            break;
        case 2:
            alpha = sin(w0) / (2.0 * _LINKWITZ_RILEY_Q_2[stage-1]);
            break;
        case 3:
            alpha = sin(w0) / (2.0 * _LINKWITZ_RILEY_Q_3[stage-1]);
            break;
        case 4:
            alpha = sin(w0) / (2.0 * _LINKWITZ_RILEY_Q_4[stage-1]);
            break;
        default:
            pa_log_error("Invalid value for num_stages.");
    }

    switch(type) {
        case ALLPASS: {
            (*bqfs).a0 = (1.0 + alpha);
            (*bqfs).a1 = (-2.0 * cos(w0)) / (*bqfs).a0;
            (*bqfs).a2 = (1.0 - alpha) / (*bqfs).a0;
            (*bqfs).b0 = (1.0 - alpha) / (*bqfs).a0;
            (*bqfs).b1 = (-2.0 * cos(w0)) / (*bqfs).a0;
            (*bqfs).b2 = (1.0 + alpha) / (*bqfs).a0;
            (*bqfs).a0 = 1.0;
            break;
        }
        case HIGHPASS: {
            (*bqfs).a0 = (1.0 + alpha);
            (*bqfs).a1 = (-2.0 * cos(w0)) / (*bqfs).a0;
            (*bqfs).a2 = (1.0 - alpha) / (*bqfs).a0;
            (*bqfs).b0 = ( (1.0 + cos(w0)) / 2.0) / (*bqfs).a0;
            (*bqfs).b1 = (-1.0 * (1.0 + cos(w0))) / (*bqfs).a0;
            (*bqfs).b2 = ( (1.0 + cos(w0)) / 2.0) / (*bqfs).a0;
            (*bqfs).a0 = 1.0;
            break;
        }
        case LOWPASS: {
            (*bqfs).a0 = (1.0 + alpha);
            (*bqfs).a1 = (-2.0 * cos(w0)) / (*bqfs).a0;
            (*bqfs).a2 = (1.0 - alpha) / (*bqfs).a0;
            (*bqfs).b0 = ( (1.0 - cos(w0)) / 2.0) / (*bqfs).a0;
            (*bqfs).b1 = ( (1.0 - cos(w0))) / (*bqfs).a0;
            (*bqfs).b2 = ( (1.0 - cos(w0)) / 2.0) / (*bqfs).a0;
            (*bqfs).a0 = 1.0;
            break;
        }
        default: {
            pa_log("%d:%s:\n\tInvalid type, %c, specified.", __LINE__, __func__, type);
            return;
        }
    }

    pa_log_info("%d:%s\n", __LINE__, __func__);
    pa_log_info("\tw0=%0.8f, alpha=%0.8f", w0, alpha);
    pa_log_info("\tbiquad_factors\n\t[b0, b1, b2]=[%0.8f, %0.8f, %0.8f]\n\t[a0, a1, a2]=[%0.8f, %0.8f, %0.8f]\n",
                (*bqfs).b0, (*bqfs).b1, (*bqfs).b2, (*bqfs).a0, (*bqfs).a1, (*bqfs).a2);

}

void pa_init_bqdt(biquad_data *bqdt, size_t num_channels) {
    size_t i;

    for (i = 0; i < num_channels; i++) {
        bqdt[i].w0 = 0.0; bqdt[i].w1 = 0.0; bqdt[i].w2 = 0.0;
        bqdt[i].y0 = 0.0; bqdt[i].y1 = 0.0; bqdt[i].y2 = 0.0;
    }
}

void pa_store_history(biquad_history *bqhist,
                                 biquad_data_element *bqdtel) {
    (*bqhist).buffer[(*bqhist).idx] = *bqdtel;
    (*bqhist).idx += 1;
    if ((*bqhist).idx >= (*bqhist).length)
        (*bqhist).idx = 0;
}

/**
 * \fn biquad_deinterleave_chunk
 * \brief Transposes a memchunk from frame oriented to channel oriented; that is, from an
 *        array of frames of channel-indexed samples to an array of channels of frame-indexed
 *        samples.
 * \param [in]  src_chunk   pointer to normal pa_memchunk
 * \param [out] dst_chunk   pointer to deinterleaved pa_memchunk
 * \param [in]  smp_spec    pointer to sample specification
 * \param [in]  len_chunk   length, in bytes, of the memchunks
 */
__attribute__((hot)) void biquad_deinterleave_chunk(pa_memchunk *src_chunk,
                                                    pa_memchunk *dst_chunk,
                                                    pa_sample_spec *smp_spec,
                                                    size_t len_chunk) {
    size_t nf, num_frames, smp_size;
    uint8_t ns, num_chans;

    /* TODO: validate input */
    pa_assert(smp_spec);
    pa_assert(src_chunk);
    pa_assert(dst_chunk);
    num_chans = (*smp_spec).channels;
    smp_size = pa_sample_size(smp_spec);
    pa_assert((len_chunk % (num_chans * smp_size)) == 0);
    num_frames = (len_chunk / num_chans);

    for (nf = 0; nf < num_frames; nf++) {
        for (ns = 0; ns < num_chans; ns++) {
            dst_chunk[(num_frames * ns) + nf] = src_chunk[(num_chans * nf) + ns];
        }
    }
}

/**
 * \fn biquad_reinterleave_chunk
 * \brief Transposes a memchunk from channel oriented to frame oriented; that is, from an
 *        array of channels of frame-indexed samples to an array of frames of channel-indexed
 *        samples.
 * \param [in]  src_chunk   pointer to deinterleaved pa_memchunk
 * \param [out] dst_chunk   pointer to normal pa_memchunk
 * \param [in]  smp_spec    pointer to sample specification
 * \param [in]  len_chunk   length, in bytes, of the memchunks
 */
__attribute__((hot)) void biquad_reinterleave_chunk(pa_memchunk *src_chunk,
                                                    pa_memchunk *dst_chunk,
                                                    pa_sample_spec *smp_spec,
                                                    size_t len_chunk) {
    size_t nf, num_frames, smp_size;
    uint8_t ns, num_chans;

    /* TODO: validate input */
    pa_assert(smp_spec);
    pa_assert(src_chunk);
    pa_assert(dst_chunk);
    num_chans = (*smp_spec).channels;
    smp_size = pa_sample_size(smp_spec);
    pa_assert((len_chunk % (num_chans * smp_size)) == 0);
    num_frames = (len_chunk / num_chans);

    for (nf = 0; nf < num_frames; nf++) {
        for (ns = 0; ns < num_chans; ns++) {
            dst_chunk[(num_frames * ns) + nf] = src_chunk[(num_chans * nf) + ns];
        }
    }
}
