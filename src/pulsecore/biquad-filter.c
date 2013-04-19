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

__attribute__((hot)) void pa_biquad_chunk_4(struct biquad_filter_map_4 *fm, float *src,
                                            float *dst, size_t num_frames, uint8_t num_chans) {
    size_t frm_idx, chan_idx;
    float *cur_sample, *dst_sample, *tmp_sample;
    float *cur_frame, *dst_frame;
    tmp_sample = calloc(1,sizeof(float));

    for (frm_idx = 0; frm_idx < num_frames; frm_idx++) {
        cur_frame = src + frm_idx * num_chans;
        dst_frame = dst + frm_idx * num_chans;
        for (chan_idx = 0; chan_idx < num_chans; chan_idx++) {
            cur_sample = cur_frame + chan_idx;
            dst_sample = dst_frame + chan_idx;
            // stage 1
            pa_biquad(fm->map[chan_idx].bqdt1,
                      fm->map[chan_idx].bqfs1,
                      fm->map[chan_idx].bqhs1,
                      cur_sample,
                      tmp_sample);
            // stage 2
            pa_biquad(fm->map[chan_idx].bqdt2,
                      fm->map[chan_idx].bqfs2,
                      fm->map[chan_idx].bqhs2,
                      tmp_sample,
                      dst_sample);
        }
    }
}

__attribute__((hot)) void pa_biquad(struct biquad_data *bqdt,
                                     struct biquad_factors *bqfs,
                                     struct biquad_history *bqhs,
                                     float *src,
                                     float *dst) {
    //#y0= (b0 * x0 + b1 * x1 + b2 * x2) − (a1 * y1 + a2 * y2);
    bqdt->w0 = (double) *src;
    bqdt->y0 = (bqdt->w0 * bqfs->b0 + bqdt->w1 * bqfs->b1 + bqdt->w2 * bqfs->b2) -
               (bqdt->y1 * bqfs->a1 +bqdt->y2 * bqfs->a2);
    bqdt->w2 = bqdt->w1;
    bqdt->w1 = bqdt->w0;
    bqdt->y2 = bqdt->y1;
    bqdt->y1 = bqdt->y0;

    /* Handle rewind history. */
    bqhs->buffer->w0 = bqdt->w0;
    bqhs->buffer->y0 = bqdt->y0;

    bqhs->idx += 1;
    if (bqhs->idx > bqhs->length) {
        bqhs->idx = 0;
    }

    *dst = (float) bqdt->y0;
}

void pa_biquad_calc_factors(pa_biquad_factors_t *bqfs, double sample_rate, double cutoff_freq,
                     pa_biquad_types type, unsigned int stage, unsigned int num_stages, double gain) {
    double w0, alpha, A;

    /* \see http://www.linkwitzlab.com/filters.htm */
    #define _SQRT_2_2       0.70710678118654757273731092936941422522068023681640625
    #define _SQRT_11875_2   0.54486236794258424698256249030237086117267608642578125
    #define _SQRT_71875_2   1.3404756618454509720095302327536046504974365234375
    const double _LINKWITZ_RILEY_Q_1[4] = {          0.5,     _SQRT_2_2,           0.5, _SQRT_11875_2};
    const double _LINKWITZ_RILEY_Q_2[4] = {    _SQRT_2_2,     _SQRT_2_2,           1.0, _SQRT_71875_2};
    const double _LINKWITZ_RILEY_Q_3[4] = {          0.5,           1.0,           1.0, _SQRT_11875_2};
    const double _LINKWITZ_RILEY_Q_4[4] = {_SQRT_11875_2, _SQRT_71875_2, _SQRT_11875_2, _SQRT_71875_2};

    /* TODO: validate input args */
    pa_assert(stage <= num_stages);
    pa_assert(num_stages < 4);
    pa_assert(stage > 0);

    A = pow(10.0,(gain/40.0));
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
        case LOWSHELF: {
            bqfs->a0 = (          ((A + 1.0) + (A - 1.0) * cos(w0) + 2.0 * sqrt(A) * alpha));
            bqfs->a2 = (   -2.0 * ((A + 1.0) + (A - 1.0) * cos(w0)                        )) / bqfs->a0;
            bqfs->a2 = (          ((A + 1.0) + (A - 1.0) * cos(w0) - 2.0 * sqrt(A) * alpha)) / bqfs->a0;
            bqfs->b0 = (      A * ((A + 1.0) - (A - 1.0) * cos(w0) + 2.0 * sqrt(A) * alpha)) / bqfs->a0;
            bqfs->b1 = (2.0 * A * ((A - 1.0) - (A + 1.0) * cos(w0)                        )) / bqfs->a0;
            bqfs->b2 = (      A * ((A + 1.0) - (A - 1.0) * cos(w0) - 2.0 * sqrt(A) * alpha)) / bqfs->a0;
            bqfs->a0 = 1.0;
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

pa_biquad_data_t *pa_init_bqdt() {
    pa_biquad_data_t *bqdt = calloc(1, sizeof(pa_biquad_data_t));

    fprintf(stderr, "JZ %s:%d:%s\n\tbqdt = %p\n", __FILE__, __LINE__, __func__, bqdt);
    bqdt->w0 = 0.0; bqdt->w1 = 0.0; bqdt->w2 = 0.0;
    bqdt->y0 = 0.0; bqdt->y1 = 0.0; bqdt->y2 = 0.0;
    return bqdt;
}

void pa_del_bqdt(pa_biquad_data_t *bqdt) {
    free(bqdt);
}

pa_biquad_factors_t *pa_init_bqfs() {
    pa_biquad_factors_t *bqfs = malloc(sizeof(pa_biquad_factors_t));
    bqfs->a0 = 0.0; bqfs->a1 = 0.0; bqfs->a2 = 0.0;
    bqfs->b0 = 0.0; bqfs->b1 = 0.0; bqfs->b2 = 0.0;
    return bqfs;
}

void pa_del_bqfs(pa_biquad_factors_t *bqfs) {
    free(bqfs);
}

pa_biquad_history_t* pa_init_biquad_history() {
    pa_biquad_history_t *bqhs = calloc(4, sizeof(size_t));
    bqhs->idx = 0;
    bqhs->start = 0;
    bqhs->length = 0;
    bqhs->buffer = NULL;
    return bqhs;
}

void pa_del_biquad_history(pa_biquad_history_t *bqhs) {
    free(bqhs);
}

pa_biquad_filter_map_4 *pa_init_biquad_filter_map_4(uint8_t num_chans) {
    size_t i = 0;
    pa_biquad_filter_map_4 *map = calloc(2, sizeof(size_t));

    map->num_chans = (size_t)num_chans;
    map->map = calloc(num_chans, sizeof(pa_biquad_map_item_4));

    for (i = 0; i < num_chans; i++) {
        map->map[i].bqdt1 = pa_init_bqdt();
        map->map[i].bqdt2 = pa_init_bqdt();
        map->map[i].bqfs1 = pa_init_bqfs();
        map->map[i].bqfs2 = pa_init_bqfs();
        map->map[i].bqhs1 = pa_init_biquad_history();
        map->map[i].bqhs2 = pa_init_biquad_history();
    }
    return map;
}

void pa_del_biquad_filter_map_4(pa_biquad_filter_map_4 *bqfm) {
    size_t i = 0;
    for (i = 0; i < bqfm->num_chans; i++) {
        pa_del_bqdt(bqfm->map[i].bqdt1);
        pa_del_bqdt(bqfm->map[i].bqdt2);
        pa_del_bqfs(bqfm->map[i].bqfs1);
        pa_del_bqfs(bqfm->map[i].bqfs2);
        pa_del_biquad_history(bqfm->map[i].bqhs1);
        pa_del_biquad_history(bqfm->map[i].bqhs2);
    }
    free(bqfm->map);
    free(bqfm);
}

__attribute__((hot)) void pa_biquad_rewind_filter(size_t rewind_frames,
                                               pa_biquad_filter_map_4 *filter_map) {
    size_t i;
    pa_biquad_map_item_4 *cmi;

    if (rewind_frames == 0) {
        pa_log_error("%d : %s : rewind_frames was 0.", __LINE__, __func__);
        return;
    }

    for (i = 0; i < filter_map->num_chans; i++) {
        cmi = &filter_map->map[i];
        if (rewind_frames > cmi->bqhs1->length) {
            pa_log_error("%d : %s : rewind_frames was too long.", __LINE__, __func__);
            break;
        }
        /* Handle stage 1*/
        /* check if we have to wrap the ring buffer backwards. */
        if (cmi->bqhs1->idx > rewind_frames) {  // no wrap needed
            cmi->bqhs1->idx -= rewind_frames;
        } else {                                // must wrap
            cmi->bqhs1->idx = cmi->bqhs1->length - (rewind_frames - cmi->bqhs1->idx);
        }

        /*  pretend loop through frames -2, -1, 0 and put history into working,
            making sure that we do not access a negative index.     */
        switch (cmi->bqhs1->idx) {
            case 0:
                cmi->bqdt1->w2 = (cmi->bqhs1->buffer[cmi->bqhs1->length-2]).w0;
                cmi->bqdt1->w1 = (cmi->bqhs1->buffer[cmi->bqhs1->length-1]).w0;
                cmi->bqdt1->w0 = (cmi->bqhs1->buffer[cmi->bqhs1->idx-0]).w0;
                cmi->bqdt1->y2 = (cmi->bqhs1->buffer[cmi->bqhs1->length-2]).y0;
                cmi->bqdt1->y1 = (cmi->bqhs1->buffer[cmi->bqhs1->length-1]).y0;
                cmi->bqdt1->y0 = (cmi->bqhs1->buffer[cmi->bqhs1->idx-0]).y0;
                break;
            case 1:
                cmi->bqdt1->w2 = (cmi->bqhs1->buffer[cmi->bqhs1->length-1]).w0;
                cmi->bqdt1->w1 = (cmi->bqhs1->buffer[cmi->bqhs1->idx-1]).w0;
                cmi->bqdt1->w0 = (cmi->bqhs1->buffer[cmi->bqhs1->idx-0]).w0;
                cmi->bqdt1->y2 = (cmi->bqhs1->buffer[cmi->bqhs1->length-1]).y0;
                cmi->bqdt1->y1 = (cmi->bqhs1->buffer[cmi->bqhs1->idx-1]).y0;
                cmi->bqdt1->y0 = (cmi->bqhs1->buffer[cmi->bqhs1->idx-0]).y0;
                break;
            default:
                cmi->bqdt1->w2 = (cmi->bqhs1->buffer[cmi->bqhs1->idx-2]).w0;
                cmi->bqdt1->w1 = (cmi->bqhs1->buffer[cmi->bqhs1->idx-1]).w0;
                cmi->bqdt1->w0 = (cmi->bqhs1->buffer[cmi->bqhs1->idx-0]).w0;
                cmi->bqdt1->y2 = (cmi->bqhs1->buffer[cmi->bqhs1->idx-2]).y0;
                cmi->bqdt1->y1 = (cmi->bqhs1->buffer[cmi->bqhs1->idx-1]).y0;
                cmi->bqdt1->y0 = (cmi->bqhs1->buffer[cmi->bqhs1->idx-0]).y0;
        }
        /* Handle stage 2*/
        /* check if we have to wrap the ring buffer backwards. */
        if (cmi->bqhs2->idx > rewind_frames) {  // no wrap needed
            cmi->bqhs2->idx -= rewind_frames;
        } else {                                // must wrap
            cmi->bqhs2->idx = cmi->bqhs2->length - (rewind_frames - cmi->bqhs2->idx);
        }

        /*  pretend loop through frames -2, -1, 0 and put history into working,
            making sure that we do not access a negative index.     */
        switch (cmi->bqhs2->idx) {
            case 0:
                cmi->bqdt2->w2 = (cmi->bqhs2->buffer[cmi->bqhs2->length-2]).w0;
                cmi->bqdt2->w1 = (cmi->bqhs2->buffer[cmi->bqhs2->length-1]).w0;
                cmi->bqdt2->w0 = (cmi->bqhs2->buffer[cmi->bqhs2->idx-0]).w0;
                cmi->bqdt2->y2 = (cmi->bqhs2->buffer[cmi->bqhs2->length-2]).y0;
                cmi->bqdt2->y1 = (cmi->bqhs2->buffer[cmi->bqhs2->length-1]).y0;
                cmi->bqdt2->y0 = (cmi->bqhs2->buffer[cmi->bqhs2->idx-0]).y0;
                break;
            case 1:
                cmi->bqdt2->w2 = (cmi->bqhs2->buffer[cmi->bqhs2->length-1]).w0;
                cmi->bqdt2->w1 = (cmi->bqhs2->buffer[cmi->bqhs2->idx-1]).w0;
                cmi->bqdt2->w0 = (cmi->bqhs2->buffer[cmi->bqhs2->idx-0]).w0;
                cmi->bqdt2->y2 = (cmi->bqhs2->buffer[cmi->bqhs2->length-1]).y0;
                cmi->bqdt2->y1 = (cmi->bqhs2->buffer[cmi->bqhs2->idx-1]).y0;
                cmi->bqdt2->y0 = (cmi->bqhs2->buffer[cmi->bqhs2->idx-0]).y0;
                break;
            default:
                cmi->bqdt2->w2 = (cmi->bqhs2->buffer[cmi->bqhs2->idx-2]).w0;
                cmi->bqdt2->w1 = (cmi->bqhs2->buffer[cmi->bqhs2->idx-1]).w0;
                cmi->bqdt2->w0 = (cmi->bqhs2->buffer[cmi->bqhs2->idx-0]).w0;
                cmi->bqdt2->y2 = (cmi->bqhs2->buffer[cmi->bqhs2->idx-2]).y0;
                cmi->bqdt2->y1 = (cmi->bqhs2->buffer[cmi->bqhs2->idx-1]).y0;
                cmi->bqdt2->y0 = (cmi->bqhs2->buffer[cmi->bqhs2->idx-0]).y0;
        }
    }
}

void pa_biquad_resize_rewind_buffer(size_t rewind_frames,
                                 pa_biquad_filter_map_4 *filter_map) {
    pa_biquad_map_item_4 *cmi;
    size_t shrinkage = 0;
    size_t i = 0;

    fprintf(stderr,"Beginning biquad_resize_rewind_buffer. %lu, %p\n", rewind_frames, filter_map);

    /* We must keep at least 3 samples to repopulate biquad_data */
    rewind_frames += 2;

    for (i = 0; i < filter_map->num_chans; i++) {
        cmi = &filter_map->map[i];
        if (rewind_frames == cmi->bqhs1->length) {
            pa_log_warn("%d:%s called without changing size.", __LINE__, __func__);
            return;
        }

        // grow histbufs to the new size
        if (rewind_frames > cmi->bqhs1->length) {
            pa_log_warn("[%d]%s all in frames\n\tcur rewind_frames = %12lu\n\treq rewind_frames = %12lu\n",
                         __LINE__, __func__, cmi->bqhs1->length, rewind_frames);
            cmi->bqhs1->buffer = realloc(cmi->bqhs1->buffer, rewind_frames*sizeof(pa_biquad_data_element_t));
            cmi->bqhs2->buffer = realloc(cmi->bqhs2->buffer, rewind_frames*sizeof(pa_biquad_data_element_t));
            cmi->bqhs1->length = rewind_frames;
            cmi->bqhs2->length = rewind_frames;
            return;
        }

        // shrink to new size
        if (rewind_frames < cmi->bqhs1->length) {
            pa_log_warn("[%d]%s all in frames\n\tcur rewind_frames = %12lu\n\treq rewind_frames = %12lu\n",
                         __LINE__, __func__, cmi->bqhs1->length, rewind_frames);
            shrinkage = cmi->bqhs1->length - rewind_frames;
            // loop through and move everybody back by shrinkage
            for (i = 0; i < shrinkage; i++)
                cmi->bqhs1->buffer[i + (cmi->bqhs1->length - shrinkage)] = cmi->bqhs1->buffer[i];
            for (i = 0; i < (cmi->bqhs1->length - shrinkage); i++)
                cmi->bqhs1->buffer[i] = cmi->bqhs1->buffer[i + shrinkage];
            cmi->bqhs1->length -= shrinkage;
            cmi->bqhs1->buffer = realloc(cmi->bqhs1->buffer,
                                         (sizeof(pa_biquad_data_element_t) * cmi->bqhs1->length));
            if (cmi->bqhs1->idx > shrinkage)
                cmi->bqhs1->idx -= shrinkage;
            else
                cmi->bqhs1->idx = cmi->bqhs1->length - (shrinkage - cmi->bqhs1->idx);
            return;
        }
    }
}
