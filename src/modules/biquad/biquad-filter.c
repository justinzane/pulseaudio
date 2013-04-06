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
#include <biquad/biquad-filter.h>

extern float filter_biquad(struct biquad_data *bqdt, struct biquad_factors bqfs, float *src) {
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


extern void filter_calc_factors(biquad_factors *bqfs, double sample_rate, double cutoff_freq,
                                char type, unsigned int stage, unsigned int num_stages) {
    double w0, alpha;

    /* TODO: validate input args */
    pa_assert(bqfs);
    pa_assert(stage <= num_stages);
    pa_assert(num_stages < 4);
    pa_assert(stage > 0);

    w0 = 2.0 * M_PI * cutoff_freq / sample_rate;
    alpha = sin(w0) / (2.0 * LINKWITZ_RILEY_Q[num_stages-1][stage-1]);

    switch(type) {
        case 'a': {
            (*bqfs).a0 = (1.0 + alpha);
            (*bqfs).a1 = (-2.0 * cos(w0)) / (*bqfs).a0;
            (*bqfs).a2 = (1.0 - alpha) / (*bqfs).a0;
            (*bqfs).b0 = (1.0 - alpha) / (*bqfs).a0;
            (*bqfs).b1 = (-2.0 * cos(w0)) / (*bqfs).a0;
            (*bqfs).b2 = (1.0 + alpha) / (*bqfs).a0;
            (*bqfs).a0 = 1.0;
            break;
        }
        case 'h': {
            (*bqfs).a0 = (1.0 + alpha);
            (*bqfs).a1 = (-2.0 * cos(w0)) / (*bqfs).a0;
            (*bqfs).a2 = (1.0 - alpha) / (*bqfs).a0;
            (*bqfs).b0 = ( (1.0 + cos(w0)) / 2.0) / (*bqfs).a0;
            (*bqfs).b1 = (-1.0 * (1.0 + cos(w0))) / (*bqfs).a0;
            (*bqfs).b2 = ( (1.0 + cos(w0)) / 2.0) / (*bqfs).a0;
            (*bqfs).a0 = 1.0;
            break;
        }
        case 'l': {
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
            pa_log("%d:%s:\n\tInvalid type, %c,  specified. Must be \'a\',\'h\' or \'l\'.",
                   __LINE__, __func__, type);
            return;
        }
    }

    pa_log_info("%d:%s\n", __LINE__, __func__);
    pa_log_info("\tw0=%0.8f, alpha=%0.8f", w0, alpha);
    pa_log_info("\tbiquad_factors\n\t[b0, b1, b2]=[%0.8f, %0.8f, %0.8f]\n\t[a0, a1, a2]=[%0.8f, %0.8f, %0.8f]\n",
                (*bqfs).b0, (*bqfs).b1, (*bqfs).b2, (*bqfs).a0, (*bqfs).a1, (*bqfs).a2);

}

extern void filter_init_bqdt(biquad_data *bqdt, size_t num_channels) {
    size_t i;

    for (i = 0; i < num_channels; i++) {
        bqdt[i].w0 = 0.0; bqdt[i].w1 = 0.0; bqdt[i].w2 = 0.0;
        bqdt[i].y0 = 0.0; bqdt[i].y1 = 0.0; bqdt[i].y2 = 0.0;
    }
}

extern void filter_store_history(biquad_history *bqhist,
                                 biquad_data_element *bqdtel) {
    (*bqhist).buffer[(*bqhist).idx] = *bqdtel;
    (*bqhist).idx += 1;
    if ((*bqhist).idx >= (*bqhist).length)
        (*bqhist).idx = 0;
}

