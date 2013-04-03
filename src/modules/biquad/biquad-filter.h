/**
 * \file		biquad-filter.h
 * \date		Apr 2, 2013
 * \author		justin
 * \copyright	justin
 * \license		GPLv3
    This file is part of pulseaudio.

    pulseaudio is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    pulseaudio is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with pulseaudio.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef BIQUAD_FILTER_H_
#define BIQUAD_FILTER_H_

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

/**
 * \struct biquad_factors: holds biquad filter coefficients/factors for a specific filter type
 * \see    http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt
 */
typedef struct biquad_factors {
    double a0; double a1; double a2;
    double b0; double b1; double b2;
} biquad_factors;

/**
 * \struct biquad_data: holds one iteration of biquad filter history data
 * \see    http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt
 */
typedef struct biquad_data {
    double y0; double y1; double y2;
    double w0; double w1; double w2;
} biquad_data;

/**
 * \struct biquad_data_element: holds one sample of biquad filter history data -- 1/3 of total.
 *                              used for the rewind buffer
 */
typedef struct biquad_data_element {
    double y0;
    double w0;
} biquad_data_element;

/**
 * \struct biquad_history  holds the rewind history of filter data
 */
typedef struct biquad_history {
    size_t idx;
    size_t start;
    size_t length;
    biquad_data_element *buffer;
} biquad_history;

/**
 * \fn      filter_biquad
 * \param   [in/out]  bqdt    the filter's working data
 * \param   [in]      bqfs    the coefficients
 * \param   [in]      src     the source sample
 * \return                    the filtered sample
 * \note    y0= (b0 * x0 + b1 * x1 + b2 * x2) − (a1 * y1 + a2 * y2);
 * \see     http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt
 */
extern float filter_biquad(struct biquad_data *bqdt,
                           struct biquad_factors bqfs,
                           float *src) __attribute__((optimize(3), hot));

/**
 * \fn                          filter_calc_factors
 * \param   [in]    sample_rate in Hz
 * \param   [in]    cutoff_freq in Hz, also called corner freq.
 * \param   [in]    type        'a' for allpass,
 *                              'h' for highpass or
 *                              'l' for lowpass
 * \param   [in]    start        1 for first stage, 2 for second
 * \return  the calculated coefficients/factors
 * \see     http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt
 * \note    Algorithm summary.
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
extern void filter_calc_factors(biquad_factors *bqfs,
                                double sample_rate,
                                double cutoff_freq,
                                char type,
                                unsigned int stage);

/**
 * \fn filter_init_bqdt
 * \brief   Sets the data elements to 0.0.
 */
extern void filter_init_bqdt(biquad_data *bqdt,
                             size_t num_channels);

/**
 * \fn filter_store_history
 * \brief store a biquad_data_element struct in the history buffer so that we
 *        can rewind without audio inconsistencies.
 */
extern void filter_store_history(biquad_history *bqhist,
                                 biquad_data_element *bqdtel) __attribute__((optimize(3), hot));

#endif /* BIQUAD_FILTER_H_ */
