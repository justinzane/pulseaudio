/**
 * \file		biquad-filter.h
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
#define MIN_CUTOFF_FREQ 20.0
#define MAX_CUTOFF_FREQ 500.0

/**
 * \enum biquad_types
 * \note allpass, highpass and lowpass are currently implemented.
 */
typedef enum biquad_types {
    LOWPASS, //!< LOWPASS
    HIGHPASS,//!< HIGHPASS
    ALLPASS  //!< ALLPASS
} biquad_types;

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
 * \struct biquad_data_element: holds one sample of biquad filter history data -- 1/3 of the
 *                              biquad_data struct.  used for the rewind buffer
 */
typedef struct biquad_data_element {
    double y0;
    double w0;
} biquad_data_element;

/**
 * \struct biquad_history  holds the rewind history of filter data
 */
typedef struct biquad_history {
    size_t idx;                      /** < \var idx     write index */
    size_t start;                    /** < \var start   index of oldest element */
    size_t length;                   /** < \var length  buffer length in samples */
    biquad_data_element *buffer;     /** < \var buffer  rewind buffer */
} biquad_history;

/**
 * \fn      filter_biquad
 * \param   [in/out]  bqdt    the filter's working data
 * \param   [in]      bqfs    the coefficients
 * \param   [in]      src     the source sample
 * \return                    the filtered sample
 * \note    y0= (b0 * w0 + b1 * w1 + b2 * w2) − (a1 * y1 + a2 * y2);
 * \see     http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt
 */
float filter_biquad(struct biquad_data *bqdt,
                           struct biquad_factors bqfs,
                           float *src) __attribute__((optimize(3), hot));

/* See below for explanation. */
#define _SQRT_2_2       0.70710678118654757273731092936941422522068023681640625
#define _SQRT_11875_2   0.54486236794258424698256249030237086117267608642578125
#define _SQRT_71875_2   1.3404756618454509720095302327536046504974365234375
const double LINKWITZ_RILEY_Q[4][4] = {{          0.5,     _SQRT_2_2,           0.5, _SQRT_11875_2},
                                       {    _SQRT_2_2,     _SQRT_2_2,           1.0, _SQRT_71875_2},
                                       {          0.5,           1.0,           1.0, _SQRT_11875_2},
                                       {_SQRT_11875_2, _SQRT_71875_2, _SQRT_11875_2, _SQRT_71875_2}};


/**
 * \fn                          filter_calc_factors
 * \param   [in]    sample_rate in Hz
 * \param   [in]    cutoff_freq in Hz, also called corner freq.
 * \param   [in]    type        'a' for allpass,
 *                              'h' for highpass or
 *                              'l' for lowpass
 * \param   [in]    stage        1 for first stage, 2 for second, etc.
 * \param   [in]    num_stages   1 for 2nd order, 2 for 4th order, 3 for 6th order, etc.
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
 * \see http://www.linkwitzlab.com/filters.htm
 *                      LR2     LR4     LR6     LR8
 *      Q0 of stage 1   0.5     0.71    0.5     0.54
 *      Q0 of stage 2           0.71    1.0     1.34
 *      Q0 of stage 3                   1.0     0.54
 *      Q0 of stage 4                           1.34
 *      dB/octave slope 12      24      36      48
 */
void filter_calc_factors(biquad_factors *bqfs,
                                double sample_rate,
                                double cutoff_freq,
                                char type,
                                unsigned int stage,
                                unsigned int num_stages);
/**
 * \fn      filter_init_bqdt
 * \brief   set data elements to 0.0
 * \param [in/out]  bqdt            biquad_data[num_channels]
 * \param [in]      num_channels
 */
void filter_init_bqdt(biquad_data *bqdt,
                             size_t num_channels);

/**
 * \fn      filter_store_history
 * \brief   store the most recent biquad element, [w0, y0], in the history buffer so that we
 *          can rewind without audio inconsistencies by restoring the filter data.
 * \param [in/out]  bqhist  the history buffer
 * \param [in]      bqdtel  the data to be stored
 */
void filter_store_history(biquad_history *bqhist,
                                 biquad_data_element *bqdtel) __attribute__((optimize(3), hot));

#endif /* BIQUAD_FILTER_H_ */
