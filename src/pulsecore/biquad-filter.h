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
//
#include <stdio.h>
#include <math.h>
//#include <time.h>
//#include <float.h>
//#include <xmmintrin.h>
//
//#include <pulse/gccmacro.h>
//#include <pulse/xmalloc.h>
//#include <pulse/def.h>
//
//#include <pulsecore/i18n.h>
//#include <pulsecore/namereg.h>
//#include <pulsecore/sink.h>
//#include <pulsecore/module.h>
#include <pulsecore/core-util.h>
#include <pulsecore/sample-util.h>
//#include <pulsecore/modargs.h>
//#include <pulsecore/log.h>
//#include <pulsecore/rtpoll.h>
//#include <pulsecore/ltdl-helper.h>

#define MEMBLOCKQ_MAXLENGTH (16*1024*1024)
#define MIN_CUTOFF_FREQ 20.0
#define MAX_CUTOFF_FREQ 500.0

/**
 * \enum biquad_types
 * \note allpass, highpass and lowpass are currently implemented.
 */
typedef enum biquad_types {
    LOWPASS,    //!< LOWPASS
    HIGHPASS,   //!< HIGHPASS
    ALLPASS     //!< ALLPASS
} biquad_types;

/**
 * \struct  biquad_filter_map
 * \brief   maps of channels to filter types, i.e.
 *          channels[0]=front-left --> ALLPASS,
 *          channels[5]=lfe        --> LOWPASS, etc.
 */
typedef struct biquad_filter_map {
    biquad_types map[PA_CHANNELS_MAX];
} biquad_filter_map;

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
 * \fn pa_biquad_array
 * \brief filters an entire array. intended to be used with individual channel arrays as
 *        provided by biquad_deinterleave.
 * \param   [in/out]  bqdt          the filter's working data
 * \param   [in]      bqfs          the coefficients
 * \param   [in]      src           the source sample array
 * \param   [out]     dst           the filtered sample array
 * \param   [in]      num_samples   the number of sample in the array
 */
__attribute__((hot)) void pa_biquad_chunk(struct biquad_data *bqdt,
                                          struct biquad_factors bqfs,
                                          float *src,
                                          float *dst,
                                          size_t num_samples);

/**
 * \fn      pa_biquad
 * \param   [in/out]  bqdt    the filter's working data
 * \param   [in]      bqfs    the coefficients
 * \param   [in]      src     the source sample
 * \return                    the filtered sample
 * \note    y0= (b0 * w0 + b1 * w1 + b2 * w2) − (a1 * y1 + a2 * y2);
 * \see     http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt
 */
float pa_biquad(struct biquad_data *bqdt,
                           struct biquad_factors bqfs,
                           float *src) __attribute__((optimize(3), hot));

/**
 * \fn                          pa_calc_factors
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
 */
void pa_calc_factors(biquad_factors *bqfs,
                                double sample_rate,
                                double cutoff_freq,
                                char type,
                                unsigned int stage,
                                unsigned int num_stages);
/**
 * \fn      pa_init_bqdt
 * \brief   set data elements to 0.0
 * \param [in/out]  bqdt            biquad_data[num_channels]
 * \param [in]      num_channels
 */
void pa_init_bqdt(biquad_data *bqdt,
                             size_t num_channels);

/**
 * \fn      pa_store_history
 * \brief   store the most recent biquad element, [w0, y0], in the history buffer so that we
 *          can rewind without audio inconsistencies by restoring the filter data.
 * \param [in/out]  bqhist  the history buffer
 * \param [in]      bqdtel  the data to be stored
 */
void pa_store_history(biquad_history *bqhist,
                      biquad_data_element *bqdtel);

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
                                                    size_t len_chunk);

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
                                                    size_t len_chunk);
#endif /* BIQUAD_FILTER_H_ */
