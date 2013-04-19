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

/* ***** Data Definitions ******************************************************************* */

/**
 * \enum biquad_types
 * \note allpass, highpass and lowpass are currently implemented.
 */
typedef enum biquad_types {
    LOWPASS,    //!< LOWPASS
    HIGHPASS,   //!< HIGHPASS
    ALLPASS,    //!< ALLPASS
    LOWSHELF,   //!< LOWSHELF
    HIGHSHELF,  //!< HIGHSHELF  currently unused
    BANDPASS,   //!< BANDPASS   currently unused
    NOTCH,      //!< NOTCH      currently unused
    PEAK,       //!< PEAK       currently unused
} pa_biquad_types;

/**
 * \struct biquad_factors: holds biquad filter coefficients/factors for a specific filter type
 * \see    http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt
 */
typedef struct biquad_factors {
    double a0; double a1; double a2;
    double b0; double b1; double b2;
} pa_biquad_factors_t;

/**
 * \struct biquad_data: holds one iteration of biquad filter history data
 * \see    http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt
 */
typedef struct biquad_data {
    double y0; double y1; double y2;
    double w0; double w1; double w2;
} pa_biquad_data_t;


/**
 * \struct biquad_data_element: holds one sample of biquad filter history data -- 1/3 of the
 *                              biquad_data struct.  used for the rewind buffer
 */
typedef struct biquad_data_element {
    double y0; double w0;
} pa_biquad_data_element_t;

/**
 * \struct biquad_history  holds the rewind history of filter data
 */
typedef struct biquad_history {
    size_t idx;                      /** < \var idx     write index */
    size_t start;                    /** < \var start   index of oldest element */
    size_t length;                   /** < \var length  buffer length in samples */
    pa_biquad_data_element_t *buffer;     /** < \var buffer  rewind buffer */
} pa_biquad_history_t;

/**
 * \struct  biquad_map_item_4
 * \brief   mapping of 4th order filter configuration and data structures to an audio channel
 */
typedef struct biquad_map_item_4 {
    pa_biquad_types    type;         /** < the all/high/low pass filter type*/
    pa_biquad_factors_t *bqfs1;      /** < the stage 1 filter factors*/
    pa_biquad_data_t    *bqdt1;      /** < the stage 1 filter data */
    pa_biquad_factors_t *bqfs2;      /** < the stage 2 filter factors*/
    pa_biquad_data_t    *bqdt2;      /** < the stage 2 filter data */
    pa_biquad_history_t *bqhs1;      /** < the stage 1 rewind history buffer */
    pa_biquad_history_t *bqhs2;      /** < the stage 2 rewind history buffer */
} pa_biquad_map_item_4;

/**
 * \struct  biquad_filter_map_4
 * \brief   mapping between filters and audio channels for a 4th order filter
 */
typedef struct biquad_filter_map_4 {
    size_t            num_chans;            /** < number of channels actually represented */
    pa_biquad_map_item_4 *map;              /** < the data */
} pa_biquad_filter_map_4;

/* ***** Functions ************************************************************************** */

/**
 * \fn pa_biquad_array
 * \brief filters an entire chunk of audio frames
 * \param   [in/out]  fm            pointer to a filter map. the history and working data within
 *                                  the map structure get updated during filtering.
 * \param   [in]      src           pointer to the source sample chunk
 * \param   [out]     dst           pointer to the filtered sample chunk
 * \param   [in]      num_frames    the number of frames in the chunks
 * \param   [in]      num_chans     the number of channels in each frame
 */
__attribute__((hot)) void pa_biquad_chunk_4(struct biquad_filter_map_4 *fm,
                                            float *src,
                                            float *dst,
                                            size_t num_frames,
                                            uint8_t num_chans);

/**
 * \fn      pa_biquad
 * \param   [in/out]  bqdt    pointer to the filter's working data
 * \param   [in]      bqfs    pointer to the coefficients
 * \param   [in]      bqhs    pointer to history
 * \param   [in]      src     pointer to the source sample
 * \param   [in]      dst     pointer to the source sample
 * \note    y0= (b0 * w0 + b1 * w1 + b2 * w2) − (a1 * y1 + a2 * y2);
 * \see     http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt
 */
__attribute__((hot)) void pa_biquad(struct biquad_data *bqdt,
                                    struct biquad_factors *bqfs,
                                    struct biquad_history *bqhs,
                                    float *src,
                                    float *dst);

/**
 * \fn      pa_calc_factors
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
 * LowShelf:   H(s) = A * (s^2 + (sqrt(A)/Q)*s + A)/(A*s^2 + (sqrt(A)/Q)*s + 1)
 *               b0 =    A*( (A+1) - (A-1)*cos(w0) + 2*sqrt(A)*alpha )
 *               b1 =  2*A*( (A-1) - (A+1)*cos(w0)                   )
 *               b2 =    A*( (A+1) - (A-1)*cos(w0) - 2*sqrt(A)*alpha )
 *               a0 =        (A+1) + (A-1)*cos(w0) + 2*sqrt(A)*alpha
 *               a1 =   -2*( (A-1) + (A+1)*cos(w0)                   )
 *               a2 =        (A+1) + (A-1)*cos(w0) - 2*sqrt(A)*alpha
 * \note:   All coefficients are normalized by dividing by the respective a0.
 * \note:   HighShelf, Notch and Bandpass are not implemented.
 */
void pa_biquad_calc_factors(pa_biquad_factors_t *bqfs,  /** < \param  [in/out] bqfs: pointer to the pa_biquad_factors struct being calculated */
                            double sample_rate,         /** < \param  [in] sample_rate: in Hz */
                            double cutoff_freq,         /** < \param  [in] cutoff_freq: in Hz, also called corner freq. */
                            pa_biquad_types type,          /** < \param  [in] type: LOWPASS, ALLPASS, LOWSHELF, etc. */
                            unsigned int stage,         /** < \param  [in] stage: 1 for first stage, 2 for second, etc. */
                            unsigned int num_stages,    /** < \param  [in] num_stages: 1 for 2nd order, 2 for 4th order, etc.*/
                            double gain);               /** < \param  [in] gain: gain, in decibels. not used by LP, HP or AP filters.*/

/**
 * \fn      pa_init_bqdt
 * \brief   allocate and set data elements to 0.0
 * \return  pointer to a biquad_data structure
 */
pa_biquad_data_t *pa_init_bqdt(void);

/**
 * \fn pa_del_bqdt
 * \brief   frees memory allocated to biquad_data structure
 * \param   [in]    bqdt
 */
void pa_del_bqdt(pa_biquad_data_t *bqdt);

/**
 * \fn      pa_init_bqfs
 * \brief   allocate and set data elements to 0.0
 * \param [in/out]  bqfs            biquad_factors
 */
pa_biquad_factors_t *pa_init_bqfs(void);

/**
 * \fn      pa_del_bqfs
 * \brief   frees memory allocated to biquad_factors structure
 * \param   [in]    bqdt
 */
void pa_del_bqfs(pa_biquad_factors_t *bqfs);

/**
 * \fn      pa_init_biquad_history
 * \brief   allocates the index, start and length elements. buffer is NULL.
 * \param   [in/out] bqhs   pointer to the history struct
 */
pa_biquad_history_t *pa_init_biquad_history(void);

/**
 * \fn      pa_del_biquad_history
 * \brief   frees memory allocated to biquad_history structure
 * \param   [in]    bqhs
 */
void pa_del_biquad_history(pa_biquad_history_t *bqhs);

/**
 * \fn pa_init_biquad_filter_map_4
 * \brief   allocates map and sents number of channels
 * \param [in/out]  map         pointer to the str
 * \param [in]      num_chans   from sample_spec.channels
 */
pa_biquad_filter_map_4 *pa_init_biquad_filter_map_4(uint8_t num_chans);

/**
 * \fn      pa_del_filter_map_4
 * \brief   frees memory allocated to biquad_filter_map_4 structure
 * \param   [in]    bqfm
 */
void pa_del_biquad_filter_map_4(pa_biquad_filter_map_4 *bqfm);

/**
 * \fn      biquad_rewind_frames
 * \brief   used by clients like modules to rewind the filters history buffer when the audio
 *          stream is rewound in pulse
 */
__attribute__((hot)) void pa_biquad_rewind_filter(size_t rewind_frames,
                                                  /** < \param [in/out] filter_map: the data structure being rewound */
                                                  pa_biquad_filter_map_4 *filter_map
                                                  /** < \param [in] rewind_frames: the number of frames to rewind.
                                                   * should never be 0 or greater than the length of the buffer */);

/**
 * \fn      biquad_resize_rewind_buffer
 * \brief   logically synonymous with pa_update_max_rewind, resizes the rewind histor buffer
 * \note    it is the author's opinion that the value of rewind_frames should be subject to
 *          minimum and maximum size constraints #defined in pulse/defs.h.
 */
void pa_biquad_resize_rewind_buffer(size_t rewind_frames,
                                    /** < \param [in] rewind_frames: the number of frames to be able to rewind.*/
                                    pa_biquad_filter_map_4 *filter_map
                                    /** < \param [in/out]  filter_map      the data structure being rewound */);
#endif /* BIQUAD_FILTER_H_ */
