/**
 * \file		biquad-filter.h
 * \date		Apr 2, 2013
 * \author      Justin Chudgar, justin@justinzane.com
 * \copyright   Justin Chudgar
 * \license		LGPLv2.1
>    This file is part of PulseAudio.
>
>    PulseAudio is free software; you can redistribute it and/or modify
>    it under the terms of the GNU Lesser General Public License as published
>    by the Free Software Foundation; either version 2.1 of the License,
>    or (at your option) any later version.
>
>    PulseAudio is distributed in the hope that it will be useful, but
>    WITHOUT ANY WARRANTY; without even the implied warranty of
>    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
>    General Public License for more details.
>
>    You should have received a copy of the GNU Lesser General Public License
>    along with PulseAudio; if not, write to the Free Software
>    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
>    USA.
 */

#ifndef BIQUAD_FILTER_H_
#define BIQUAD_FILTER_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <math.h>
#include <pulsecore/core-util.h>
#include <pulsecore/sample-util.h>

#define MEMBLOCKQ_MAXLENGTH (16*1024*1024)
/**
 * \def     BIQUAD_MIN_LP_FREQ
 * \note    the choice of 20.0 Hz was somewhat arbitrarily chosen by the standard minimum of human
 *          hearing. another factor in this choice is that very few loudspeakers are able to
 *          reproduce sounds below 20.0 Hz effecively.
 */
#define BIQUAD_MIN_LP_FREQ 20.0
/**
 * \def     BIQUAD_MAX_LP_FREQ
 * \note    the choice of 500.0 Hz was very arbitrarily chosen by the author based on the
 *          assumption that it is silly to send anything over 500Hz to a dedicate subwoofer.
 *          if this is a bogus assumption, please let the author know.
 */
#define BIQUAD_MAX_LP_FREQ 500.0
/**
 * \def     BIQUAD_MAX_HP_FREQ
 * \note    the choice of 20.0 KHz was somewhat arbitrarily chosen by the standard minimum of human
 *          hearing. another factor in this choice is that very few loudspeakers are able to
 *          reproduce sounds above 20.0 Hz effecively.
 */
#define BIQUAD_MAX_HP_FREQ 20000.0
/**
 * \def     BIQUAD_MIN_HP_FREQ
 * \note    the choice of 500.0 Hz was very arbitrarily chosen by the author based on the
 *          assumption that it is silly to send anything under 500Hz to a dedicated tweeter.
 *          if this is a bogus assumption, please let the author know.
 */
#define BIQUAD_MIN_HP_FREQ 500.0
#define BIQUAD_MAX_STAGES 4

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
 * \struct biquad_factors
 * \brief  holds biquad filter coefficients/factors for a specific filter type
 * \see    http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt
 */
typedef struct biquad_factors {
    double a0; double a1; double a2;
    double b0; double b1; double b2;
} pa_biquad_factors_t;

/**
 * \struct biquad_data
 * \brief  holds one iteration of biquad filter history data
 * \see    http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt
 */
typedef struct biquad_data {
    double y0; double y1; double y2;
    double w0; double w1; double w2;
} pa_biquad_data_t;


/**
 * \struct biquad_data_element
 * \brief  holds one sample of biquad filter history data -- 1/3 of the
 *         biquad_data struct.  used for the rewind buffer
 */
typedef struct biquad_data_element {
    double y0; double w0;
} pa_biquad_data_element_t;

/**
 * \struct  biquad_history
 * \brief   holds the rewind history of filter data
 */
typedef struct biquad_history {
    size_t idx;                         /**< write index */
    size_t start;                       /**< index of oldest element */
    size_t length;                      /**< buffer length in samples */
    pa_biquad_data_element_t *buffer;   /**< rewind buffer */
} pa_biquad_history_t;

/**
 * \struct  biquad_map_item_4
 * \brief   mapping of 4th order filter configuration and data structures to an audio channel
 */
typedef struct biquad_map_item_4 {
    pa_biquad_types    type;         /**< the all/high/low pass filter type*/
    pa_biquad_factors_t *bqfs1;      /**< the stage 1 filter factors*/
    pa_biquad_data_t    *bqdt1;      /**< the stage 1 filter data */
    pa_biquad_factors_t *bqfs2;      /**< the stage 2 filter factors*/
    pa_biquad_data_t    *bqdt2;      /**< the stage 2 filter data */
    pa_biquad_history_t *bqhs1;      /**< the stage 1 rewind history buffer */
    pa_biquad_history_t *bqhs2;      /**< the stage 2 rewind history buffer */
} pa_biquad_map_item_4;

/**
 * \struct  biquad_filter_map_4
 * \brief   mapping between filters and audio channels for a 4th order filter
 */
typedef struct biquad_filter_map_4 {
    size_t               num_chans;         /**< number of channels actually represented */
    pa_biquad_map_item_4 *map;              /**< the data */
} pa_biquad_filter_map_4;

/* ***** Functions ************************************************************************** */

/**
 * \fn      pa_biquad_chunk_4
 * \brief   filters an entire chunk [pa_memchunk] of audio frames. requires that the source and
 *          destination chunks are floats. working data and history data are updated during
 *          filtering.
 */
__attribute__((hot)) void pa_biquad_chunk_4(struct biquad_filter_map_4 *fm,
                                            /**< [in/out] pointer to a filter map.
                                             * the history and working data within the map
                                             * structure get updated during filtering. */
                                            float *src,
                                            /**< [in] pointer to the source sample chunk */
                                            float *dst,
                                            /**< [out] pointer to the filtered sample chunk */
                                            size_t num_frames,
                                            /**< [in] the number of frames in the chunks */
                                            uint8_t num_chans
                                            /**< [in] the number of channels in each frame */);

/**
 * \fn      pa_biquad
 * \brief   the core filtering function. performs a biquad filter on source data based on
 *          supplied coefficients
 * \note    \f[ y0= (b0 * w0 + b1 * w1 + b2 * w2) − (a1 * y1 + a2 * y2);\f]
 * \see     http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt
 */
__attribute__((hot)) void pa_biquad(struct biquad_data *bqdt,
                                    /**< [in/out] pointer to the filter's working data */
                                    struct biquad_factors *bqfs,
                                    /**< [in] pointer to the coefficients */
                                    struct biquad_history *bqhs,
                                    /**< [in] pointer to history */
                                    float *src,
                                    /**< [in] pointer to the source sample */
                                    float *dst
                                    /**< [in] pointer to the source sample */);

/**
 * \fn      pa_biquad_calc_factors
 * \brief   calculates the factors/coefficients used by the biquad function based on the type
 *          of filter and given parameters
 * \details Algorithm summary.
LPF:        H(s) = 1 / (s^2 + s/Q + 1)
--------------------------------------
              a0 =   1 + alpha
              a1 =  -2*cos(w0)
              a2 =   1 - alpha
              b0 =  (1 - cos(w0))/2
              b1 =  (1 - cos(w0))
              b2 =  (1 - cos(w0))/2
HPF:        H(s) = s^2 / (s^2 + s/Q + 1)
----------------------------------------
              a0 =   1 + alpha
              a1 =  -2*cos(w0)
              a2 =   1 - alpha
              b0 =  (1 + cos(w0))/2
              b1 = -(1 + cos(w0))
              b2 =  (1 + cos(w0))/2
APF:        H(s) = (s^2 - s/Q + 1) / (s^2 + s/Q + 1)
----------------------------------------------------
              a0 =   1 + alpha
              a1 =  -2*cos(w0)
              a2 =   1 - alpha
              b0 =   1 - alpha
              b1 =  -2*cos(w0)
              b2 =   1 + alpha
LowShelf:   H(s) = A * (s^2 + (sqrt(A)/Q)*s + A)/(A*s^2 + (sqrt(A)/Q)*s + 1)
----------------------------------------------------------------------------
              b0 =    A*( (A+1) - (A-1)*cos(w0) + 2*sqrt(A)*alpha )
              b1 =  2*A*( (A-1) - (A+1)*cos(w0)                   )
              b2 =    A*( (A+1) - (A-1)*cos(w0) - 2*sqrt(A)*alpha )
              a0 =        (A+1) + (A-1)*cos(w0) + 2*sqrt(A)*alpha
              a1 =   -2*( (A-1) + (A+1)*cos(w0)                   )
              a2 =        (A+1) + (A-1)*cos(w0) - 2*sqrt(A)*alpha
 * \see     http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt
 * \note:   All coefficients are normalized by dividing by the respective a0.
 * \note:   HighShelf, Notch and Bandpass are not implemented.
 */
void pa_biquad_calc_factors(pa_biquad_factors_t *bqfs,
                            /**< [in/out] pointer to the pa_biquad_factors struct being calculated */
                            double sample_rate,
                            /**< [in]  in Hz */
                            double cutoff_freq,
                            /**< [in] in Hz, also called corner freq. */
                            pa_biquad_types type,
                            /**< [in] LOWPASS, ALLPASS, LOWSHELF, etc. */
                            unsigned int stage,
                            /**< [in] 1 for first stage, 2 for second, etc. */
                            unsigned int num_stages,
                            /**< [in] 1 for 2nd order, 2 for 4th order, etc.*/
                            double gain
                            /**< [in] gain, in dB. not used by LP, HP or AP filters.*/);

/**
 * \fn      pa_init_biquad_data
 * \brief   allocate and set data elements to 0.0
 * \return  pointer to a biquad_data structure
 */
pa_biquad_data_t *pa_init_biquad_data(void);

/**
 * \fn pa_del_biquad_data
 * \brief   frees memory allocated to biquad_data structure
 */
void pa_del_biquad_data(pa_biquad_data_t *bqdt
                        /**< pointer to structure to be freed */);

/**
 * \fn      pa_init_biquad_factors
 * \brief   allocate and set data elements to 0.0
 * \return  pointer to a biquad_factors structure
 */
pa_biquad_factors_t *pa_init_biquad_factors(void);

/**
 * \fn      pa_del_biquad_factors
 * \brief   frees memory allocated to biquad_factors structure
 */
void pa_del_biquad_factors(pa_biquad_factors_t *bqfs
                           /**< pointer to structure to be freed */);

/**
 * \fn      pa_init_biquad_history
 * \brief   allocates the index, start and length elements. buffer is NULL.
 * \return  pointer to a biquad_history struct
 */
pa_biquad_history_t *pa_init_biquad_history(void);

/**
 * \fn      pa_del_biquad_history
 * \brief   frees memory allocated to biquad_history structure
 */
void pa_del_biquad_history(pa_biquad_history_t *bqhs
                           /**< pointer to structure to be freed */);

/**
 * \fn      pa_init_biquad_filter_map_4
 * \brief   allocates map and sents number of channels
 * \returns pointer to biquad_filter_map_4
 */
pa_biquad_filter_map_4 *pa_init_biquad_filter_map_4(uint8_t num_chans
                                                    /**< [in] from sample_spec.channels */);

/**
 * \fn      pa_del_biquad_filter_map_4
 * \brief   frees memory allocated to biquad_filter_map_4 structure
 */
void pa_del_biquad_filter_map_4(pa_biquad_filter_map_4 *bqfm
                                /**< [in] bqfm  pointer to structure to be freed */);

/**
 * \fn      pa_biquad_rewind_filter
 * \brief   used by clients like modules to rewind the filters history buffer when the audio
 *          stream is rewound in pulse
 */
__attribute__((hot)) void pa_biquad_rewind_filter(size_t rewind_frames,
                                                  /**< [in/out] filter_map: the data structure being rewound */
                                                  pa_biquad_filter_map_4 *filter_map
                                                  /**< [in] rewind_frames: the number of frames to rewind.
                                                   * should never be 0 or greater than the length of the buffer */);

/**
 * \fn      pa_biquad_resize_rewind_buffer
 * \brief   logically synonymous with pa_update_max_rewind, resizes the rewind histor buffer
 * \note    it is the author's opinion that the value of rewind_frames should be subject to
 *          minimum and maximum size constraints #defined in pulse/defs.h.
 */
void pa_biquad_resize_rewind_buffer(size_t rewind_frames,
                                    /**< [in] rewind_frames: the number of frames to be able to rewind.*/
                                    pa_biquad_filter_map_4 *filter_map
                                    /**< [in/out]  filter_map      the data structure being rewound */);
#endif /* BIQUAD_FILTER_H_ */
