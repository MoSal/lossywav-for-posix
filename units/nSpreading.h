/**===========================================================================

    lossyWAV: Added noise WAV bit reduction method by David Robinson;
              Noise shaping coefficients by Sebastian Gesemann;

    Copyright (C) 2007-2016 Nick Currie, Copyleft.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: lossywav <at> hotmail <dot> co <dot> uk

==============================================================================
    Initial translation to C++ from Delphi
    Copyright (C) Tyge Løvset (tycho), Aug. 2012
===========================================================================**/

#ifndef nSpreading_h_
#define nSpreading_h_

#include "nCore.h"

static const int32_t SKEWING_AMPLITUDE = 36;
static const int32_t SPREADING_STEPS = 8190;
static const double  SPREADING_STEPS_RECIP = 1.0 / SPREADING_STEPS;

static const int32_t THRESHOLD_INDEX_SPREAD = 64;
static const int32_t THRESHOLD_INDEX_SPREAD_RANGE = 256;
static const int32_t THRESHOLD_INDEX_SPREAD_MAX = THRESHOLD_INDEX_SPREAD_RANGE * THRESHOLD_INDEX_SPREAD;

//============================================================================
// Frequency range over which to check for lowest amplitude signal
//============================================================================
static const int32_t SPREAD_ZONES = 7;

const unsigned char SPREADING_FUNCTION_ARRAY[PRECALC_ANALYSES + 1][SPREAD_ZONES + 1] =
{
    { 1, 1, 1, 1, 1, 1, 1, 1 }, // index 0 not used // TY
    { 1, 2, 2, 2, 2, 2, 2, 2 },
    { 1, 2, 2, 2, 2, 2, 2, 2 },
    { 1, 2, 2, 2, 2, 2, 2, 2 },
    { 1, 2, 2, 2, 2, 2, 2, 3 },
    { 1, 2, 2, 2, 2, 2, 2, 3 },
    { 1, 2, 2, 2, 2, 2, 2, 4 },
    { 1, 1, 2, 2, 2, 2, 3, 4 }
};

static const double SPREADING_FUNCTION_WIDTHS[PRECALC_ANALYSES + 1] =
{
    0.000, // index 0 not used // TY
    sqrt(1./8),
    0.250, // sqrt(1./16),
    sqrt(1./32),
    0.250, // sqrt(1./16),
    sqrt(1./8),
    0.500, // sqrt(1./4),
    sqrt(1./2)
};


struct Spreading_rec
{
    uint16_t width;
    uint16_t fractint;
};


typedef Spreading_rec sprec_array[1];
typedef Spreading_rec* sprec_array_ptr;


// Vars

extern double Frequency_Limits[SPREAD_ZONES + 2]    __attribute__ ((aligned(16)));

extern struct spreading_type
{
    float Widths[SPREADING_STEPS + 2]    __attribute__ ((aligned(16)));
    float divisors[SPREADING_STEPS + 2]    __attribute__ ((aligned(16)));

    double alt_average __attribute__ ((aligned(16)));
    double old_minimum __attribute__ ((aligned(16)));
    double new_minimum __attribute__ ((aligned(16)));

    sprec_array_ptr averages_ptr[PRECALC_ANALYSES + 1]  __attribute__ ((aligned(16))); // base 1.
    double* Bark_Value[PRECALC_ANALYSES + 1]            __attribute__ ((aligned(16))); // base 1.
//============================================================================
    // For each FFT bit length store bits_to_remove with respect to calculated minimum dB of FFT result.
    //============================================================================
    uint16_t threshold_index[THRESHOLD_INDEX_SPREAD_MAX]    __attribute__ ((aligned(16)));

    //============================================================================
    // FFT Bin associated with relevant frequencies for each FFT length.
    //============================================================================
    int32_t Frequency_Bins[PRECALC_ANALYSES + 1][SPREAD_ZONES + 2] __attribute__ ((aligned(16)));

    //============================================================================
    // FFT bins for defining calculation limits.Spreading.
    struct
    {
        int32_t Lower[MAX_FFT_BIT_LENGTH + 1]   __attribute__ ((aligned(16)));
        int32_t Upper[MAX_FFT_BIT_LENGTH + 1]   __attribute__ ((aligned(16)));
        int32_t Mid  [MAX_FFT_BIT_LENGTH + 1]   __attribute__ ((aligned(16)));
        double  Recip[MAX_FFT_BIT_LENGTH + 1]   __attribute__ ((aligned(16)));
    }
    Bins;
}
spreading    __attribute__ ((aligned(16)));

//============================================================================
// skewing function low frequency bin bin attenuation lookup table.
//============================================================================
extern double Skewing_Gain[(MAX_FFT_LENGTH / 2) + 2];

void Spreading_Function_Alt(int);

void Spreading_Function(Results_Type* this_result);

void nSpreading_Init();

void nSpreading_Cleanup();

#endif
