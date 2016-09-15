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

#ifndef nSGNS_h_
#define nSGNS_h_

#include <algorithm>
#include <cmath>
#include <cstring>

#include "nCore.h"
#include "nMaths.h"
#ifndef _no_fftw_
#include "fftw_interface.h"
#endif // _no_fftw_
#include "nFillFFT.h"
#include "nFFT.h"

//============================================================================
// Maximum order of FIR filter for Adaptive Noise Shaping
//============================================================================
static const int32_t MAX_FILTER_ORDER = 256;

//============================================================================

typedef double  tFFT_Array_Double[MAX_FFT_LENGTH_HALF + 2]     __attribute__ ((aligned(16)));
typedef int32_t     tFFT_Array_Integer[MAX_FFT_LENGTH_HALF + 2]    __attribute__ ((aligned(16)));
typedef bool    tFFT_Array_Boolean[MAX_FFT_LENGTH_HALF + 2]    __attribute__ ((aligned(16)));
typedef double  tFLT_Array[MAX_FILTER_ORDER + 2]              __attribute__ ((aligned(16)));

struct LD_type
{
    tFLT_Array A    __attribute__ ((aligned(16)));
    tFLT_Array B    __attribute__ ((aligned(16)));
    tFLT_Array r    __attribute__ ((aligned(16)));
};

struct Filter_Rec
{
    LD_type LD                  __attribute__ ((aligned(16)));
    tFLT_Array k                __attribute__ ((aligned(16)));
    tFLT_Array State            __attribute__ ((aligned(16)));
    tFLT_Array xState           __attribute__ ((aligned(16)));
    double e;                   __attribute__ ((aligned(16)));
    double One_Over_Minus_C;    __attribute__ ((aligned(16)));
    bool Valid;                 __attribute__ ((aligned(16)));
};

//============================================================================


void nSGNS_Initialise();
void nSGNS_Cleanup();

double Make_Filter(int);

void Warped_Lattice_Filter_Init(int);
double Warped_Lattice_Filter_Evaluate(int);
void Warped_Lattice_Filter_Update(int, double);

bool Filter_Valid(int);
double Filter_Error(int);

#endif // nSGNS_h_
