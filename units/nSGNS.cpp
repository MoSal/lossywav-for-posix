//==============================================================================
//
//    lossyWAV: Added noise WAV bit reduction method by David Robinson;
//              Noise shaping coefficients by Sebastian Gesemann;
//
//    Copyright (C) 2007-2013 Nick Currie, Copyleft.
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//    Contact: lossywav <at> hotmail <dot> co <dot> uk
//
//==============================================================================
//    Initial translation to C++ from Delphi
//    by Tyge Løvset (tycho), Aug. 2012
//==============================================================================

#include "nSGNS.h"

namespace { // anonymous

static const int32_t _factor_lookup_shift  __attribute__ ((aligned(16))) = (LONG_FFT_BIT_LENGTH - SHORT_FFT_BIT_LENGTH);
static const int32_t _factor_lookup_mask   __attribute__ ((aligned(16))) = (1 << _factor_lookup_shift) - 1;

//============================================================================
// Frequency and Gain data for creation of fixed noise shaping filter.
//============================================================================
static const float FIXED_GAIN_Freq[37] __attribute__ ((aligned(16))) =
{0.00000, 10.00000, 12.58925, 15.84893, 19.95262, 25.11886, 31.62278, 39.81072, 50.11872, 63.09573, 79.43282,
          100.0000, 125.8925, 158.4893, 199.5262, 251.1886, 316.2278, 398.1072, 501.1872, 630.9573, 794.3282,
          1000.000, 1258.925, 1584.893, 1995.262, 2511.886, 3162.278, 3981.072, 5011.872, 6309.573, 7943.282,
          10000.00, 12589.25, 15848.93, 19952.62, 25118.86,
          1048576};

static const float FNS_Gain[37]  __attribute__ ((aligned(16))) =
{ 15.6000d, 15.5448d, 15.5305d, 15.5126d, 15.4899d, 15.4614d, 15.4255d, 15.3802d, 15.3232d, 15.2515d, 15.1610d, 15.0471d, 14.9034d, 14.7226d, 14.4948d, 14.2076d, 13.8463d, 13.3919d, 12.8226d, 12.1137d, 11.2395d, 10.1635d,  8.8176d,  7.1221d,  4.9766d, 2.1116d, 0.2329d, 1.2827d, 5.4487d, 11.1622d, 16.6788d, 22.2524d, 30.2921d, 60.0341d, 70.0000d, 70.0000d, 70.0000d};

static const float FNS_Spline[37]  __attribute__ ((aligned(16))) =
{  0.0000d,  0.0204d, -0.0018d, -0.0024d, -0.0029d, -0.0037d, -0.0047d, -0.0059d, -0.0073d, -0.0094d, -0.0117d, -0.0149d, -0.0185d, -0.0235d, -0.0297d, -0.0371d, -0.0465d, -0.0574d, -0.0698d, -0.0826d, -0.1009d, -0.1350d, -0.1748d, -0.2250d, -0.3598d, 0.4932d, 1.4643d, 1.5581d, 0.7738d, -0.0985d,  0.0285d,  1.2331d, 10.8512d, -9.8881d, -4.9830d,  0.0000d,  0.0000d};


struct SGNS_type
{
    tFFT_Array_Double  Shape_Curve          __attribute__ ((aligned(16)));
    tFFT_Array_Double  Shape_Total          __attribute__ ((aligned(16)));

    tFFT_Array_Integer Warp_Int             __attribute__ ((aligned(16)));

    tFFT_Array_Double  Warp_Frac            __attribute__ ((aligned(16)));
    tFFT_Array_Double  Warp_1_M_Frac        __attribute__ ((aligned(16)));
    tFFT_Array_Double  Warp_Correction      __attribute__ ((aligned(16)));

    tFFT_Array_Double  Gain                 __attribute__ ((aligned(16)));

    tFFT_Array_Double  Average_Ratio        __attribute__ ((aligned(16)));
    tFFT_Array_Double  Length_Factor        __attribute__ ((aligned(16)));
    tFFT_Array_Double  Length_1_M_Factor    __attribute__ ((aligned(16)));

    double hi_fact[1 << _factor_lookup_shift] __attribute__ ((aligned(16)));
    double lo_fact[1 << _factor_lookup_shift] __attribute__ ((aligned(16)));

    Filter_Rec         Filters[MAX_CHANNELS + 1]  __attribute__ ((aligned(16)));

    double             Lambda               __attribute__ ((aligned(16)));
    double             MinusLambda          __attribute__ ((aligned(16)));
    double             OneMinusLambda2      __attribute__ ((aligned(16)));

    int32_t            Filter_Order         __attribute__ ((aligned(16)));
    int32_t            Filter_Order_M1      __attribute__ ((aligned(16)));

    void*              FFTW_Plan_Inverse    __attribute__ ((aligned(16)));

    FFT_Proc_Rec       nFFT_plan            __attribute__ ((aligned(16)));

    int32_t            Limit_Freq_Bin       __attribute__ ((aligned(16)));
    int32_t            Upper_Freq_Bin       __attribute__ ((aligned(16)));

    bool               Use_Average          __attribute__ ((aligned(16)));

    struct
    {
       void (*Process_Stored_Results)();
       double (*Fill_FFT_With_Warped_Spectrum)();
    } Control;
} SGNS;

} // namespace

#define _filter_rolled_loop_

void Warped_Lattice_Filter_Init(int32_t this_channel)
{
    Filter_Rec* this_Filter = &SGNS.Filters[this_channel];

    double o = 1;
    double u = 1;

    int32_t count=0;

    while (SGNS.Filter_Order-count>63)
    {
        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        double A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;
    }

    if (SGNS.Filter_Order-count>31)
    {
        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        double A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;
    }

    if (SGNS.Filter_Order-count>15)
    {
        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        double A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;
    }

    if (SGNS.Filter_Order-count>7)
    {
        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        double A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;
    }

    if (SGNS.Filter_Order-count>3)
    {
        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        double A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;
    }

    if (SGNS.Filter_Order-count>1)
    {
        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        double A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;

        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;
    }

    if ((SGNS.Filter_Order-count)>0)
    {
        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        double A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;
    }

    this_Filter->One_Over_Minus_C = -1.0d / o;
}


double Warped_Lattice_Filter_Evaluate(int32_t this_channel)
{
    Filter_Rec* this_Filter = &SGNS.Filters[this_channel];

    int32_t count=0;
    double o = 0;
    double u = 0;

    while ((SGNS.Filter_Order - count) > 63)
    {
        double A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        double B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;
    }

    if ((SGNS.Filter_Order-count)>31)
    {
        double A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        double B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;
    }

    if ((SGNS.Filter_Order-count)>15)
    {
        double A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        double B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;
    }

    if ((SGNS.Filter_Order-count)>7)
    {
        double A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        double B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;
    }

    if ((SGNS.Filter_Order-count)>3)
    {
        double A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        double B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;
    }

    if ((SGNS.Filter_Order-count)>1)
    {
        double A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        double B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;

        A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;
    }

    if ((SGNS.Filter_Order-count)>0)
    {
        double A = this_Filter->State[count];
        this_Filter->State[count] = u + SGNS.Lambda * A;
        u = SGNS.OneMinusLambda2 * A - SGNS.Lambda * u;
        double B = o;
        o += this_Filter->k[count] * u;
        u += this_Filter->k[count] * B;
        count++;
    }

    return o * this_Filter->One_Over_Minus_C;
}


void Warped_Lattice_Filter_Update(int32_t this_channel, double fqe)
{
    int32_t count=0;

    Filter_Rec* this_Filter = &SGNS.Filters[this_channel];

    while ((SGNS.Filter_Order-count)>63)
    {
        this_Filter->State[count] += fqe * this_Filter->xState[count];
        this_Filter->State[count + 0x01] += fqe * this_Filter->xState[count + 0x01];
        this_Filter->State[count + 0x02] += fqe * this_Filter->xState[count + 0x02];
        this_Filter->State[count + 0x03] += fqe * this_Filter->xState[count + 0x03];
        this_Filter->State[count + 0x04] += fqe * this_Filter->xState[count + 0x04];
        this_Filter->State[count + 0x05] += fqe * this_Filter->xState[count + 0x05];
        this_Filter->State[count + 0x06] += fqe * this_Filter->xState[count + 0x06];
        this_Filter->State[count + 0x07] += fqe * this_Filter->xState[count + 0x07];
        this_Filter->State[count + 0x08] += fqe * this_Filter->xState[count + 0x08];
        this_Filter->State[count + 0x09] += fqe * this_Filter->xState[count + 0x09];
        this_Filter->State[count + 0x0a] += fqe * this_Filter->xState[count + 0x0a];
        this_Filter->State[count + 0x0b] += fqe * this_Filter->xState[count + 0x0b];
        this_Filter->State[count + 0x0c] += fqe * this_Filter->xState[count + 0x0c];
        this_Filter->State[count + 0x0d] += fqe * this_Filter->xState[count + 0x0d];
        this_Filter->State[count + 0x0e] += fqe * this_Filter->xState[count + 0x0e];
        this_Filter->State[count + 0x0f] += fqe * this_Filter->xState[count + 0x0f];
        this_Filter->State[count + 0x10] += fqe * this_Filter->xState[count + 0x10];
        this_Filter->State[count + 0x11] += fqe * this_Filter->xState[count + 0x11];
        this_Filter->State[count + 0x12] += fqe * this_Filter->xState[count + 0x12];
        this_Filter->State[count + 0x13] += fqe * this_Filter->xState[count + 0x13];
        this_Filter->State[count + 0x14] += fqe * this_Filter->xState[count + 0x14];
        this_Filter->State[count + 0x15] += fqe * this_Filter->xState[count + 0x15];
        this_Filter->State[count + 0x16] += fqe * this_Filter->xState[count + 0x16];
        this_Filter->State[count + 0x17] += fqe * this_Filter->xState[count + 0x17];
        this_Filter->State[count + 0x18] += fqe * this_Filter->xState[count + 0x18];
        this_Filter->State[count + 0x19] += fqe * this_Filter->xState[count + 0x19];
        this_Filter->State[count + 0x1a] += fqe * this_Filter->xState[count + 0x1a];
        this_Filter->State[count + 0x1b] += fqe * this_Filter->xState[count + 0x1b];
        this_Filter->State[count + 0x1c] += fqe * this_Filter->xState[count + 0x1c];
        this_Filter->State[count + 0x1d] += fqe * this_Filter->xState[count + 0x1d];
        this_Filter->State[count + 0x1e] += fqe * this_Filter->xState[count + 0x1e];
        this_Filter->State[count + 0x1f] += fqe * this_Filter->xState[count + 0x1f];
        this_Filter->State[count + 0x20] += fqe * this_Filter->xState[count + 0x20];
        this_Filter->State[count + 0x21] += fqe * this_Filter->xState[count + 0x21];
        this_Filter->State[count + 0x22] += fqe * this_Filter->xState[count + 0x22];
        this_Filter->State[count + 0x23] += fqe * this_Filter->xState[count + 0x23];
        this_Filter->State[count + 0x24] += fqe * this_Filter->xState[count + 0x24];
        this_Filter->State[count + 0x25] += fqe * this_Filter->xState[count + 0x25];
        this_Filter->State[count + 0x26] += fqe * this_Filter->xState[count + 0x26];
        this_Filter->State[count + 0x27] += fqe * this_Filter->xState[count + 0x27];
        this_Filter->State[count + 0x28] += fqe * this_Filter->xState[count + 0x28];
        this_Filter->State[count + 0x29] += fqe * this_Filter->xState[count + 0x29];
        this_Filter->State[count + 0x2a] += fqe * this_Filter->xState[count + 0x2a];
        this_Filter->State[count + 0x2b] += fqe * this_Filter->xState[count + 0x2b];
        this_Filter->State[count + 0x2c] += fqe * this_Filter->xState[count + 0x2c];
        this_Filter->State[count + 0x2d] += fqe * this_Filter->xState[count + 0x2d];
        this_Filter->State[count + 0x2e] += fqe * this_Filter->xState[count + 0x2e];
        this_Filter->State[count + 0x2f] += fqe * this_Filter->xState[count + 0x2f];
        this_Filter->State[count + 0x30] += fqe * this_Filter->xState[count + 0x30];
        this_Filter->State[count + 0x31] += fqe * this_Filter->xState[count + 0x31];
        this_Filter->State[count + 0x32] += fqe * this_Filter->xState[count + 0x32];
        this_Filter->State[count + 0x33] += fqe * this_Filter->xState[count + 0x33];
        this_Filter->State[count + 0x34] += fqe * this_Filter->xState[count + 0x34];
        this_Filter->State[count + 0x35] += fqe * this_Filter->xState[count + 0x35];
        this_Filter->State[count + 0x36] += fqe * this_Filter->xState[count + 0x36];
        this_Filter->State[count + 0x37] += fqe * this_Filter->xState[count + 0x37];
        this_Filter->State[count + 0x38] += fqe * this_Filter->xState[count + 0x38];
        this_Filter->State[count + 0x39] += fqe * this_Filter->xState[count + 0x39];
        this_Filter->State[count + 0x3a] += fqe * this_Filter->xState[count + 0x3a];
        this_Filter->State[count + 0x3b] += fqe * this_Filter->xState[count + 0x3b];
        this_Filter->State[count + 0x3c] += fqe * this_Filter->xState[count + 0x3c];
        this_Filter->State[count + 0x3d] += fqe * this_Filter->xState[count + 0x3d];
        this_Filter->State[count + 0x3e] += fqe * this_Filter->xState[count + 0x3e];
        this_Filter->State[count + 0x3f] += fqe * this_Filter->xState[count + 0x3f];
        count += 64;
    }

    if ((SGNS.Filter_Order-count)>31)
    {
        this_Filter->State[count] += fqe * this_Filter->xState[count];
        this_Filter->State[count + 0x01] += fqe * this_Filter->xState[count + 0x01];
        this_Filter->State[count + 0x02] += fqe * this_Filter->xState[count + 0x02];
        this_Filter->State[count + 0x03] += fqe * this_Filter->xState[count + 0x03];
        this_Filter->State[count + 0x04] += fqe * this_Filter->xState[count + 0x04];
        this_Filter->State[count + 0x05] += fqe * this_Filter->xState[count + 0x05];
        this_Filter->State[count + 0x06] += fqe * this_Filter->xState[count + 0x06];
        this_Filter->State[count + 0x07] += fqe * this_Filter->xState[count + 0x07];
        this_Filter->State[count + 0x08] += fqe * this_Filter->xState[count + 0x08];
        this_Filter->State[count + 0x09] += fqe * this_Filter->xState[count + 0x09];
        this_Filter->State[count + 0x0a] += fqe * this_Filter->xState[count + 0x0a];
        this_Filter->State[count + 0x0b] += fqe * this_Filter->xState[count + 0x0b];
        this_Filter->State[count + 0x0c] += fqe * this_Filter->xState[count + 0x0c];
        this_Filter->State[count + 0x0d] += fqe * this_Filter->xState[count + 0x0d];
        this_Filter->State[count + 0x0e] += fqe * this_Filter->xState[count + 0x0e];
        this_Filter->State[count + 0x0f] += fqe * this_Filter->xState[count + 0x0f];
        this_Filter->State[count + 0x10] += fqe * this_Filter->xState[count + 0x10];
        this_Filter->State[count + 0x11] += fqe * this_Filter->xState[count + 0x11];
        this_Filter->State[count + 0x12] += fqe * this_Filter->xState[count + 0x12];
        this_Filter->State[count + 0x13] += fqe * this_Filter->xState[count + 0x13];
        this_Filter->State[count + 0x14] += fqe * this_Filter->xState[count + 0x14];
        this_Filter->State[count + 0x15] += fqe * this_Filter->xState[count + 0x15];
        this_Filter->State[count + 0x16] += fqe * this_Filter->xState[count + 0x16];
        this_Filter->State[count + 0x17] += fqe * this_Filter->xState[count + 0x17];
        this_Filter->State[count + 0x18] += fqe * this_Filter->xState[count + 0x18];
        this_Filter->State[count + 0x19] += fqe * this_Filter->xState[count + 0x19];
        this_Filter->State[count + 0x1a] += fqe * this_Filter->xState[count + 0x1a];
        this_Filter->State[count + 0x1b] += fqe * this_Filter->xState[count + 0x1b];
        this_Filter->State[count + 0x1c] += fqe * this_Filter->xState[count + 0x1c];
        this_Filter->State[count + 0x1d] += fqe * this_Filter->xState[count + 0x1d];
        this_Filter->State[count + 0x1e] += fqe * this_Filter->xState[count + 0x1e];
        this_Filter->State[count + 0x1f] += fqe * this_Filter->xState[count + 0x1f];
        count += 32;
    }

    if ((SGNS.Filter_Order-count)>15)
    {
        this_Filter->State[count] += fqe * this_Filter->xState[count];
        this_Filter->State[count +  1] += fqe * this_Filter->xState[count +  1];
        this_Filter->State[count +  2] += fqe * this_Filter->xState[count +  2];
        this_Filter->State[count +  3] += fqe * this_Filter->xState[count +  3];
        this_Filter->State[count +  4] += fqe * this_Filter->xState[count +  4];
        this_Filter->State[count +  5] += fqe * this_Filter->xState[count +  5];
        this_Filter->State[count +  6] += fqe * this_Filter->xState[count +  6];
        this_Filter->State[count +  7] += fqe * this_Filter->xState[count +  7];
        this_Filter->State[count +  8] += fqe * this_Filter->xState[count +  8];
        this_Filter->State[count +  9] += fqe * this_Filter->xState[count +  9];
        this_Filter->State[count + 10] += fqe * this_Filter->xState[count + 10];
        this_Filter->State[count + 11] += fqe * this_Filter->xState[count + 11];
        this_Filter->State[count + 12] += fqe * this_Filter->xState[count + 12];
        this_Filter->State[count + 13] += fqe * this_Filter->xState[count + 13];
        this_Filter->State[count + 14] += fqe * this_Filter->xState[count + 14];
        this_Filter->State[count + 15] += fqe * this_Filter->xState[count + 15];
        count += 16;
    }

    if ((SGNS.Filter_Order-count)>7)
    {
        this_Filter->State[count] += fqe * this_Filter->xState[count];
        this_Filter->State[count + 1] += fqe * this_Filter->xState[count + 1];
        this_Filter->State[count + 2] += fqe * this_Filter->xState[count + 2];
        this_Filter->State[count + 3] += fqe * this_Filter->xState[count + 3];
        this_Filter->State[count + 4] += fqe * this_Filter->xState[count + 4];
        this_Filter->State[count + 5] += fqe * this_Filter->xState[count + 5];
        this_Filter->State[count + 6] += fqe * this_Filter->xState[count + 6];
        this_Filter->State[count + 7] += fqe * this_Filter->xState[count + 7];
        count += 8;
    }

    if ((SGNS.Filter_Order-count)>3)
    {
        this_Filter->State[count] += fqe * this_Filter->xState[count];
        this_Filter->State[count + 1] += fqe * this_Filter->xState[count + 1];
        this_Filter->State[count + 2] += fqe * this_Filter->xState[count + 2];
        this_Filter->State[count + 3] += fqe * this_Filter->xState[count + 3];
        count += 4;
    }

    if ((SGNS.Filter_Order-count)>1)
    {
        this_Filter->State[count] += fqe * this_Filter->xState[count];
        this_Filter->State[count + 1] += fqe * this_Filter->xState[count + 1];
        count += 2;
    }

    if ((SGNS.Filter_Order-count)>0)
    {
        this_Filter->State[count] += fqe * this_Filter->xState[count];
        count++;
    }
}


void Process_Stored_Results_Cubic_AltFilter()
{
    double running_total = 0.0d;

    CubicInterp_Rec CI = {0.0d, 0.0d, 0.0d, 0.0d, 0.0d, 0.0d};

    CI.y2 = SGNSFFT.FFT_results_short_root[0];
    CI.y3 = SGNSFFT.FFT_results_short_root[1];
    CI.y1 = CI.y2 + CI.y2 - CI.y3;

    //==========================================================================
    // Fill desired shape curve.
    //==========================================================================
    for (int32_t ts_i = 0; ts_i <= SGNSFFT.FFT.length_half; ts_i++)
    {
        int32_t short_index = ts_i >> _factor_lookup_shift;

        int32_t factor_lookup = ts_i & _factor_lookup_mask;

        CI.x = SGNS.hi_fact[factor_lookup];

        if (factor_lookup == 0)
        {
            CI.y0 = CI.y1;
            CI.y1 = CI.y2;
            CI.y2 = CI.y3;

            if (short_index < history.FFT.length_half_m1)
                CI.y3 = SGNSFFT.FFT_results_short_root[short_index+2];
            else
                CI.y3 = CI.y2 + CI.y2 - CI.y1;
        }

        double result_short = CubicInterp(&CI);
        result_short *= ((result_short >= 0.0d) * (result_short * result_short));

        double result_combined = nroot((result_short * SGNS.Length_Factor[ts_i]) + (SGNSFFT.FFT_results_long[ts_i] * SGNSFFT.FFT_results_long_root[ts_i] * SGNS.Length_1_M_Factor[ts_i])) * SGNS.Gain[ts_i] + parameters.shaping.extra;

        running_total += result_combined;
        SGNS.Shape_Curve[ts_i] = result_combined;
        SGNS.Shape_Total[ts_i] = running_total;
    }
}


void Process_Stored_Results_AltFilter()
{
    double running_total = 0.0d;

    //==========================================================================
    // Fill desired shape curve.
    //==========================================================================
    for (int32_t ts_i = 0; ts_i <= SGNSFFT.FFT.length_half; ts_i++)
    {
        int32_t factor_lookup = ts_i & _factor_lookup_mask;
        int32_t short_index = ts_i >> _factor_lookup_shift;

        double result_short = (SGNSFFT.FFT_results_short_root[short_index] * SGNS.lo_fact[factor_lookup]) + (SGNSFFT.FFT_results_short_root[short_index + 1] * SGNS.hi_fact[factor_lookup]);
        result_short *= ((result_short >= 0.0d) * (result_short * result_short));

        double result_combined = nroot((result_short * SGNS.Length_Factor[ts_i]) + (SGNSFFT.FFT_results_long[ts_i] * SGNSFFT.FFT_results_long_root[ts_i] * SGNS.Length_1_M_Factor[ts_i])) * SGNS.Gain[ts_i] + parameters.shaping.extra;

        running_total += result_combined;
        SGNS.Shape_Curve[ts_i] = result_combined;
        SGNS.Shape_Total[ts_i] = running_total;
    }
}

void Process_Stored_Results_Cubic()
{
    double running_total = 0.0d;

    CubicInterp_Rec CI = {0.0d, 0.0d, 0.0d, 0.0d, 0.0d, 0.0d};

    CI.y3 = SGNSFFT.FFT_results_short_root[1];
    CI.y2 = SGNSFFT.FFT_results_short_root[0];
    CI.y1 = CI.y2 + CI.y2 - CI.y3;

    //==========================================================================
    // Fill desired shape curve.
    //==========================================================================
    for (int32_t ts_i = 0; ts_i <= SGNSFFT.FFT.length_half; ts_i++)
    {
        int32_t short_index = ts_i >> _factor_lookup_shift;

        int32_t factor_lookup = ts_i & _factor_lookup_mask;

        CI.x = SGNS.hi_fact[factor_lookup];

//                std::cerr << short_index << "; " << factor_lookup << "; " << CI.x << "; " << std::endl;

        if (factor_lookup == 0)
        {
            CI.y0 = CI.y1;
            CI.y1 = CI.y2;
            CI.y2 = CI.y3;

            if (short_index < history.FFT.length_half_m1)
                CI.y3 = SGNSFFT.FFT_results_short_root[short_index+2];
            else
                CI.y3 = CI.y2 + CI.y2 - CI.y1;
        }

        double result_short = CubicInterp(&CI);
        result_short *= ((result_short >= 0.0d) * (result_short));

        double result_combined = nroot((result_short * SGNS.Length_Factor[ts_i]) + (SGNSFFT.FFT_results_long[ts_i] * SGNS.Length_1_M_Factor[ts_i])) * SGNS.Gain[ts_i] + parameters.shaping.extra;

        running_total += result_combined;
        SGNS.Shape_Curve[ts_i] = result_combined;
        SGNS.Shape_Total[ts_i] = running_total;
    }
}

void Process_Stored_Results_Default()
{
    double running_total = 0.0d;

    //==========================================================================
    // Fill desired shape curve.
    //==========================================================================
    for (int32_t ts_i = 0; ts_i <= SGNSFFT.FFT.length_half; ts_i++)
    {
        int32_t factor_lookup = ts_i & _factor_lookup_mask;

        int32_t short_index = ts_i >> _factor_lookup_shift;

        double result_short = (SGNSFFT.FFT_results_short[short_index] * SGNS.lo_fact[factor_lookup]) + (SGNSFFT.FFT_results_short[short_index + 1] * SGNS.hi_fact[factor_lookup]);

        double result_combined = nroot((result_short * SGNS.Length_Factor[ts_i]) + (SGNSFFT.FFT_results_long[ts_i] * SGNS.Length_1_M_Factor[ts_i])) * SGNS.Gain[ts_i] + parameters.shaping.extra;

        running_total += result_combined;
        SGNS.Shape_Curve[ts_i] = result_combined;
        SGNS.Shape_Total[ts_i] = running_total;
    }


}


void Modify_Shape_Using_Average()
{
    if (SGNS.Use_Average)
    {
        double lower_average = SGNS.Shape_Total[SGNS.Upper_Freq_Bin] * OneOver[SGNS.Upper_Freq_Bin + 1];

        //==========================================================================
        // Calculate average up to Limit_Freq_Bin
        //==========================================================================

        for (int32_t ts_i = SGNS.Upper_Freq_Bin + 1; ts_i <= SGNSFFT.FFT.length_half; ts_i++)
        {
            double ts_x = SGNS.Shape_Curve[ts_i];

            SGNS.Shape_Curve[ts_i] = ts_x + ((ts_x < lower_average) * (lower_average - ts_x)) * SGNS.Average_Ratio[ts_i];

            SGNS.Shape_Total[ts_i] = SGNS.Shape_Total[ts_i - 1] + SGNS.Shape_Curve[ts_i];
        }
    }
}


double Fill_FFT_With_Warped_Spectrum_Cubic()
{
    double ts_x = 0;
    double ts_t = 0;
    CubicInterp_Rec CI;

    for (int32_t ts_i = 0; ts_i <= SGNSFFT.FFT.length_half; ts_i++)
    {
        CI.x = SGNS.Warp_Frac[ts_i];

        int32_t ts_j = SGNS.Warp_Int[ts_i];

        bool _is_not_zero = (ts_j > 0);
        bool _is_not_max = (ts_j < SGNSFFT.FFT.length_half_m1);

        if (_is_not_zero && _is_not_max)
        {
            CI.y0 = SGNS.Shape_Total[ts_j - 1];
            CI.y1 = SGNS.Shape_Total[ts_j];
            CI.y2 = SGNS.Shape_Total[ts_j + 1];
            CI.y3 = SGNS.Shape_Total[ts_j + 2];
        }
        else
            if (!_is_not_zero)
            {
                CI.y1 = SGNS.Shape_Total[0];
                CI.y2 = SGNS.Shape_Total[1];
                CI.y3 = SGNS.Shape_Total[2];
                CI.y0 = CI.y1 + CI.y1 - CI.y2;
                CI.y0 -= (CI.y0) * (CI.y0 < 0);
            }
            else
                if (ts_j == SGNSFFT.FFT.length_half_m1)
                {
                    CI.y0 = SGNS.Shape_Total[SGNSFFT.FFT.length_half_m1 - 1];
                    CI.y1 = SGNS.Shape_Total[SGNSFFT.FFT.length_half_m1];
                    CI.y2 = SGNS.Shape_Total[SGNSFFT.FFT.length_half_m1 + 1];
                    CI.y3 = CI.y2 + CI.y2 - CI.y1;
                    CI.y0 -= (CI.y3) * (CI.y3 < 0);
                }
                else
                {
                    CI.y0 = SGNS.Shape_Total[SGNSFFT.FFT.length_half_m1];
                    CI.y1 = SGNS.Shape_Total[SGNSFFT.FFT.length_half_m1 + 1];
                    CI.y2 = CI.y1 + CI.y1 - CI.y0;
                    CI.y2 -= (CI.y2) * (CI.y2 < 0);
                    CI.y3 = CI.y2 + CI.y2 - CI.y1;
                    CI.y3 -= (CI.y3) * (CI.y3 < 0);
                }

        double ts_z = CubicInterpMaxLimit(&CI);

        double ts_y = (ts_z - ts_x) * SGNS.Warp_Correction[ts_i];

        ts_y *= (ts_y > 0);

        ts_t += ts_y * ts_y;

        FFT_Array.DComplex[ts_i + ((SGNSFFT.FFT.length - 2 * ts_i) * ((ts_i & SGNSFFT.FFT.length_half_m1) > 0))] = real(ts_y);

        FFT_Array.DComplex[ts_i] = real(ts_y);

        ts_x = ts_z;
    }

    return ts_t;
}


double Fill_FFT_With_Warped_Spectrum_Default()
{
    double ts_x = 0;
    double ts_t = 0;

    for (int32_t ts_i = 0; ts_i <= SGNSFFT.FFT.length_half; ts_i++)
    {
        int32_t ts_j = SGNS.Warp_Int[ts_i];

        double ts_z = SGNS.Warp_1_M_Frac[ts_i] * SGNS.Shape_Total[ts_j]
                    + SGNS.Warp_Frac[ts_i] * SGNS.Shape_Total[ts_j + 1];

        double ts_y = (ts_z - ts_x) * SGNS.Warp_Correction[ts_i];

        ts_t += ts_y * ts_y;

        FFT_Array.DComplex[ts_i + ((SGNSFFT.FFT.length - ts_i - ts_i) * ((ts_i & SGNSFFT.FFT.length_half_m1) > 0))] = real(ts_y);

        FFT_Array.DComplex[ts_i] = real(ts_y);

        ts_x = ts_z;
    }

    return ts_t;
}


static bool Levinson(Filter_Rec* this_Filter)
{
    this_Filter->LD.A[0] = 1.0d;

    this_Filter->e = this_Filter->LD.r[0];

    int32_t m = 0;

    while ((this_Filter->e != 0) && (m < SGNS.Filter_Order))
    {
        m++;

        double err = this_Filter->LD.r[m];

        for (int32_t J = 1; J < m; J++)
        {
            int32_t K = m - J;
            this_Filter->LD.B[K] = this_Filter->LD.A[J];
            err += this_Filter->LD.A[J] * this_Filter->LD.r[K];
        }

        double km = -err / this_Filter->e;

        for (int32_t J = 1; J < m; J++)
        {
            this_Filter->LD.A[J] += km * this_Filter->LD.B[J];
        }

        this_Filter->LD.A[m] = km;

        this_Filter->k[m - 1] = km;

        this_Filter->e *= (1.0d - km * km);
    }

    this_Filter->Valid = (m == SGNS.Filter_Order);

    return this_Filter->Valid;
}


double Make_Filter(int32_t this_channel)
{
    SGNS.nFFT_plan.NumberOfBitsNeeded = SGNSFFT.FFT.bit_length;
    SGNS.nFFT_plan.FFT_Array = &FFT_Array;

    if ((parameters.shaping.active) && (!parameters.shaping.fixed))
    {
        SGNS.Control.Process_Stored_Results();

        if (SGNS.Use_Average)
          Modify_Shape_Using_Average();

        if (SGNS.Control.Fill_FFT_With_Warped_Spectrum())
        {
            if (FFTW_Initialised())
            {
                FFTW.Execute_C2C_New_Array(SGNS.FFTW_Plan_Inverse, &SGNS.nFFT_plan.DReal[0],  &SGNS.nFFT_plan.DReal[0]);
            }
            else
            {
                IFFT_DIT_Complex(&SGNS.nFFT_plan);
            }

            for (int32_t ts_i = 0; ts_i <= SGNS.Filter_Order; ts_i++)
            {
                SGNS.Filters[this_channel].LD.r[ts_i] = FFT_Array.DComplex[ts_i].Re;
            }

            if (Levinson(&SGNS.Filters[this_channel]))
            {
                return SGNS.Filters[this_channel].e;
            }

        }

        Stats.Skipped_Filters++;
    }

    SGNS.Filters[this_channel] = SGNS.Filters[MAX_CHANNELS];
    return SGNS.Filters[this_channel].e;
}


bool Filter_Valid(int32_t this_channel)
{
    return SGNS.Filters[this_channel].Valid;
}


double Filter_Error(int32_t this_channel)
{
    return SGNS.Filters[this_channel].e;
}


double warped_frequency(double Lambda, double Freq)
{
    tDComplex y = complex_exp(Freq);
    tDComplex z = (y + Lambda) / ((y * Lambda) + 1.0d);

    return std::atan2(z.Im, z.Re);
}


void nSGNS_Initialise()
{
    int32_t si_i, si_j, si_k;

    double average_factor;
    double average_power;
    double gain_amplitude;
    double gain_power;
    double Bins_to_Freq = Global.sample_rate * SGNSFFT.FFT.length_recip;


    //============================================================================
    // Set pointer to relevant FFTW plan (if available).
    //============================================================================
    if (FFTW_Initialised())
    {
        SGNS.FFTW_Plan_Inverse = FFTW.Plan_DFT_1d(SGNSFFT.FFT.length, FFT_Array.DReal, FFT_Array.DReal, FFTW_BACKWARD, 0x69);
    }


    //============================================================================
    // Set shaping factor for fixed noise shaping.
    //============================================================================
    double fixed_noise_shaping_factor = 1.0d;

    if (parameters.shaping.scale == -99)
    {
        parameters.shaping.scale = 1.0d;

        if (parameters.shaping.fixed)
            fixed_noise_shaping_factor = settings.fixed_noise_shaping_factor;
    }
    else
        if (parameters.shaping.fixed)
        {
            fixed_noise_shaping_factor = parameters.shaping.scale;
            parameters.shaping.scale = 1.0d;
        }
        else
            fixed_noise_shaping_factor = parameters.shaping.scale;


    //============================================================================
    // Use high frequency average tweak?
    // Set amplitude of default gain curve.
    //============================================================================
    SGNS.Use_Average  = true;

    average_power     = 2.0d;

    gain_amplitude    = 4.5d;
    gain_power        = 2.0d;

    if ((parameters.shaping.altfilter) && (parameters.shaping.interp_cubic))
    {
        average_factor    = 0.13d;

        SGNS.Control.Process_Stored_Results = Process_Stored_Results_Cubic_AltFilter;

        if (parameters.shaping.warp)
            SGNS.Control.Fill_FFT_With_Warped_Spectrum = Fill_FFT_With_Warped_Spectrum_Cubic;
        else
            SGNS.Control.Fill_FFT_With_Warped_Spectrum = Fill_FFT_With_Warped_Spectrum_Default;
    }
    else
        if (parameters.shaping.interp_cubic)
        {
            average_factor    = 1.087d;

            SGNS.Control.Process_Stored_Results = Process_Stored_Results_Cubic;

            if (parameters.shaping.warp)
                SGNS.Control.Fill_FFT_With_Warped_Spectrum = Fill_FFT_With_Warped_Spectrum_Cubic;
            else
                SGNS.Control.Fill_FFT_With_Warped_Spectrum = Fill_FFT_With_Warped_Spectrum_Default;
        }
        else
            if (parameters.shaping.altfilter)
            {
                average_factor    = 0.10515d;

                SGNS.Control.Process_Stored_Results = Process_Stored_Results_AltFilter;

                SGNS.Control.Fill_FFT_With_Warped_Spectrum = Fill_FFT_With_Warped_Spectrum_Default;
            }
            else
            {
                average_factor    = 1.0d;

                SGNS.Control.Process_Stored_Results = Process_Stored_Results_Default;

                SGNS.Control.Fill_FFT_With_Warped_Spectrum = Fill_FFT_With_Warped_Spectrum_Default;
            }

    //============================================================================
    // Zero number of skipped filters.
    //============================================================================
    Stats.Skipped_Filters = 0;


    //============================================================================
    // Initialise filters and create channel offset tables.
    //============================================================================
    std::memset((void*) &SGNS.Filters[0], 0, sizeof(SGNS.Filters));

    for (si_i = 0; si_i < MAX_CHANNELS; si_i++)
    {
        SGNS.Filters[si_i].One_Over_Minus_C = 1;
    }

    //============================================================================
    // Set filter order
    //============================================================================
    if (parameters.shaping.taps == -1)
    {
        SGNS.Filter_Order = 64;
    }
    else
    {
        SGNS.Filter_Order = parameters.shaping.taps;
    }

    SGNS.Filter_Order = std::min(MAX_FILTER_ORDER, std::min(SGNS.Filter_Order, SGNSFFT.FFT.length - 16));
    SGNS.Filter_Order_M1 = std::max(0, SGNS.Filter_Order - 1);

    //============================================================================
    // Set bin range points for "Average Ratio" curve.
    //============================================================================
    SGNS.Upper_Freq_Bin = nRoundEvenInt32(SGNSFFT.FFT.length * std::min(double(Global.upper_freq_limit) / Global.sample_rate, 46 * OneOver[96]));
    SGNS.Limit_Freq_Bin = nRoundEvenInt32(SGNSFFT.FFT.length * std::min((HUMAN_FREQ_LIMIT + 2000.0d * Global.Codec_Block.bit_shift) / Global.sample_rate, 46 * OneOver[96]));

    //============================================================================
    // Define "average ratio" curve.
    //============================================================================
    double this_average_ratio = (32 + (6 * Global.Codec_Block.bit_shift)) * OneOver[96] * average_factor;

    si_j = 0;
    si_k = SGNS.Upper_Freq_Bin;

    for (si_i = si_j; si_i <= si_k; si_i++)
    {
        SGNS.Average_Ratio[si_i] = 0.0d;
    }

    si_j = si_k;
    si_k = SGNS.Limit_Freq_Bin;

    for (si_i = si_j; si_i <= si_k; si_i++)
    {
        SGNS.Average_Ratio[si_i] = npowerx(double(si_i - si_j) * OneOver[std::max(si_k - si_j, 1)], average_power) * this_average_ratio;
    }

    si_j = si_k;
    si_k = (SGNSFFT.FFT.length_half - (SGNSFFT.FFT.length_half - SGNS.Limit_Freq_Bin) * OneOver[std::max(2, Global.Codec_Block.bit_shift + 2)]);

    for (si_i = si_j; si_i <= si_k; si_i++)
    {
        SGNS.Average_Ratio[si_i] = this_average_ratio;
    }

    si_j = si_k;
    si_k = SGNSFFT.FFT.length_half;

    for (si_i = si_j; si_i <= si_k; si_i++)
    {
        SGNS.Average_Ratio[si_i] = npowerx(double(si_k - si_i) * OneOver[std::max(si_k - si_j, 1)], average_power) * this_average_ratio;
    }

    //============================================================================
    // Calculate lambda value to use in frequency warping.
    //============================================================================
    double this_lambda = 0.0d;
    double this_freq = 0.0d;

    if (parameters.shaping.warp)
    {
        double this_delta = 0.5d;
        double this_freq_limit = std::min(double(HUMAN_FREQ_LIMIT), Global.sample_rate * 46.0d * OneOver[96]);
        this_lambda = 1.0d;
        double this_target = 46.0d * OneOver[96];

        do
        {
            this_freq = warped_frequency(this_lambda, this_target * TwoPi) * OneOverTwoPi * Global.sample_rate;

            if (this_freq < this_freq_limit)
            {
                this_lambda -= this_delta;
            }

            if (this_freq > this_freq_limit)
            {
                this_delta *= 0.50d;
                this_lambda += this_delta;
            }
        }
        while (!((this_delta < 1e-12d) || (abs(this_freq - Global.upper_freq_limit) < 1e-6d) || (this_lambda >= 1.0d) || (this_lambda < 0.0d) || (this_freq == this_freq_limit)));
    }

    SGNS.Lambda = (std::min(1.0d, std::max(0.0d, this_lambda)));
    SGNS.OneMinusLambda2 = 1.0d - nsqrd(SGNS.Lambda);
    SGNS.MinusLambda = -SGNS.Lambda;

    //============================================================================
    // Calculate lambda value to use in frequency warping.
    //============================================================================
    double last_warp = 0.0d;

    double this_total = 0.0d;

    for (si_i = 0; si_i <= SGNSFFT.FFT.length_half; si_i++)
    {
        //====================================================================
        // Zero stored spectral shapes.
        //====================================================================
        SGNSFFT.FFT_results_long[si_i] = 0.0d;
        SGNSFFT.FFT_results_short[si_i] = 0.0d;
        //====================================================================

        //====================================================================
        // Calculate unwarped and warped frequencies for each bin.
        //====================================================================
        double this_freq = si_i * SGNSFFT.FFT.length_recip;
        double this_warp = warped_frequency(SGNS.Lambda, this_freq * TwoPi) * OneOverTwoPi;
        //====================================================================

        double warped_frequency_bin = this_warp * SGNSFFT.FFT.length;

        SGNS.Warp_Int[si_i] = int(warped_frequency_bin);

        SGNS.Warp_Frac[si_i] = warped_frequency_bin - SGNS.Warp_Int[si_i];
        SGNS.Warp_1_M_Frac[si_i] = 1 - SGNS.Warp_Frac[si_i];

        SGNS.Length_Factor[si_i] = std::pow(SGNSFFT.FFT.length_half_recip * si_i, 1.5d);
        SGNS.Length_1_M_Factor[si_i] = 1.0d - SGNS.Length_Factor[si_i];
        SGNS.Length_Factor[si_i] = SGNS.Length_Factor[si_i] * PowersOf.Two[4];

        //====================================================================
        // Calculate warped-bin-width correction factors
        //====================================================================

        double warp_bin_delta = this_warp - last_warp;

        if (warp_bin_delta > 0.0d)
        {
            SGNS.Warp_Correction[si_i] = SGNSFFT.FFT.length_recip / warp_bin_delta;
        }
        else
        {
            SGNS.Warp_Correction[si_i] = 1.0d;
        }

        last_warp = this_warp;
        //====================================================================

        double this_frequency = si_i * Bins_to_Freq;
        double this_gain = 0;

        //====================================================================
        // Calculate unwarped gain curve to use with adaptive shaping method.
        //====================================================================
        if (this_frequency <= 16000.0d)
        {
            this_gain = 0.0d;
        }
        else
        {
            if (this_frequency >= 20000.0d)
            {
                this_gain = 1.0d;
            }
            else
            {
                this_gain = npowerx((this_frequency - 16000.0d) * OneOver[4000], gain_power);
            }

        }
        SGNS.Gain[si_i] = dB_Amplitude_Ratio(this_gain * gain_amplitude);
        //====================================================================


        //====================================================================
        // Calculate Fixed Noise Shaping curve
        //====================================================================
        si_j = 0;

        while (FIXED_GAIN_Freq[si_j + 1] < this_frequency)
        {
            si_j++;
        }

        double this_factor = ((this_frequency - FIXED_GAIN_Freq[si_j]) / (FIXED_GAIN_Freq[si_j + 1] - FIXED_GAIN_Freq[si_j]));
        double this_factor_1_M = 1.0d - this_factor;
        double this_factor_p3m1 = ((this_factor * this_factor * this_factor) - this_factor);
        double this_factor_1_M_p3m1 = ((this_factor_1_M * this_factor_1_M * this_factor_1_M) - this_factor_1_M);


        this_gain =  (this_factor_1_M * FNS_Gain[si_j]
                     + this_factor * FNS_Gain[si_j+1]
                     + this_factor_1_M_p3m1 * FNS_Spline[si_j]
                     + this_factor_p3m1 * FNS_Spline[si_j+1]) * 0.5d;

        this_total += dB_Amplitude_Ratio(this_gain * fixed_noise_shaping_factor);

        SGNS.Shape_Total[si_i] = this_total;

    }


//    if (parameters.shaping.interp_cubic)
//    {
//        SGNS.Warp_Correction[0] = (2.0d * (SGNS.Warp_Correction[1] - SGNS.Warp_Correction[2]));
//        SGNS.Warp_Correction[SGNSFFT.FFT.length_half] = 0.0d;
//    }


    for (int32_t ts_i = 0; ts_i < (1 << _factor_lookup_shift); ts_i++)
    {
        double hi_fact = ts_i * OneOver[(1 << _factor_lookup_shift)];
        double lo_fact = 1.0d - hi_fact;

        SGNS.hi_fact[ts_i] = hi_fact;
        SGNS.lo_fact[ts_i] = lo_fact;
    }


    //============================================================================
    // Use the fixed shaping filter when selected (i.e. skip "make_filter") and
    // for cases where the calculated adaptive shaping filter is invalid.
    // (skip Process_Stored_Results stage as desired curve in SGNS.Shape_Total)
    //============================================================================
    bool temp_bool = SGNS.Use_Average;
    SGNS.Use_Average = false;

    SGNS.Control.Fill_FFT_With_Warped_Spectrum();

    SGNS.nFFT_plan.NumberOfBitsNeeded = SGNSFFT.FFT.bit_length;
    SGNS.nFFT_plan.FFT_Array = &FFT_Array;

    if (FFTW_Initialised())
    {
        FFTW.Execute_C2C_New_Array(SGNS.FFTW_Plan_Inverse, &SGNS.nFFT_plan.DReal[0],  &SGNS.nFFT_plan.DReal[0]);
    }
    else
    {
        IFFT_DIT_Complex(&SGNS.nFFT_plan);
    }

    for (int32_t ts_i = 0; ts_i <= SGNS.Filter_Order; ts_i++)
    {
        SGNS.Filters[MAX_CHANNELS].LD.r[ts_i] = FFT_Array.DComplex[ts_i].Re;
    }

    Levinson(&SGNS.Filters[MAX_CHANNELS]);

    for (si_i = 0; si_i < MAX_CHANNELS; si_i++)
    {
        SGNS.Filters[si_i] = SGNS.Filters[MAX_CHANNELS];
    }

    SGNS.Use_Average = temp_bool;
}


void nSGNS_Cleanup()
{
    if ((FFTW_Initialised()) && (SGNS.FFTW_Plan_Inverse != nullptr))
    {
        FFTW.Destroy_Plan(SGNS.FFTW_Plan_Inverse);
    }
}
