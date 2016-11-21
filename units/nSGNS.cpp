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

#include "nSGNS.h"
#include "nParameter.h"

namespace { // anonymous

static const int32_t _factor_lookup_shift  __attribute__ ((aligned(16))) = (LONG_ANALYSIS - SHORT_ANALYSIS);
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
{ 15.6000, 15.5448, 15.5305, 15.5126, 15.4899, 15.4614, 15.4255, 15.3802, 15.3232, 15.2515, 15.1610, 15.0471, 14.9034, 14.7226, 14.4948, 14.2076, 13.8463, 13.3919, 12.8226, 12.1137, 11.2395, 10.1635,  8.8176,  7.1221,  4.9766, 2.1116, 0.2329, 1.2827, 5.4487, 11.1622, 16.6788, 22.2524, 30.2921, 60.0341, 70.0000, 70.0000, 70.0000};

static const float FNS_Spline[37]  __attribute__ ((aligned(16))) =
{  0.0000,  0.0204, -0.0018, -0.0024, -0.0029, -0.0037, -0.0047, -0.0059, -0.0073, -0.0094, -0.0117, -0.0149, -0.0185, -0.0235, -0.0297, -0.0371, -0.0465, -0.0574, -0.0698, -0.0826, -0.1009, -0.1350, -0.1748, -0.2250, -0.3598, 0.4932, 1.4643, 1.5581, 0.7738, -0.0985,  0.0285,  1.2331, 10.8512, -9.8881, -4.9830,  0.0000,  0.0000};

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
    FFT_Data_Rec       FFT                  __attribute__ ((aligned(16)));
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

void Warped_Lattice_Filter_Init(int32_t this_channel)
{
    Filter_Rec* this_Filter = &SGNS.Filters[this_channel];

    double o = 1;
    double u = 1;

    int32_t count=0;

    while (SGNS.Filter_Order > count)
    {
        this_Filter->State[count] = 0;
        this_Filter->xState[count] = u;
        u *= SGNS.MinusLambda;
        double A = u + this_Filter->k[count] * o;
        o += this_Filter->k[count] * u;
        u = A;
        count++;
    }

    this_Filter->One_Over_Minus_C = -1.0 / o;
}


double Warped_Lattice_Filter_Evaluate(int32_t this_channel)
{
    Filter_Rec* this_Filter = &SGNS.Filters[this_channel];

    double o = 0;
    double u = 0;

    int32_t count=0;

    while (SGNS.Filter_Order > count)
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
    Filter_Rec* this_Filter = &SGNS.Filters[this_channel];

    int32_t count=0;

    while (SGNS.Filter_Order > count)
    {
        this_Filter->State[count] += fqe * this_Filter->xState[count];
        count++;
    }
}


void Process_Stored_Results_Cubic_AltFilter()
{
    double running_total = 0.0;

    Results_Type* results_short = &results.WAVE[SHORT_ANALYSIS][Current.Channel];
    Results_Type* results_long  = &results.WAVE[LONG_ANALYSIS][Current.Channel];

    CubicInterp_Rec CI = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    CI.y2 = results_short->SGNSRoot[0];
    CI.y3 = results_short->SGNSRoot[1];
    CI.y1 = CI.y2 + CI.y2 - CI.y3;

    //==========================================================================
    // Fill desired shape curve.
    //==========================================================================
    for (int32_t ts_i = 0; ts_i <= SGNS.FFT.length_half; ts_i++)
    {
        int32_t short_index = ts_i >> _factor_lookup_shift;

        int32_t factor_lookup = ts_i & _factor_lookup_mask;

        CI.x = SGNS.hi_fact[factor_lookup];

        if (factor_lookup == 0)
        {
            CI.y0 = CI.y1;
            CI.y1 = CI.y2;
            CI.y2 = CI.y3;

            if (short_index < settings.analysis[SHORT_ANALYSIS].FFT.length_half_m1)
                CI.y3 = results_short->SGNSRoot[short_index+2];
            else
                CI.y3 = CI.y2 + CI.y2 - CI.y1;
        }

        double result_short = CubicInterp(&CI);
        result_short *= ((result_short >= 0.0) * (result_short * result_short));

        double result_combined = nroot((result_short * SGNS.Length_Factor[ts_i]) + (results_long->SGNSUnity[ts_i] * results_long->SGNSRoot[ts_i] * SGNS.Length_1_M_Factor[ts_i])) * SGNS.Gain[ts_i] + parameters.shaping.extra;

        running_total += result_combined;
        SGNS.Shape_Curve[ts_i] = result_combined;
        SGNS.Shape_Total[ts_i] = running_total;
    }
}


void Process_Stored_Results_AltFilter()
{
    double running_total = 0.0;

    Results_Type* results_short = &results.WAVE[SHORT_ANALYSIS][Current.Channel];
    Results_Type* results_long  = &results.WAVE[LONG_ANALYSIS][Current.Channel];

    //==========================================================================
    // Fill desired shape curve.
    //==========================================================================
    for (int32_t ts_i = 0; ts_i <= SGNS.FFT.length_half; ts_i++)
    {
        int32_t factor_lookup = ts_i & _factor_lookup_mask;
        int32_t short_index = ts_i >> _factor_lookup_shift;

        double result_short = (results_short->SGNSRoot[short_index] * SGNS.lo_fact[factor_lookup]) + (results_short->SGNSRoot[short_index + 1] * SGNS.hi_fact[factor_lookup]);
        result_short *= ((result_short >= 0.0) * (result_short * result_short));

        double result_combined = nroot((result_short * SGNS.Length_Factor[ts_i]) + (results_long->SGNSUnity[ts_i] * results_long->SGNSRoot[ts_i] * SGNS.Length_1_M_Factor[ts_i])) * SGNS.Gain[ts_i] + parameters.shaping.extra;

        running_total += result_combined;
        SGNS.Shape_Curve[ts_i] = result_combined;
        SGNS.Shape_Total[ts_i] = running_total;
    }
}

void Process_Stored_Results_Cubic()
{
    double running_total = 0.0;

    Results_Type* results_short = &results.WAVE[SHORT_ANALYSIS][Current.Channel];
    Results_Type* results_long  = &results.WAVE[LONG_ANALYSIS][Current.Channel];

    CubicInterp_Rec CI = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    CI.y3 = results_short->SGNSRoot[1];
    CI.y2 = results_short->SGNSRoot[0];
    CI.y1 = CI.y2 + CI.y2 - CI.y3;

    //==========================================================================
    // Fill desired shape curve.
    //==========================================================================
    for (int32_t ts_i = 0; ts_i <= SGNS.FFT.length_half; ts_i++)
    {
        int32_t short_index = ts_i >> _factor_lookup_shift;

        int32_t factor_lookup = ts_i & _factor_lookup_mask;

        CI.x = SGNS.hi_fact[factor_lookup];

        if (factor_lookup == 0)
        {
            CI.y0 = CI.y1;
            CI.y1 = CI.y2;
            CI.y2 = CI.y3;

            if (short_index < settings.analysis[SHORT_ANALYSIS].FFT.length_half_m1)
                CI.y3 = results_short->SGNSRoot[short_index+2];
            else
                CI.y3 = CI.y2 + CI.y2 - CI.y1;
        }

        double result_short = CubicInterp(&CI);
        result_short *= ((result_short >= 0.0) * (result_short));

        double result_combined = nroot((result_short * SGNS.Length_Factor[ts_i]) + (results_long->SGNSUnity[ts_i] * SGNS.Length_1_M_Factor[ts_i])) * SGNS.Gain[ts_i] + parameters.shaping.extra;

        running_total += result_combined;
        SGNS.Shape_Curve[ts_i] = result_combined;
        SGNS.Shape_Total[ts_i] = running_total;
    }
}

void Process_Stored_Results_Default()
{
    double running_total = 0.0;

    Results_Type* results_short = &results.WAVE[SHORT_ANALYSIS][Current.Channel];
    Results_Type* results_long  = &results.WAVE[LONG_ANALYSIS][Current.Channel];

    //==========================================================================
    // Fill desired shape curve.
    //==========================================================================
    for (int32_t ts_i = 0; ts_i <= SGNS.FFT.length_half; ts_i++)
    {
        int32_t factor_lookup = ts_i & _factor_lookup_mask;

        int32_t short_index = ts_i >> _factor_lookup_shift;

        double result_short = (results_short->SGNSUnity[short_index] * SGNS.lo_fact[factor_lookup]) + (results_short->SGNSUnity[short_index + 1] * SGNS.hi_fact[factor_lookup]);

        double result_combined = nroot((result_short * SGNS.Length_Factor[ts_i]) + (results_long->SGNSUnity[ts_i] * SGNS.Length_1_M_Factor[ts_i])) * SGNS.Gain[ts_i] + parameters.shaping.extra;

        running_total += result_combined;
        SGNS.Shape_Curve[ts_i] = result_combined;
        SGNS.Shape_Total[ts_i] = running_total;
    }
}


void Process_Stored_Results_Hybrid()
{
    {
        Results_Type* this_result = &results.WAVE[LONG_ANALYSIS][Current.Channel];

        for (int32_t ts_i = 0; ts_i <= SGNS.FFT.length_half; ts_i++)
        {
            SGNS.Shape_Curve[ts_i] = this_result->SGNSHybrid[ts_i] * LONG_ANALYSIS;
        }
    }

    int32_t analyses_sum = LONG_ANALYSIS;

    for (int32_t this_analysis = 1; this_analysis < PRECALC_ANALYSES; this_analysis++)
    {
        if (settings.analysis[this_analysis].active)
        {
            Results_Type* this_result = &results.WAVE[this_analysis][Current.Channel];

            analyses_sum += this_analysis;

            int32_t this_factor_shift = LONG_ANALYSIS - this_analysis;

            int32_t this_factor_mask = PowersOf.TwoM1[this_factor_shift];

            double this_factor_recip = PowersOf.TwoX[TWO_OFFSET + -this_factor_shift];

            for (int32_t ts_i = 0; ts_i <= SGNS.FFT.length_half; ts_i++)
            {
                double factor = (ts_i & this_factor_mask) * this_factor_recip;

                int32_t index = ts_i >> this_factor_shift;

                SGNS.Shape_Curve[ts_i] += ((this_result->SGNSHybrid[index] * (1.0 - factor)) + (this_result->SGNSHybrid[index + 1] * factor)) * this_analysis;
            }
        }
    }

    double this_divisor = OneOver[analyses_sum];

    double running_total = 0.0;

    for (int32_t ts_i = 0; ts_i <= SGNS.FFT.length_half; ts_i++)
    {
        double result =  nroot(SGNS.Shape_Curve[ts_i] * this_divisor) * SGNS.Gain[ts_i] + parameters.shaping.extra;

        running_total += result;

        SGNS.Shape_Curve[ts_i] = result;

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

        for (int32_t ts_i = SGNS.Upper_Freq_Bin + 1; ts_i <= SGNS.FFT.length_half; ts_i++)
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

    for (int32_t ts_i = 0; ts_i <= SGNS.FFT.length_half; ts_i++)
    {
        CI.x = SGNS.Warp_Frac[ts_i];

        int32_t ts_j = SGNS.Warp_Int[ts_i];

        bool _is_not_zero = (ts_j > 0);
        bool _is_not_max = (ts_j < SGNS.FFT.length_half_m1);

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
                if (ts_j == SGNS.FFT.length_half_m1)
                {
                    CI.y0 = SGNS.Shape_Total[SGNS.FFT.length_half_m1 - 1];
                    CI.y1 = SGNS.Shape_Total[SGNS.FFT.length_half_m1];
                    CI.y2 = SGNS.Shape_Total[SGNS.FFT.length_half_m1 + 1];
                    CI.y3 = CI.y2 + CI.y2 - CI.y1;
                    CI.y3 -= (CI.y3) * (CI.y3 < 0);
                }
                else
                {
                    CI.y0 = SGNS.Shape_Total[SGNS.FFT.length_half_m1];
                    CI.y1 = SGNS.Shape_Total[SGNS.FFT.length_half_m1 + 1];
                    CI.y2 = CI.y1 + CI.y1 - CI.y0;
                    CI.y2 -= (CI.y2) * (CI.y2 < 0);
                    CI.y3 = CI.y2 + CI.y2 - CI.y1;
                    CI.y3 -= (CI.y3) * (CI.y3 < 0);
                }

        double ts_z = CubicInterpMaxLimit(&CI);

        double ts_y = (ts_z - ts_x) * SGNS.Warp_Correction[ts_i];

        ts_y *= (ts_y > 0);

        ts_t += ts_y * ts_y;

        FFT_Array.DComplex[ts_i + ((SGNS.FFT.length - 2 * ts_i) * ((ts_i & SGNS.FFT.length_half_m1) > 0))] = DComplex(ts_y, 0.);

        FFT_Array.DComplex[ts_i] = DComplex(ts_y, 0.);

        ts_x = ts_z;
    }

    return ts_t;
}


double Fill_FFT_With_Warped_Spectrum_Default()
{
    double ts_x = 0;
    double ts_t = 0;

    for (int32_t ts_i = 0; ts_i <= SGNS.FFT.length_half; ts_i++)
    {
        int32_t ts_j = SGNS.Warp_Int[ts_i];

        double ts_z = SGNS.Warp_1_M_Frac[ts_i] * SGNS.Shape_Total[ts_j]
                    + SGNS.Warp_Frac[ts_i] * SGNS.Shape_Total[ts_j + 1];

        double ts_y = (ts_z - ts_x) * SGNS.Warp_Correction[ts_i];

        ts_t += ts_y * ts_y;

        FFT_Array.DComplex[ts_i + ((SGNS.FFT.length - ts_i - ts_i) * ((ts_i & SGNS.FFT.length_half_m1) > 0))] = DComplex(ts_y, 0.);

        FFT_Array.DComplex[ts_i] = DComplex(ts_y, 0.);

        ts_x = ts_z;
    }

    return ts_t;
}


static bool Levinson(Filter_Rec* this_Filter)
{
    this_Filter->LD.A[0] = 1.0;

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

        this_Filter->e *= (1.0 - km * km);
    }

    this_Filter->Valid = (m == SGNS.Filter_Order);

    return this_Filter->Valid;
}


double Make_Filter(int32_t this_channel)
{
    SGNS.nFFT_plan.NumberOfBitsNeeded = SGNS.FFT.bit_length;
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
    tDComplex z = (y + Lambda) / ((y * Lambda) + 1.0);

    return std::atan2(z.Im, z.Re);
}


void nSGNS_Initialise()
{
    SGNS.FFT = FFT_PreCalc_Data_Rec[settings.analysis[LONG_ANALYSIS].FFT.bit_length];

    int32_t si_i, si_j, si_k;

    double average_factor;
    double average_power;
    double gain_amplitude;
    double gain_power;
    double Bins_to_Freq = Global.sample_rate * SGNS.FFT.length_recip;


    //============================================================================
    // Set pointer to relevant FFTW plan (if available).
    //============================================================================
    if (FFTW_Initialised())
    {
        SGNS.FFTW_Plan_Inverse = FFTW.Plan_DFT_1d(SGNS.FFT.length, FFT_Array.DReal, FFT_Array.DReal, FFTW_BACKWARD, 0x69);
    }


    //============================================================================
    // Set shaping factor for fixed noise shaping.
    //============================================================================
    double fixed_noise_shaping_factor = 1.0;

    if (parameters.shaping.scale == -99)
    {
        parameters.shaping.scale = 1.0;

        if (parameters.shaping.fixed)
            fixed_noise_shaping_factor = settings.fixed_noise_shaping_factor;
    }
    else
        if (parameters.shaping.fixed)
        {
            fixed_noise_shaping_factor = parameters.shaping.scale;
            parameters.shaping.scale = 1.0;
        }
        else
            fixed_noise_shaping_factor = parameters.shaping.scale;


    //============================================================================
    // Use high frequency average tweak?
    // Set amplitude of default gain curve.
    //============================================================================
    SGNS.Use_Average  = true;

    average_power     = 2.0;

    gain_amplitude    = 4.5;
    gain_power        = 2.0;

    if (parameters.shaping.hybrid)
    {
        SGNS.Control.Process_Stored_Results = Process_Stored_Results_Hybrid;

        average_factor    = 1.0;

        if (parameters.shaping.interp_cubic)
            SGNS.Control.Fill_FFT_With_Warped_Spectrum = Fill_FFT_With_Warped_Spectrum_Cubic;
        else
            SGNS.Control.Fill_FFT_With_Warped_Spectrum = Fill_FFT_With_Warped_Spectrum_Default;
    }
    else
        if (parameters.shaping.altfilter)
        {
            SGNS.Control.Process_Stored_Results = Process_Stored_Results_Cubic_AltFilter;

            if (parameters.shaping.interp_cubic)
            {
                SGNS.Control.Process_Stored_Results = Process_Stored_Results_Cubic_AltFilter;
                SGNS.Control.Fill_FFT_With_Warped_Spectrum = Fill_FFT_With_Warped_Spectrum_Cubic;
                average_factor    = 0.13;
            }
            else
            {
                SGNS.Control.Process_Stored_Results = Process_Stored_Results_AltFilter;
                SGNS.Control.Fill_FFT_With_Warped_Spectrum = Fill_FFT_With_Warped_Spectrum_Default;
                average_factor    = 0.10515;
            }
        }
        else
            if (parameters.shaping.interp_cubic)
            {
                SGNS.Control.Process_Stored_Results = Process_Stored_Results_Cubic;

                SGNS.Control.Fill_FFT_With_Warped_Spectrum = Fill_FFT_With_Warped_Spectrum_Cubic;

                average_factor    = 1.087;
            }
            else
            {
                SGNS.Control.Process_Stored_Results = Process_Stored_Results_Default;

                SGNS.Control.Fill_FFT_With_Warped_Spectrum = Fill_FFT_With_Warped_Spectrum_Default;

                average_factor    = 1.0;
            }

    if (parameters.shaping.average != -99)
    {
        average_factor = parameters.shaping.average;
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

    SGNS.Filter_Order = std::min(MAX_FILTER_ORDER, std::min(SGNS.Filter_Order, SGNS.FFT.length - 16));
    SGNS.Filter_Order_M1 = std::max(0, SGNS.Filter_Order - 1);

    //============================================================================
    // Set bin range points for "Average Ratio" curve.
    //============================================================================
    SGNS.Upper_Freq_Bin = nRoundEvenInt32(SGNS.FFT.length * std::min(double(Global.upper_freq_limit) / Global.sample_rate, 46 * OneOver[96]));
    SGNS.Limit_Freq_Bin = nRoundEvenInt32(SGNS.FFT.length * std::min((HUMAN_FREQ_LIMIT + 2000.0 * Global.Codec_Block.bit_shift) / Global.sample_rate, 46 * OneOver[96]));

    //============================================================================
    // Define "average ratio" curve.
    //============================================================================
    double this_average_ratio = (32 + (6 * Global.Codec_Block.bit_shift)) * OneOver[96] * average_factor;

    si_j = 0;
    si_k = SGNS.Upper_Freq_Bin;

    for (si_i = si_j; si_i <= si_k; si_i++)
    {
        SGNS.Average_Ratio[si_i] = 0.0;
    }

    si_j = si_k;
    si_k = SGNS.Limit_Freq_Bin;

    for (si_i = si_j; si_i <= si_k; si_i++)
    {
        SGNS.Average_Ratio[si_i] = npowerx(double(si_i - si_j) * OneOver[std::max(si_k - si_j, 1)], average_power) * this_average_ratio;
    }

    si_j = si_k;
    si_k = (SGNS.FFT.length_half - (SGNS.FFT.length_half - SGNS.Limit_Freq_Bin) * OneOver[std::max(2, Global.Codec_Block.bit_shift + 2)]);

    for (si_i = si_j; si_i <= si_k; si_i++)
    {
        SGNS.Average_Ratio[si_i] = this_average_ratio;
    }

    si_j = si_k;
    si_k = SGNS.FFT.length_half;

    for (si_i = si_j; si_i <= si_k; si_i++)
    {
        SGNS.Average_Ratio[si_i] = npowerx(double(si_k - si_i) * OneOver[std::max(si_k - si_j, 1)], average_power) * this_average_ratio;
    }

    //============================================================================
    // Calculate lambda value to use in frequency warping.
    //============================================================================
    double this_lambda = 0.0;
    double this_freq = 0.0;

    if (parameters.shaping.warp)
    {
        double this_delta = 0.5;
        double this_freq_limit = std::min(double(HUMAN_FREQ_LIMIT), Global.sample_rate * 46.0 * OneOver[96]);
        this_lambda = 1.0;
        double this_target = 46.0 * OneOver[96];

        do
        {
            this_freq = warped_frequency(this_lambda, this_target * TwoPi) * OneOverTwoPi * Global.sample_rate;

            if (this_freq < this_freq_limit)
            {
                this_lambda -= this_delta;
            }

            if (this_freq > this_freq_limit)
            {
                this_delta *= 0.50;
                this_lambda += this_delta;
            }
        }
        while (!((this_delta < 1e-12) || (abs(this_freq - Global.upper_freq_limit) < 1e-6) || (this_lambda >= 1.0) || (this_lambda < 0.0) || (this_freq == this_freq_limit)));
    }

    SGNS.Lambda = (std::min(1.0, std::max(0.0, this_lambda)));
    SGNS.OneMinusLambda2 = 1.0 - nsqrd(SGNS.Lambda);
    SGNS.MinusLambda = -SGNS.Lambda;

    //============================================================================
    // Calculate lambda value to use in frequency warping.
    //============================================================================
    double last_warp = 0.0;

    double this_total = 0.0;

    for (si_i = 0; si_i <= SGNS.FFT.length_half; si_i++)
    {
        //====================================================================
        // Calculate unwarped and warped frequencies for each bin.
        //====================================================================
        double this_freq = si_i * SGNS.FFT.length_recip;
        double this_warp = warped_frequency(SGNS.Lambda, this_freq * TwoPi) * OneOverTwoPi;
        //====================================================================

        double warped_frequency_bin = this_warp * SGNS.FFT.length;

        SGNS.Warp_Int[si_i] = int(warped_frequency_bin);

        SGNS.Warp_Frac[si_i] = warped_frequency_bin - SGNS.Warp_Int[si_i];
        SGNS.Warp_1_M_Frac[si_i] = 1 - SGNS.Warp_Frac[si_i];

        SGNS.Length_Factor[si_i] = std::pow(SGNS.FFT.length_half_recip * si_i, 1.5);
        SGNS.Length_1_M_Factor[si_i] = 1.0 - SGNS.Length_Factor[si_i];
        SGNS.Length_Factor[si_i] = SGNS.Length_Factor[si_i] * PowersOf.TwoX[TWO_OFFSET + 4];

        //====================================================================
        // Calculate warped-bin-width correction factors
        //====================================================================

        double warp_bin_delta = this_warp - last_warp;

        if (warp_bin_delta > 0.0)
        {
            SGNS.Warp_Correction[si_i] = SGNS.FFT.length_recip / warp_bin_delta;
        }
        else
        {
            SGNS.Warp_Correction[si_i] = 1.0;
        }

        last_warp = this_warp;
        //====================================================================

        double this_frequency = si_i * Bins_to_Freq;
        double this_gain = 0;

        //====================================================================
        // Calculate unwarped gain curve to use with adaptive shaping method.
        //====================================================================
        if (this_frequency <= 16000.0)
        {
            this_gain = 0.0;
        }
        else
        {
            if (this_frequency >= 20000.0)
            {
                this_gain = 1.0;
            }
            else
            {
                this_gain = npowerx((this_frequency - 16000.0) * OneOver[4000], gain_power);
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
        double this_factor_1_M = 1.0 - this_factor;
        double this_factor_p3m1 = ((this_factor * this_factor * this_factor) - this_factor);
        double this_factor_1_M_p3m1 = ((this_factor_1_M * this_factor_1_M * this_factor_1_M) - this_factor_1_M);


        this_gain =  (this_factor_1_M * FNS_Gain[si_j]
                     + this_factor * FNS_Gain[si_j+1]
                     + this_factor_1_M_p3m1 * FNS_Spline[si_j]
                     + this_factor_p3m1 * FNS_Spline[si_j+1]) * 0.5;

        this_total += dB_Amplitude_Ratio(this_gain * fixed_noise_shaping_factor);

        SGNS.Shape_Total[si_i] = this_total;
    }


    for (int32_t ts_i = 0; ts_i < (1 << _factor_lookup_shift); ts_i++)
    {
        double hi_fact = ts_i * OneOver[(1 << _factor_lookup_shift)];
        double lo_fact = 1.0 - hi_fact;

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

    SGNS.nFFT_plan.NumberOfBitsNeeded = SGNS.FFT.bit_length;
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
