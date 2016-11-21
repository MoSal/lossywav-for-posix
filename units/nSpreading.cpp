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

#include <cmath>

#include "nMaths.h"
#include "nSpreading.h"

// Globals
spreading_type spreading    __attribute__ ((aligned(16)));
double Skewing_Gain[MAX_FFT_LENGTH_HALF + 2];
double Frequency_Limits[SPREAD_ZONES + 2] = { 20, 1378.125, 3445.3125, 5512.5, 8268.75, 10335.9375, 12403.125, 14470.3125, 16000 };

const double _2x = -0.670485;
const double _3x =  0.780926; //-0.073165;
const double _4x =  1.543765;
double altspread_factor[PRECALC_ANALYSES+1] = {0,_4x,_3x,_2x,0,0,0,_2x};

void Spreading_Function(Results_Type* this_result)
{
    double alt_average = 0;

    spreading.new_minimum = Max_dB;
    spreading.old_minimum = Max_dB;

    for (int32_t sc_i = spreading.Bins.Lower[Current.FFT.bit_length]; sc_i <= spreading.Bins.Upper[Current.FFT.bit_length]; ++sc_i)
    {
        Spreading_rec sc_s = spreading.averages_ptr[Current.Analysis.number][sc_i];

        alt_average += this_result->Skewed[sc_i];

        if (!parameters.altspread)
        {
            double old_value = 0;

            for (int32_t sc_j = 0; sc_j < sc_s.width; ++sc_j)
            {
                old_value += this_result->Skewed[sc_i - sc_j];
            }
            old_value *= OneOver[sc_s.width];

            if (spreading.old_minimum > old_value)
            {
                spreading.old_minimum = old_value;
                process.old_min_bin = sc_i;
            }
        }

        double new_value = ((this_result->Skewed[sc_i - 1] + this_result->Skewed[sc_i + 1]) * spreading.Widths[sc_s.fractint] + this_result->Skewed[sc_i]) * spreading.divisors[sc_s.fractint];

        if (spreading.new_minimum > new_value)
        {
            spreading.new_minimum = new_value;
            process.new_min_bin = sc_i;
        }
    }

    spreading.alt_average = (nlog2(alt_average * spreading.Bins.Recip[Current.FFT.bit_length]) + Current.FFT.threshold_shift) * log10_2x20 + settings.noise_threshold_shift_average;
    spreading.old_minimum = (nlog2(spreading.old_minimum) + Current.FFT.threshold_shift) * log10_2x20 + settings.noise_threshold_shift_minimum;
    spreading.new_minimum = (nlog2(spreading.new_minimum) + Current.FFT.threshold_shift) * log10_2x20 + settings.noise_threshold_shift_minimum;

    if (parameters.altspread)
    {
        spreading.alt_average += altspread_factor[Current.Analysis.number] * log10_2x20;
    }
    else
    {
        if ((Current.FFT.length < PRECALC_ANALYSES_LENGTHS[SHORT_ANALYSIS]) && (parameters.fft.analyses > 3))
        {
            if (Current.FFT.length == PRECALC_ANALYSES_LENGTHS[IMPULSE_ANALYSIS])
                spreading.alt_average += 0.3072 * log10_2x20;
            else
                spreading.alt_average += 0.73065 * log10_2x20;
        }
    }
}


void Spreading_Memory_Allocate()
{
    for (int32_t sa_i = 1; sa_i <= PRECALC_ANALYSES; ++sa_i)
    {
        spreading.averages_ptr[sa_i] = new Spreading_rec[Current.Analysis.length[sa_i] + 1];
        spreading.Bark_Value[sa_i] = new double[(Current.Analysis.length[sa_i] >> 1) + 1];
    }
}


void Skewing_Curve_Init()
{
    double sa_lfb, sa_mfb, sa_tfb, sa_amp;

    sa_lfb = log10(LOWER_FREQ_LIMIT);
    sa_mfb = log10(MID_FREQ_LIMIT);

    if (parameters.fft.dccorrect)
        sa_amp = SKEWING_AMPLITUDE * 0.50f;
    else
        sa_amp = SKEWING_AMPLITUDE;

    for (int32_t sa_i = 0; sa_i <= PowersOf.TwoM1[MAX_FFT_BIT_LENGTH - 1] + 1; ++sa_i)
    {
        sa_tfb = log10((sa_i * PowersOf.TwoX[TWO_OFFSET + - MAX_FFT_BIT_LENGTH] * Global.sample_rate));

        if ((parameters.skewing) && (sa_tfb <= sa_mfb))
            Skewing_Gain[sa_i] = std::pow(2,(std::pow(std::sin(HalfPi * std::max(0.0, sa_tfb - sa_lfb) / (sa_mfb - sa_lfb)), 0.75) - 1) * sa_amp * OneOver[20] * log2_10);
        else
            Skewing_Gain[sa_i] = 1;
    }
}


void Spreading_Widths_Init()
{
    for (int32_t sa_i = 0; sa_i <= SPREADING_STEPS; ++sa_i)
    {
        spreading.Widths[sa_i] = float(sa_i * SPREADING_STEPS_RECIP);
        spreading.divisors[sa_i] = 1. / (1. + (2. * spreading.Widths[sa_i]));
    }
}


void Threshold_Index_Init()
{
    double this_reference_threshold;
    int32_t last_filled = 0;

    for (int32_t this_bit_to_remove = 1; this_bit_to_remove <= BITS_TO_CALCULATE; ++this_bit_to_remove)
    {
        // Constant below is used to achieve "closest" approximation to empirical noise calculations.

        this_reference_threshold = nlog2(3.0 * OneOver[64])
                                 + Global.Codec_Block.bits
                                 + this_bit_to_remove * 2.0
                                 + nlog2(1 + PowersOf.TwoX[TWO_OFFSET + 1 - this_bit_to_remove * 2])
                                 -0.00182241629983919;  // Empirical constant

        this_reference_threshold *= log10_2x20 * 0.5 * THRESHOLD_INDEX_SPREAD;

        while (last_filled < this_reference_threshold)
        {
            spreading.threshold_index[last_filled] = this_bit_to_remove - 1;
            ++ last_filled;
        }
    }

    while (last_filled < THRESHOLD_INDEX_SPREAD_MAX)
    {
        spreading.threshold_index[last_filled] = BITS_TO_CALCULATE - 1;
        ++ last_filled;
    }
}


void nSpreading_Init()
{
    Spreading_Memory_Allocate();
    Spreading_Widths_Init();
    Skewing_Curve_Init();
    Threshold_Index_Init();

    for (int32_t sa_i = 1; sa_i <= MAX_FFT_BIT_LENGTH; ++sa_i)
    {
        spreading.Bins.Lower[sa_i] = std::max(1, nRoundEvenInt32(double(Global.lower_freq_limit) / Global.sample_rate * PowersOf.TwoX[TWO_OFFSET + sa_i]) - 1);
        spreading.Bins.Upper[sa_i] = std::max(1, nRoundEvenInt32(double(Global.upper_freq_limit) / Global.sample_rate * PowersOf.TwoX[TWO_OFFSET + sa_i]) - 1);
        spreading.Bins.Recip[sa_i] = OneOver[spreading.Bins.Upper[sa_i] - spreading.Bins.Lower[sa_i] + 1];
        spreading.Bins.Mid[sa_i]   = std::max(1, nRoundEvenInt32(double(MID_FREQ_LIMIT) / Global.sample_rate * PowersOf.TwoX[TWO_OFFSET + sa_i]) - 1);
    }

    for (int32_t this_analysis_number = 1; this_analysis_number <= PRECALC_ANALYSES; ++this_analysis_number)
    {
        double this_spreading_function_width = SPREADING_FUNCTION_WIDTHS[this_analysis_number];

        if (parameters.altspread == true)
        {
            this_spreading_function_width = sqrt(1./PowersOf.TwoInt32[8-this_analysis_number]) * parameters.altspread_value;
        }

        Current.Analysis.number = this_analysis_number;
        Current.FFT = FFT_PreCalc_Data_Rec[Current.Analysis.bits[Current.Analysis.number]];

        for (int32_t sa_i = 0; sa_i <= SPREAD_ZONES + 1; ++sa_i)
        {
            spreading.Frequency_Bins[Current.Analysis.number][sa_i] = std::min(Current.FFT.length / 2 - 1, std::max(1, nRoundEvenInt32(Frequency_Limits[sa_i] / Global.sample_rate * Current.FFT.length)));
            if (spreading.Frequency_Bins[Current.Analysis.number][sa_i] > Current.FFT.length / 2 - 1)
                lossyWAVError(std::string("frequency too high : ") + NumToStr(Frequency_Limits[sa_i]* OneOver[1000], 3) + "kHz", 0x21);
        }

        for (int32_t sa_j = SPREAD_ZONES; sa_j >= 0; --sa_j)
            for (int32_t sa_i = spreading.Frequency_Bins[Current.Analysis.number][sa_j]; sa_i <= spreading.Frequency_Bins[Current.Analysis.number][sa_j + 1]; sa_i ++)
                spreading.averages_ptr[Current.Analysis.number][sa_i].width = std::min(uint16_t(sa_i), uint16_t(SPREADING_FUNCTION_ARRAY[Current.Analysis.number][sa_j]));

        for (int32_t sa_i = 0; sa_i <= spreading.Bins.Mid[Current.FFT.bit_length]; ++sa_i)
            spreading.averages_ptr[Current.Analysis.number][sa_i].fractint = 0;

        for (int32_t sa_i = spreading.Bins.Mid[Current.FFT.bit_length]; sa_i <= Current.FFT.length_half_m1; ++sa_i)
            spreading.averages_ptr[Current.Analysis.number][sa_i].fractint = uint16_t(SPREADING_STEPS * std::max(0.0, std::min(double(UPPER_FREQ_LIMIT - MID_FREQ_LIMIT), sa_i * double(Global.sample_rate) / Current.FFT.length - MID_FREQ_LIMIT)) / (UPPER_FREQ_LIMIT - MID_FREQ_LIMIT) * this_spreading_function_width);
    }
}

void nSpreading_Cleanup()
{
    for (int32_t sa_i = 1; sa_i <= PRECALC_ANALYSES; ++sa_i)
        if (spreading.averages_ptr[sa_i] != nullptr)
            delete[] spreading.averages_ptr[sa_i];
}
