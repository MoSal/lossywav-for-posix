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

#include "nCore.h"
#include "nMaths.h"
#include "nFFT.h"
#include "nFillFFT.h"

double* window_function[MAX_FFT_BIT_LENGTH + 1];

struct fill_fft_lookup_type
{
    int32_t block;
    int32_t offset;
}
fill_fft_lookup [MAX_BLOCK_SIZE * 4]; // 4 copies - Prev, Last, This, Next.


double Apply_Window_Function(FFT_Proc_Rec* this_FFT_plan, double Filled_Average)
{
    double ff_l = 0;
    double* this_window_function = window_function[this_FFT_plan->FFT->bit_length];
    Filled_Average*=this_FFT_plan->FFT->length_recip * settings.dccorrect_multiplier;

    for (int32_t ff_i = 0; ff_i < this_FFT_plan->FFT->length; ff_i++)
    {
        double ff_m = this_FFT_plan->DReal[ff_i] - Filled_Average;
        ff_l+= ff_m * ff_m;
        this_FFT_plan->DReal[ff_i] = ff_m * this_window_function[ff_i];
    }
    return nlog2(ff_l * this_FFT_plan->FFT->length_recip) * 0.50f;
}


double FillFFT_Input_From_WAVE(FFT_Proc_Rec* this_FFT_plan)
{
    double ff_k = 0;
    int32_t ff_n = (Global.Codec_Block.Size << 1) + this_FFT_plan->Task.block_start;

    for (int32_t ff_i = 0; ff_i<this_FFT_plan->FFT->length; ff_i++)
    {
        int32_t ff_j = ff_i + ff_n;
        double ff_m = AudioData.WAVEPTR[fill_fft_lookup[ff_j].block][Current.Channel][fill_fft_lookup[ff_j].offset].Integers[0] * settings.scaling_factor;
        ff_k+= ff_m;
        this_FFT_plan->DReal[ff_i] = ff_m;
    }

    return Apply_Window_Function(this_FFT_plan, ff_k);
}


double FillFFT_Input_From_BTRD(FFT_Proc_Rec* this_FFT_plan)
{
    double ff_m;
    double ff_k = 0;

    for (int32_t ff_i = 0; ff_i<this_FFT_plan->FFT->length; ff_i++)
    {
        int32_t ff_j = ff_i + this_FFT_plan->Task.block_start;

        if (ff_j < Global.Codec_Block.Size)
        {
            ff_j += Global.Codec_Block.Size << 1;
            ff_m = AudioData.BTRDPTR[fill_fft_lookup[ff_j].block][Current.Channel][fill_fft_lookup[ff_j].offset].Integers[0];
        }
        else
            ff_m = AudioData.WAVEPTR[NEXT_CODEC_BLOCK][Current.Channel][ff_j - Global.Codec_Block.Size].Integers[0] * settings.scaling_factor;

        ff_k += ff_m;
        this_FFT_plan->DReal[ff_i] = ff_m;
    }

    return Apply_Window_Function(this_FFT_plan, ff_k);
}


double FillFFT_Input_From_CORR(FFT_Proc_Rec* this_FFT_plan)
{
    double ff_m;
    double ff_k = 0;

    for (int32_t ff_i = 0; ff_i<this_FFT_plan->FFT->length; ff_i++)
    {
        int32_t ff_j = ff_i + this_FFT_plan->Task.block_start;

        if (ff_j < Global.Codec_Block.Size)
        {
            ff_j += Global.Codec_Block.Size << 1;
            ff_m = AudioData.CORRPTR[fill_fft_lookup[ff_j].block][Current.Channel][fill_fft_lookup[ff_j].offset].Integers[0] * settings.scaling_factor;
        }
        else
            ff_m = AudioData.WAVEPTR[NEXT_CODEC_BLOCK][Current.Channel][ff_j - Global.Codec_Block.Size].Integers[0] * settings.scaling_factor;

        ff_k += ff_m;
        this_FFT_plan->DReal[ff_i] = ff_m;
    }

    return Apply_Window_Function(this_FFT_plan, ff_k);
}


void nFillFFT_Init()
{
    //=========================================================================================================================
    // Create lookup table to allow FFT arrays to be quickly filled with samples.
    //=========================================================================================================================
    for (int32_t ff_i = 0; ff_i < Global.Codec_Block.Size * 4;  ++ff_i)
    {
        fill_fft_lookup[ff_i].block = (ff_i / Global.Codec_Block.Size);
        fill_fft_lookup[ff_i].offset = (ff_i % Global.Codec_Block.Size);
    }
    //=========================================================================================================================

    //=========================================================================================================================
    // Create window function lookup arrays.
    //=========================================================================================================================
    for (int32_t nf_i = 1; nf_i <= MAX_FFT_BIT_LENGTH; nf_i++)
    {
        window_function[nf_i] = new double[sizeof(double) << nf_i];
    }
    //=========================================================================================================================
    // Calculate window function values for MAX_FFT_LENGTH array.
    //=========================================================================================================================
    for (int32_t nf_j = 0; nf_j < MAX_FFT_LENGTH; nf_j++)
    {
        window_function[MAX_FFT_BIT_LENGTH][nf_j] = Hann(nf_j, MAX_FFT_BIT_LENGTH);
    }
    //=========================================================================================================================
    // Backfill shorter window function arrays quickly using shifted indices.
    //=========================================================================================================================
    for (int32_t nf_i=1; nf_i<MAX_FFT_BIT_LENGTH; nf_i++)
    {
        for (int32_t nf_j = 0; nf_j < 1 << nf_i; nf_j++)
        {
            window_function[nf_i][nf_j] = window_function[MAX_FFT_BIT_LENGTH][nf_j << (MAX_FFT_BIT_LENGTH - nf_i)];
        }
    }
    //=========================================================================================================================
}

void nFillFFT_Cleanup()
{
    //=========================================================================================================================
    // Delete window function arrays.
    //=========================================================================================================================
    for (int32_t ff_i = 0; ff_i <= MAX_FFT_BIT_LENGTH; ++ff_i)
    {
        if (window_function[ff_i] != nullptr)
        {
            delete[] window_function[ff_i];
        }
    }
    //=========================================================================================================================
}

