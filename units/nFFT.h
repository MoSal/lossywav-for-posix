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

#ifndef nFFT_h_
#define nFFT_h_

#include "nCore.h"

//==============================================================================
// declarations used in FFT testing in the nNoiseCalc unit.
//==============================================================================
extern void (*  FFT_DIT[32+1])(FFT_Proc_Rec*);
extern void (* IFFT_DIT[32+1])(FFT_Proc_Rec*);
//==============================================================================
void Shuffle_in_Place_DComplex(FFT_Proc_Rec*);

void Radix_16_DIT(FFT_Proc_Rec*);
void Radix_08_DIT(FFT_Proc_Rec*);
void Radix_04_DIT(FFT_Proc_Rec*);
void Radix_02_DIT(FFT_Proc_Rec*);
//==============================================================================

//====================================================================

inline long double Hann(int32_t Index, int32_t fft_bit_length)
{
    return 0.5 * (1.0 - std::cos(Index * TwoPi * PowersOf.TwoX[TWO_OFFSET - fft_bit_length]));
}

//==============================================================================
// Set the number of samples for FFT analysis;
//==============================================================================
int nFFT_Valid_Length(FFT_Proc_Rec*);

//==============================================================================
//  Calculates the (in place) n Fourier Transform of the array of complex
//  numbers represented by FFT_Input to produce the output complex
//  numbers in FFT_Array.
//==============================================================================
void FFT_DIT_Complex(FFT_Proc_Rec*);

//==============================================================================
//  Calculates the (in place) n Fourier Transform of the packed array of real
//  numbers represented by FFT_Input to produce the output complex
//  numbers in FFT_Array, then post-processes the results to unpack the
//  Real output in the full FFT_Array.
//==============================================================================
void FFT_DIT_Real(FFT_Proc_Rec*);

//==============================================================================
//  Calculates the (in place) Inverse n Fourier Transform of the array of
//  complex numbers represented by FFT_Input to produce the output
//  complex numbers in FFT_Array.
//==============================================================================
void IFFT_DIT_Complex(FFT_Proc_Rec*);

//==============================================================================
// Initialise lookup tables for internal FFT routines.
//==============================================================================
void nFFT_Init(int32_t desired_max_fft_bit_length);

//==============================================================================
// Free memory used by FFT related arrays
//==============================================================================
void nFFT_Cleanup();

#endif // nFFT_h
