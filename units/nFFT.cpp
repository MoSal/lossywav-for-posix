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
    along with this program.  If not, see <http:www.gnu.org/licenses/>.

    Contact: lossywav <at> hotmail <dot> co <dot> uk

==============================================================================

     fourier.pas  -  Don Cross <dcross <at> intersrv <dot> com>

     This is a Turbo Pascal Unit for calculating the n Fourier Transform
     (FFT) and the Inverse n Fourier Transform (IFFT).
     Visit the following URL for the latest version of this code.
     This page also has a C/C++ version, and a brief discussion of the
     theory behind the FFT algorithm.

        http:www.intersrv.com/~dcross/fft.html#pascal [!!Dead link!!]

     Revision history [most recent first]:

 2012 November [Nick Currie]
      Further optimisation of the FFT code using complex operators;
      Radix Sixteen FFT forward and inverse added;
      One step Radix 2, 4, 8 and 16 FFT forward added;
      Bit reversal arrays shortened to bit-reversal lists containing only
      those elements that require to be exchanged - saves jumps at the
      expense of extra memory footprint.

 2012 October [Nick Currie]
      Optimisation of FFT code;
      Radix Four and Eight added to inverse FFT;
      Bit-reversal arrays added for each FFT length.

 2012 Aug [Tyge Løvset (tycho), Aug. 2012]
     Translated to C++ from Delphi

 2009 May [Nick Currie]
     Radix Eight added.

 2009 April 24 [Nick Currie]
     Radix Four added.

 2008 September 29 [Nick Currie]
     All code transcoded back to Delphi from IA-32/x87.

 2008 May 12 [Nick Currie]
     All code transcoded to IA-32/x87. FFT_Real implemented to allow calculation
     of a packed Real array of length 2N to be calculated in a Complex array of
     length N, with subsequent untangling. ReversedBits, Ai & Ar lookup tables
     implemented. Calc_Frequency, FFT_Integer and FFT_Integer_Cleanup removed
     as not used in lossyWAV.

 1996 December 11 [Don Cross]
     Improved documentation of the procedure CalcFrequency.
     Fixed some messed up comments in procedure IFFT.

 1996 December 6 [Don Cross]
     Made procedure 'fft_integer' more efficient when buffer size changes
     in successive calls:  the buffer is now only resized when the input
     has more samples, not a differing number of samples.
     Also changed the way 'fft_integer_cleanup' works so that it is
     more "bullet-proof".

 1996 December 4 [Don Cross]
     Adding the procedure 'CalcFrequency', which calculates the FFT
     at a specific frequency index p=0..n-1, instead of the whole
     FFT.  This is O(n^2) instead of O(n*log(n)).

 1996 November 30 [Don Cross]
     Adding a routine to allow FFT of an input array of integers.
     It is called 'fft_integer'.

 1996 November 18 [Don Cross]
     Added some comments.

 1996 November 17 [Don Cross]
     Wrote and debugged first version.

==============================================================================**/
#include "math.h"
#include "nMaths.h"
#include "nFFT.h"
#include "nCore.h"

//#define FACTOR_16
//#define FACTOR_8_4
#define FACTOR_4_8
//#define FACTOR_4_2
//#define FACTOR_2

static tDComplex* A_Arr [32+1]   __attribute__ ((aligned(16)));
static tDComplex* A_Arr_Conj [32+1]   __attribute__ ((aligned(16)));

static int32_t* RevBits [32+1]   __attribute__ ((aligned(16)));

static double cos_pi_over_32[32] __attribute__ ((aligned(16)));

static int32_t nFFT_MAX_FFT_BIT_LENGTH  __attribute__ ((aligned(16))) = 0;

static const unsigned char BitReversedLookupTable[] __attribute__ ((aligned(16))) =
{ 0x00, 0x80, 0x40, 0xC0, 0x20, 0xA0, 0x60, 0xE0, 0x10, 0x90, 0x50, 0xD0, 0x30, 0xB0, 0x70, 0xF0,
  0x08, 0x88, 0x48, 0xC8, 0x28, 0xA8, 0x68, 0xE8, 0x18, 0x98, 0x58, 0xD8, 0x38, 0xB8, 0x78, 0xF8,
  0x04, 0x84, 0x44, 0xC4, 0x24, 0xA4, 0x64, 0xE4, 0x14, 0x94, 0x54, 0xD4, 0x34, 0xB4, 0x74, 0xF4,
  0x0C, 0x8C, 0x4C, 0xCC, 0x2C, 0xAC, 0x6C, 0xEC, 0x1C, 0x9C, 0x5C, 0xDC, 0x3C, 0xBC, 0x7C, 0xFC,
  0x02, 0x82, 0x42, 0xC2, 0x22, 0xA2, 0x62, 0xE2, 0x12, 0x92, 0x52, 0xD2, 0x32, 0xB2, 0x72, 0xF2,
  0x0A, 0x8A, 0x4A, 0xCA, 0x2A, 0xAA, 0x6A, 0xEA, 0x1A, 0x9A, 0x5A, 0xDA, 0x3A, 0xBA, 0x7A, 0xFA,
  0x06, 0x86, 0x46, 0xC6, 0x26, 0xA6, 0x66, 0xE6, 0x16, 0x96, 0x56, 0xD6, 0x36, 0xB6, 0x76, 0xF6,
  0x0E, 0x8E, 0x4E, 0xCE, 0x2E, 0xAE, 0x6E, 0xEE, 0x1E, 0x9E, 0x5E, 0xDE, 0x3E, 0xBE, 0x7E, 0xFE,
  0x01, 0x81, 0x41, 0xC1, 0x21, 0xA1, 0x61, 0xE1, 0x11, 0x91, 0x51, 0xD1, 0x31, 0xB1, 0x71, 0xF1,
  0x09, 0x89, 0x49, 0xC9, 0x29, 0xA9, 0x69, 0xE9, 0x19, 0x99, 0x59, 0xD9, 0x39, 0xB9, 0x79, 0xF9,
  0x05, 0x85, 0x45, 0xC5, 0x25, 0xA5, 0x65, 0xE5, 0x15, 0x95, 0x55, 0xD5, 0x35, 0xB5, 0x75, 0xF5,
  0x0D, 0x8D, 0x4D, 0xCD, 0x2D, 0xAD, 0x6D, 0xED, 0x1D, 0x9D, 0x5D, 0xDD, 0x3D, 0xBD, 0x7D, 0xFD,
  0x03, 0x83, 0x43, 0xC3, 0x23, 0xA3, 0x63, 0xE3, 0x13, 0x93, 0x53, 0xD3, 0x33, 0xB3, 0x73, 0xF3,
  0x0B, 0x8B, 0x4B, 0xCB, 0x2B, 0xAB, 0x6B, 0xEB, 0x1B, 0x9B, 0x5B, 0xDB, 0x3B, 0xBB, 0x7B, 0xFB,
  0x07, 0x87, 0x47, 0xC7, 0x27, 0xA7, 0x67, 0xE7, 0x17, 0x97, 0x57, 0xD7, 0x37, 0xB7, 0x77, 0xF7,
  0x0F, 0x8F, 0x4F, 0xCF, 0x2F, 0xAF, 0x6F, 0xEF, 0x1F, 0x9F, 0x5F, 0xDF, 0x3F, 0xBF, 0x7F, 0xFF};

//============================================================================

int nFFT_Valid_Length(FFT_Proc_Rec* this_FFT_plan)
{
    if (((this_FFT_plan->NumberOfBitsNeeded <= nFFT_MAX_FFT_BIT_LENGTH) && (this_FFT_plan->NumberOfBitsNeeded > 0)) && (this_FFT_plan->DComplex != nullptr))
    {
        this_FFT_plan->FFT = &FFT_PreCalc_Data_Rec[this_FFT_plan->NumberOfBitsNeeded];
        this_FFT_plan->BlockBitLen = 0;

        return this_FFT_plan->FFT->length;
    }
    else
        return 0;
}


//=====================================================================================================================
// Radix 2 FFT
//=====================================================================================================================
void Radix_02FS_DIT(FFT_Proc_Rec* this_FFT_plan)
{
    this_FFT_plan->BlockBitLen++;

    tDComplex V1 = this_FFT_plan->DComplex[1];
    this_FFT_plan->DComplex[1] = this_FFT_plan->DComplex[0] - V1;
    this_FFT_plan->DComplex[0] += V1;
}


void Radix_02RS_DIT(FFT_Proc_Rec* this_FFT_plan)
{
    this_FFT_plan->BlockBitLen++;

    tDComplex V1 = this_FFT_plan->DComplex[1];
    this_FFT_plan->DComplex[1] = this_FFT_plan->DComplex[0] - V1;
    this_FFT_plan->DComplex[0] += V1;
}


void Radix_02F_DIT(FFT_Proc_Rec* this_FFT_plan)
{
    int32_t DataStride = 1 << this_FFT_plan->BlockBitLen;

    this_FFT_plan->BlockBitLen++;

    int32_t BlockLength = 1 << this_FFT_plan->BlockBitLen;

    if (this_FFT_plan->BlockBitLen > nFFT_MAX_FFT_BIT_LENGTH)
        throw(-1);

    tDComplex* this_A_Arr = A_Arr[this_FFT_plan->BlockBitLen];

    for (int32_t i = 0; i < this_FFT_plan->FFT->length; i += BlockLength)
    {
        int32_t D0 = i;
        int32_t D1 = i + DataStride;

        tDComplex V1 = this_FFT_plan->DComplex[D1];
        this_FFT_plan->DComplex[D1++] = this_FFT_plan->DComplex[D0] - V1;
        this_FFT_plan->DComplex[D0++] += V1;

        for (int32_t n = 1; n < DataStride; n++)
        {
            tDComplex V1 = (this_FFT_plan->DComplex[D1] * this_A_Arr[n]);
            this_FFT_plan->DComplex[D1++] = this_FFT_plan->DComplex[D0] - V1;
            this_FFT_plan->DComplex[D0++] += V1;
        }
    }
}


void Radix_02R_DIT(FFT_Proc_Rec* this_FFT_plan)
{
    int32_t DataStride = 1 << this_FFT_plan->BlockBitLen;

    this_FFT_plan->BlockBitLen++;

    int32_t BlockLength = 1 << this_FFT_plan->BlockBitLen;

    if (this_FFT_plan->BlockBitLen > nFFT_MAX_FFT_BIT_LENGTH)
        throw(-1);

    tDComplex* this_A_Arr = A_Arr_Conj[this_FFT_plan->BlockBitLen];

    for (int32_t i = 0; i < this_FFT_plan->FFT->length; i += BlockLength)
    {
        int32_t D0 = i;
        int32_t D1 = i + DataStride;

        tDComplex V1 = this_FFT_plan->DComplex[D1];
        this_FFT_plan->DComplex[D1++] = this_FFT_plan->DComplex[D0] - V1;
        this_FFT_plan->DComplex[D0++] += V1;

        for (int32_t n = 1; n < DataStride; n++)
        {
            tDComplex V1 = (this_FFT_plan->DComplex[D1] * this_A_Arr[n]);
            this_FFT_plan->DComplex[D1++] = this_FFT_plan->DComplex[D0] - V1;
            this_FFT_plan->DComplex[D0++] += V1;
        }
    }
}


//=====================================================================================================================
// Radix 4 FFT
//=====================================================================================================================
void Radix_04FS_DIT(FFT_Proc_Rec* this_FFT_plan)
{
    this_FFT_plan->BlockBitLen += 2;

    tDComplex Y0 = (this_FFT_plan->DComplex[0] + this_FFT_plan->DComplex[2]);
    tDComplex Y1 = (this_FFT_plan->DComplex[1] + this_FFT_plan->DComplex[3]);
    tDComplex Y2 = (this_FFT_plan->DComplex[0] - this_FFT_plan->DComplex[2]);
    tDComplex Y3 = (this_FFT_plan->DComplex[1] - this_FFT_plan->DComplex[3]).divided_by_i();

    this_FFT_plan->DComplex[0] = (Y0 + Y1);
    this_FFT_plan->DComplex[1] = (Y2 + Y3);
    this_FFT_plan->DComplex[2] = (Y0 - Y1);
    this_FFT_plan->DComplex[3] = (Y2 - Y3);
}


void Radix_04RS_DIT(FFT_Proc_Rec* this_FFT_plan)
{
    this_FFT_plan->BlockBitLen += 2;

    tDComplex Y0 = (this_FFT_plan->DComplex[0] + this_FFT_plan->DComplex[2]);
    tDComplex Y1 = (this_FFT_plan->DComplex[1] + this_FFT_plan->DComplex[3]);
    tDComplex Y2 = (this_FFT_plan->DComplex[0] - this_FFT_plan->DComplex[2]);
    tDComplex Y3 = (this_FFT_plan->DComplex[1] - this_FFT_plan->DComplex[3]).multiplied_by_i();

    this_FFT_plan->DComplex[0] = (Y0 + Y1);
    this_FFT_plan->DComplex[1] = (Y2 + Y3);
    this_FFT_plan->DComplex[2] = (Y0 - Y1);
    this_FFT_plan->DComplex[3] = (Y2 - Y3);
}


void Radix_04F_DIT(FFT_Proc_Rec* this_FFT_plan)
{
    int32_t DataStride = 1 << this_FFT_plan->BlockBitLen;
    int32_t DataStride_shl_1 = DataStride << 1;

    this_FFT_plan->BlockBitLen += 2;

    int32_t BlockLength = 1 << this_FFT_plan->BlockBitLen;

    if (this_FFT_plan->BlockBitLen > nFFT_MAX_FFT_BIT_LENGTH)
        throw(-1);

    tDComplex* this_A_Arr = A_Arr[this_FFT_plan->BlockBitLen];

    for (int32_t i = 0; i < this_FFT_plan->FFT->length; i += BlockLength)
    {
        int32_t D0 = i;
        int32_t D1 = D0 + DataStride;
        int32_t D2 = D0 + DataStride_shl_1;
        int32_t D3 = D1 + DataStride_shl_1;

        tDComplex V0 = this_FFT_plan->DComplex[D0];
        tDComplex V1 = this_FFT_plan->DComplex[D2];
        tDComplex V2 = this_FFT_plan->DComplex[D1];
        tDComplex V3 = this_FFT_plan->DComplex[D3];

        tDComplex Y0 = (V0 + V2);
        tDComplex Y1 = (V1 + V3);
        tDComplex Y2 = (V0 - V2);
        tDComplex Y3 = (V1 - V3).divided_by_i();

        this_FFT_plan->DComplex[D0++] = (Y0 + Y1);
        this_FFT_plan->DComplex[D1++] = (Y2 + Y3);
        this_FFT_plan->DComplex[D2++] = (Y0 - Y1);
        this_FFT_plan->DComplex[D3++] = (Y2 - Y3);

        for (int32_t n = 1; n < DataStride; n++)
        {
            int32_t P = n;

            tDComplex V0 = this_FFT_plan->DComplex[D0];
            tDComplex V1 = (this_FFT_plan->DComplex[D2] * this_A_Arr[P]);
            tDComplex V2 = (this_FFT_plan->DComplex[D1] * this_A_Arr[P+=n]);
            tDComplex V3 = (this_FFT_plan->DComplex[D3] * this_A_Arr[P+=n]);

            tDComplex Y0 = (V0 + V2);
            tDComplex Y1 = (V1 + V3);
            tDComplex Y2 = (V0 - V2);
            tDComplex Y3 = (V1 - V3).divided_by_i();

            this_FFT_plan->DComplex[D0++] = (Y0 + Y1);
            this_FFT_plan->DComplex[D1++] = (Y2 + Y3);
            this_FFT_plan->DComplex[D2++] = (Y0 - Y1);
            this_FFT_plan->DComplex[D3++] = (Y2 - Y3);
        }
    }
}


void Radix_04R_DIT(FFT_Proc_Rec* this_FFT_plan)
{
    int32_t DataStride = 1 << this_FFT_plan->BlockBitLen;
    int32_t DataStride_shl_1 = DataStride << 1;

    this_FFT_plan->BlockBitLen += 2;

    int32_t BlockLength = 1 << this_FFT_plan->BlockBitLen;

    if (this_FFT_plan->BlockBitLen > nFFT_MAX_FFT_BIT_LENGTH)
        throw(-1);

    tDComplex* this_A_Arr = A_Arr_Conj[this_FFT_plan->BlockBitLen];

    for (int32_t i = 0; i < this_FFT_plan->FFT->length; i += BlockLength)
    {
        int32_t D0 = i;
        int32_t D1 = D0 + DataStride;
        int32_t D2 = D0 + DataStride_shl_1;
        int32_t D3 = D1 + DataStride_shl_1;

        tDComplex V0 = this_FFT_plan->DComplex[D0];
        tDComplex V1 = this_FFT_plan->DComplex[D2];
        tDComplex V2 = this_FFT_plan->DComplex[D1];
        tDComplex V3 = this_FFT_plan->DComplex[D3];

        tDComplex Y0 = (V0 + V2);
        tDComplex Y1 = (V1 + V3);
        tDComplex Y2 = (V0 - V2);
        tDComplex Y3 = (V1 - V3).multiplied_by_i();

        this_FFT_plan->DComplex[D0++] = (Y0 + Y1);
        this_FFT_plan->DComplex[D1++] = (Y2 + Y3);
        this_FFT_plan->DComplex[D2++] = (Y0 - Y1);
        this_FFT_plan->DComplex[D3++] = (Y2 - Y3);

        for (int32_t n = 1; n < DataStride; n++)
        {
            int32_t P = n;

            tDComplex V0 =  this_FFT_plan->DComplex[D0];
            tDComplex V1 = (this_FFT_plan->DComplex[D2] * this_A_Arr[P]);
            tDComplex V2 = (this_FFT_plan->DComplex[D1] * this_A_Arr[P+=n]);
            tDComplex V3 = (this_FFT_plan->DComplex[D3] * this_A_Arr[P+=n]);

            tDComplex Y0 = (V0 + V2);
            tDComplex Y1 = (V1 + V3);
            tDComplex Y2 = (V0 - V2);
            tDComplex Y3 = (V1 - V3).multiplied_by_i();

            this_FFT_plan->DComplex[D0++] = (Y0 + Y1);
            this_FFT_plan->DComplex[D1++] = (Y2 + Y3);
            this_FFT_plan->DComplex[D2++] = (Y0 - Y1);
            this_FFT_plan->DComplex[D3++] = (Y2 - Y3);
        }
    }
}


//=====================================================================================================================
// Radix 8 FFT
//=====================================================================================================================
void Radix_08FS_DIT(FFT_Proc_Rec* this_FFT_plan)
{
    this_FFT_plan->BlockBitLen += 3;

    tDComplex X0 = (this_FFT_plan->DComplex[0] + this_FFT_plan->DComplex[4]);
    tDComplex X1 = (this_FFT_plan->DComplex[1] + this_FFT_plan->DComplex[7]);
    tDComplex X2 = (this_FFT_plan->DComplex[2] + this_FFT_plan->DComplex[6]);
    tDComplex X3 = (this_FFT_plan->DComplex[3] + this_FFT_plan->DComplex[5]);
    tDComplex X4 = (this_FFT_plan->DComplex[0] - this_FFT_plan->DComplex[4]);
    tDComplex X5 = (this_FFT_plan->DComplex[1] - this_FFT_plan->DComplex[7]);
    tDComplex X6 = (this_FFT_plan->DComplex[2] - this_FFT_plan->DComplex[6]);
    tDComplex X7 = (this_FFT_plan->DComplex[3] - this_FFT_plan->DComplex[5]);

    tDComplex Y0 = (X0 + X2);
    tDComplex Y1 = (X1 + X3);
    tDComplex Y2 = (X1 - X3) * cos_pi_over_32[8];
    tDComplex Y3 = (X5 + X7) * cos_pi_over_32[8];

//    tDComplex Z0 = (Y0 + Y1);
    tDComplex Z1 = (X4 + Y2);
    tDComplex Z2 = (X0 - X2);
    tDComplex Z3 = (X4 - Y2);
//    tDComplex Z4 = (Y0 - Y1);
    tDComplex Z5 = (Y3 + X6).divided_by_i();
    tDComplex Z6 = (X5 - X7).divided_by_i();
    tDComplex Z7 = (Y3 - X6).divided_by_i();

    this_FFT_plan->DComplex[0] = (Y0 + Y1); // Z0;
    this_FFT_plan->DComplex[1] = (Z1 + Z5);
    this_FFT_plan->DComplex[2] = (Z2 + Z6);
    this_FFT_plan->DComplex[3] = (Z3 + Z7);
    this_FFT_plan->DComplex[4] = (Y0 - Y1); // Z4;
    this_FFT_plan->DComplex[5] = (Z3 - Z7);
    this_FFT_plan->DComplex[6] = (Z2 - Z6);
    this_FFT_plan->DComplex[7] = (Z1 - Z5);
}


void Radix_08RS_DIT(FFT_Proc_Rec* this_FFT_plan)
{
    this_FFT_plan->BlockBitLen += 3;

    tDComplex X0 = (this_FFT_plan->DComplex[0] + this_FFT_plan->DComplex[4]);
    tDComplex X1 = (this_FFT_plan->DComplex[1] + this_FFT_plan->DComplex[7]);
    tDComplex X2 = (this_FFT_plan->DComplex[2] + this_FFT_plan->DComplex[6]);
    tDComplex X3 = (this_FFT_plan->DComplex[3] + this_FFT_plan->DComplex[5]);
    tDComplex X4 = (this_FFT_plan->DComplex[0] - this_FFT_plan->DComplex[4]);
    tDComplex X5 = (this_FFT_plan->DComplex[1] - this_FFT_plan->DComplex[7]);
    tDComplex X6 = (this_FFT_plan->DComplex[2] - this_FFT_plan->DComplex[6]);
    tDComplex X7 = (this_FFT_plan->DComplex[3] - this_FFT_plan->DComplex[5]);

    tDComplex Y0 = (X0 + X2);
    tDComplex Y1 = (X1 + X3);
    tDComplex Y2 = (X1 - X3) * cos_pi_over_32[8];
    tDComplex Y3 = (X5 + X7) * cos_pi_over_32[8];

//    tDComplex Z0 = (Y0 + Y1);
    tDComplex Z1 = (X4 + Y2);
    tDComplex Z2 = (X0 - X2);
    tDComplex Z3 = (X4 - Y2);
//    tDComplex Z4 = (Y0 - Y1);
    tDComplex Z5 = (Y3 + X6).multiplied_by_i();
    tDComplex Z6 = (X5 - X7).multiplied_by_i();
    tDComplex Z7 = (Y3 - X6).multiplied_by_i();

    this_FFT_plan->DComplex[0] = (Y0 + Y1); // Z0;
    this_FFT_plan->DComplex[1] = (Z1 + Z5);
    this_FFT_plan->DComplex[2] = (Z2 + Z6);
    this_FFT_plan->DComplex[3] = (Z3 + Z7);
    this_FFT_plan->DComplex[4] = (Y0 - Y1); // Z4;
    this_FFT_plan->DComplex[5] = (Z3 - Z7);
    this_FFT_plan->DComplex[6] = (Z2 - Z6);
    this_FFT_plan->DComplex[7] = (Z1 - Z5);
}


void Radix_08F_DIT(FFT_Proc_Rec* this_FFT_plan)
{
    int32_t DataStride = 1 << this_FFT_plan->BlockBitLen;
    int32_t DataStride_shl_1 = DataStride << 1;
    int32_t DataStride_shl_2 = DataStride << 2;

    this_FFT_plan->BlockBitLen += 3;

    int32_t BlockLength = 1 << this_FFT_plan->BlockBitLen;

    if (this_FFT_plan->BlockBitLen > nFFT_MAX_FFT_BIT_LENGTH)
        throw(-1);

    tDComplex* this_A_Arr = A_Arr[this_FFT_plan->BlockBitLen];

    for (int32_t i = 0; i < this_FFT_plan->FFT->length; i += BlockLength)
    {
        int32_t D0 = i;
        int32_t D1 = D0 + DataStride;
        int32_t D2 = D0 + DataStride_shl_1;
        int32_t D3 = D1 + DataStride_shl_1;
        int32_t D4 = D0 + DataStride_shl_2;
        int32_t D5 = D1 + DataStride_shl_2;
        int32_t D6 = D2 + DataStride_shl_2;
        int32_t D7 = D3 + DataStride_shl_2;

        tDComplex V0 = this_FFT_plan->DComplex[D0];
        tDComplex V1 = this_FFT_plan->DComplex[D4];
        tDComplex V2 = this_FFT_plan->DComplex[D2];
        tDComplex V3 = this_FFT_plan->DComplex[D6];
        tDComplex V4 = this_FFT_plan->DComplex[D1];
        tDComplex V5 = this_FFT_plan->DComplex[D5];
        tDComplex V6 = this_FFT_plan->DComplex[D3];
        tDComplex V7 = this_FFT_plan->DComplex[D7];

        tDComplex X0 = (V0 + V4);
        tDComplex X1 = (V1 + V7);
        tDComplex X2 = (V2 + V6);
        tDComplex X3 = (V3 + V5);
        tDComplex X4 = (V0 - V4);
        tDComplex X5 = (V1 - V7);
        tDComplex X6 = (V2 - V6);
        tDComplex X7 = (V3 - V5);

        tDComplex Y0 = (X0 + X2);
        tDComplex Y1 = (X1 + X3);
        tDComplex Y2 = (X1 - X3) * cos_pi_over_32[8];
        tDComplex Y3 = (X5 + X7) * cos_pi_over_32[8];

    //    tDComplex Z0 = (Y0 + Y1);
        tDComplex Z1 = (X4 + Y2);
        tDComplex Z2 = (X0 - X2);
        tDComplex Z3 = (X4 - Y2);
    //    tDComplex Z4 = (Y0 - Y1);
        tDComplex Z5 = (Y3 + X6).divided_by_i();
        tDComplex Z6 = (X5 - X7).divided_by_i();
        tDComplex Z7 = (Y3 - X6).divided_by_i();

        this_FFT_plan->DComplex[D0++] = (Y0 + Y1); // Z0;
        this_FFT_plan->DComplex[D1++] = (Z1 + Z5);
        this_FFT_plan->DComplex[D2++] = (Z2 + Z6);
        this_FFT_plan->DComplex[D3++] = (Z3 + Z7);
        this_FFT_plan->DComplex[D4++] = (Y0 - Y1); // Z4;
        this_FFT_plan->DComplex[D5++] = (Z3 - Z7);
        this_FFT_plan->DComplex[D6++] = (Z2 - Z6);
        this_FFT_plan->DComplex[D7++] = (Z1 - Z5);

        for (int32_t n = 1; n < DataStride; n++)
        {
            int32_t P=n;

            tDComplex V0 = this_FFT_plan->DComplex[D0];
            tDComplex V1 = (this_FFT_plan->DComplex[D4] * this_A_Arr[P]);
            tDComplex V2 = (this_FFT_plan->DComplex[D2] * this_A_Arr[P+=n]);
            tDComplex V3 = (this_FFT_plan->DComplex[D6] * this_A_Arr[P+=n]);
            tDComplex V4 = (this_FFT_plan->DComplex[D1] * this_A_Arr[P+=n]);
            tDComplex V5 = (this_FFT_plan->DComplex[D5] * this_A_Arr[P+=n]);
            tDComplex V6 = (this_FFT_plan->DComplex[D3] * this_A_Arr[P+=n]);
            tDComplex V7 = (this_FFT_plan->DComplex[D7] * this_A_Arr[P+=n]);

            tDComplex X0 = (V0 + V4);
            tDComplex X1 = (V1 + V7);
            tDComplex X2 = (V2 + V6);
            tDComplex X3 = (V3 + V5);
            tDComplex X4 = (V0 - V4);
            tDComplex X5 = (V1 - V7);
            tDComplex X6 = (V2 - V6);
            tDComplex X7 = (V3 - V5);

            tDComplex Y0 = (X0 + X2);
            tDComplex Y1 = (X1 + X3);
            tDComplex Y2 = (X1 - X3) * cos_pi_over_32[8];
            tDComplex Y3 = (X5 + X7) * cos_pi_over_32[8];

        //    tDComplex Z0 = (Y0 + Y1);
            tDComplex Z1 = (X4 + Y2);
            tDComplex Z2 = (X0 - X2);
            tDComplex Z3 = (X4 - Y2);
        //    tDComplex Z4 = (Y0 - Y1);
            tDComplex Z5 = (Y3 + X6).divided_by_i();
            tDComplex Z6 = (X5 - X7).divided_by_i();
            tDComplex Z7 = (Y3 - X6).divided_by_i();

            this_FFT_plan->DComplex[D0++] = (Y0 + Y1); // Z0;
            this_FFT_plan->DComplex[D1++] = (Z1 + Z5);
            this_FFT_plan->DComplex[D2++] = (Z2 + Z6);
            this_FFT_plan->DComplex[D3++] = (Z3 + Z7);
            this_FFT_plan->DComplex[D4++] = (Y0 - Y1); // Z4;
            this_FFT_plan->DComplex[D5++] = (Z3 - Z7);
            this_FFT_plan->DComplex[D6++] = (Z2 - Z6);
            this_FFT_plan->DComplex[D7++] = (Z1 - Z5);
        }
    }
}


void Radix_08R_DIT(FFT_Proc_Rec* this_FFT_plan)
{
    int32_t DataStride = 1 << this_FFT_plan->BlockBitLen;
    int32_t DataStride_shl_1 = DataStride << 1;
    int32_t DataStride_shl_2 = DataStride << 2;

    this_FFT_plan->BlockBitLen += 3;

    int32_t BlockLength = 1 << this_FFT_plan->BlockBitLen;

    if (this_FFT_plan->BlockBitLen > nFFT_MAX_FFT_BIT_LENGTH)
        throw(-1);

    tDComplex* this_A_Arr = A_Arr_Conj[this_FFT_plan->BlockBitLen];

    for (int32_t i = 0; i < this_FFT_plan->FFT->length; i += BlockLength)
    {
        int32_t D0 = i;
        int32_t D1 = D0 + DataStride;
        int32_t D2 = D0 + DataStride_shl_1;
        int32_t D3 = D1 + DataStride_shl_1;
        int32_t D4 = D0 + DataStride_shl_2;
        int32_t D5 = D1 + DataStride_shl_2;
        int32_t D6 = D2 + DataStride_shl_2;
        int32_t D7 = D3 + DataStride_shl_2;

        tDComplex V0 = this_FFT_plan->DComplex[D0];
        tDComplex V1 = this_FFT_plan->DComplex[D4];
        tDComplex V2 = this_FFT_plan->DComplex[D2];
        tDComplex V3 = this_FFT_plan->DComplex[D6];
        tDComplex V4 = this_FFT_plan->DComplex[D1];
        tDComplex V5 = this_FFT_plan->DComplex[D5];
        tDComplex V6 = this_FFT_plan->DComplex[D3];
        tDComplex V7 = this_FFT_plan->DComplex[D7];

        tDComplex X0 = (V0 + V4);
        tDComplex X1 = (V1 + V7);
        tDComplex X2 = (V2 + V6);
        tDComplex X3 = (V3 + V5);
        tDComplex X4 = (V0 - V4);
        tDComplex X5 = (V1 - V7);
        tDComplex X6 = (V2 - V6);
        tDComplex X7 = (V3 - V5);

        tDComplex Y0 = (X0 + X2);
        tDComplex Y1 = (X1 + X3);
        tDComplex Y2 = (X1 - X3) * cos_pi_over_32[8];
        tDComplex Y3 = (X5 + X7) * cos_pi_over_32[8];

        tDComplex Z1 = (X4 + Y2);
        tDComplex Z2 = (X0 - X2);
        tDComplex Z3 = (X4 - Y2);

        tDComplex Z5 = (Y3 + X6).multiplied_by_i();
        tDComplex Z6 = (X5 - X7).multiplied_by_i();
        tDComplex Z7 = (Y3 - X6).multiplied_by_i();

        this_FFT_plan->DComplex[D0++] = (Y0 + Y1);
        this_FFT_plan->DComplex[D1++] = (Z1 + Z5);
        this_FFT_plan->DComplex[D2++] = (Z2 + Z6);
        this_FFT_plan->DComplex[D3++] = (Z3 + Z7);
        this_FFT_plan->DComplex[D4++] = (Y0 - Y1);
        this_FFT_plan->DComplex[D5++] = (Z3 - Z7);
        this_FFT_plan->DComplex[D6++] = (Z2 - Z6);
        this_FFT_plan->DComplex[D7++] = (Z1 - Z5);

        for (int32_t n = 1; n < DataStride; n++)
        {
            int32_t P=n;

            V0 = this_FFT_plan->DComplex[D0];
            V1 = (this_FFT_plan->DComplex[D4] * this_A_Arr[P]);
            V2 = (this_FFT_plan->DComplex[D2] * this_A_Arr[P+=n]);
            V3 = (this_FFT_plan->DComplex[D6] * this_A_Arr[P+=n]);
            V4 = (this_FFT_plan->DComplex[D1] * this_A_Arr[P+=n]);
            V5 = (this_FFT_plan->DComplex[D5] * this_A_Arr[P+=n]);
            V6 = (this_FFT_plan->DComplex[D3] * this_A_Arr[P+=n]);
            V7 = (this_FFT_plan->DComplex[D7] * this_A_Arr[P+=n]);

            tDComplex X0 = (V0 + V4);
            tDComplex X1 = (V1 + V7);
            tDComplex X2 = (V2 + V6);
            tDComplex X3 = (V3 + V5);
            tDComplex X4 = (V0 - V4);
            tDComplex X5 = (V1 - V7);
            tDComplex X6 = (V2 - V6);
            tDComplex X7 = (V3 - V5);

            tDComplex Y0 = (X0 + X2);
            tDComplex Y1 = (X1 + X3);
            tDComplex Y2 = (X1 - X3) * cos_pi_over_32[8];
            tDComplex Y3 = (X5 + X7) * cos_pi_over_32[8];

            tDComplex Z1 = (X4 + Y2);
            tDComplex Z2 = (X0 - X2);
            tDComplex Z3 = (X4 - Y2);
            tDComplex Z5 = (Y3 + X6).multiplied_by_i();
            tDComplex Z6 = (X5 - X7).multiplied_by_i();
            tDComplex Z7 = (Y3 - X6).multiplied_by_i();

            this_FFT_plan->DComplex[D0++] = (Y0 + Y1);
            this_FFT_plan->DComplex[D1++] = (Z1 + Z5);
            this_FFT_plan->DComplex[D2++] = (Z2 + Z6);
            this_FFT_plan->DComplex[D3++] = (Z3 + Z7);
            this_FFT_plan->DComplex[D4++] = (Y0 - Y1);
            this_FFT_plan->DComplex[D5++] = (Z3 - Z7);
            this_FFT_plan->DComplex[D6++] = (Z2 - Z6);
            this_FFT_plan->DComplex[D7++] = (Z1 - Z5);
        }
    }
}


//=====================================================================================================================
// Radix 16 FFT
//=====================================================================================================================
void Radix_16FS_DIT(FFT_Proc_Rec* this_FFT_plan)
{
    this_FFT_plan->BlockBitLen += 4;

    tDComplex x08  = (this_FFT_plan->DComplex[0] + this_FFT_plan->DComplex[8]);
    tDComplex x0_8 = (this_FFT_plan->DComplex[0] - this_FFT_plan->DComplex[8]);

    tDComplex x1f  = (this_FFT_plan->DComplex[1] + this_FFT_plan->DComplex[15]);
    tDComplex x1_f = (this_FFT_plan->DComplex[1] - this_FFT_plan->DComplex[15]);

    tDComplex x2e  = (this_FFT_plan->DComplex[2] + this_FFT_plan->DComplex[14]);
    tDComplex x2_e = (this_FFT_plan->DComplex[2] - this_FFT_plan->DComplex[14]);

    tDComplex x3d  = (this_FFT_plan->DComplex[3] + this_FFT_plan->DComplex[13]);
    tDComplex x3_d = (this_FFT_plan->DComplex[3] - this_FFT_plan->DComplex[13]);

    tDComplex x4c  = (this_FFT_plan->DComplex[4] + this_FFT_plan->DComplex[12]);
    tDComplex x4_c = (this_FFT_plan->DComplex[4] - this_FFT_plan->DComplex[12]);

    tDComplex xb5  = (this_FFT_plan->DComplex[11] + this_FFT_plan->DComplex[5]);
    tDComplex xb_5 = (this_FFT_plan->DComplex[11] - this_FFT_plan->DComplex[5]);

    tDComplex x6a  = (this_FFT_plan->DComplex[6] + this_FFT_plan->DComplex[10]);
    tDComplex x6_a = (this_FFT_plan->DComplex[6] - this_FFT_plan->DComplex[10]);

    tDComplex x97  = (this_FFT_plan->DComplex[9] + this_FFT_plan->DComplex[7]);
    tDComplex x9_7 = (this_FFT_plan->DComplex[9] - this_FFT_plan->DComplex[7]);

    tDComplex x1f97  = (x1f + x97);
    tDComplex x1f_97 = (x1f - x97);

    tDComplex x3db5  = (x3d + xb5);
    tDComplex x3d_b5 = (x3d - xb5);

    tDComplex x2e6a  =  (x2e + x6a);
    tDComplex x2e_6a =  (x2e - x6a);

    tDComplex y91_b3 = (x2e_6a * cos_pi_over_32[8]);

    tDComplex y91 = (x0_8 + y91_b3);
    tDComplex yb3 = (x0_8 - y91_b3);

    tDComplex y3_b = ((x1f_97 * cos_pi_over_32[12]) - (x3d_b5 * cos_pi_over_32[4]));
    tDComplex y1_9 = ((x1f_97 * cos_pi_over_32[4]) + (x3d_b5 * cos_pi_over_32[12]));

    tDComplex y1 = (y91 + y1_9);
    tDComplex y9 = (y91 - y1_9);

    tDComplex y3 = (yb3 + y3_b);
    tDComplex yb = (yb3 - y3_b);

    tDComplex y62 = (x08 - x4c);

    tDComplex y2_6 = ((x1f97 - x3db5)* cos_pi_over_32[8]);

    tDComplex y2 = (y62 + y2_6);
    tDComplex y6 = (y62 - y2_6);

    tDComplex y048 = (x08 + x4c);

    tDComplex y4  = (y048 - x2e6a);
    tDComplex y08 = (y048 + x2e6a);

    tDComplex y0_8 = (x1f97 + x3db5);

//    tDComplex y0 = (y08 + y0_8);
//    tDComplex y8 = (y08 - y0_8);

    tDComplex x19_f7 = (x1_f + x9_7);
    tDComplex x17_f9 = (x1_f - x9_7);

    tDComplex x3b_d5 = (x3_d + xb_5);
    tDComplex x35_db = (x3_d - xb_5);

    tDComplex x26_ea = (x2_e + x6_a);
    tDComplex x2a_e6 = (x2_e - x6_a);

    tDComplex zb931 = (x26_ea * cos_pi_over_32[8]);

    tDComplex zb3 = (zb931 - x4_c);
    tDComplex z91 = (zb931 + x4_c);

    tDComplex z3_b = ((x17_f9 * cos_pi_over_32[4]) - (x35_db * cos_pi_over_32[12]));
    tDComplex z1_9 = ((x17_f9 * cos_pi_over_32[12]) + (x35_db * cos_pi_over_32[4]));

    tDComplex z62 = ((x19_f7 + x3b_d5) * cos_pi_over_32[8]);

    tDComplex z1 = (z91 + z1_9).divided_by_i();
    tDComplex z2 = (z62 + x2a_e6).divided_by_i();
    tDComplex z3 = (zb3 + z3_b).divided_by_i();
    tDComplex z4 = (x19_f7 - x3b_d5).divided_by_i();
    tDComplex z6 = (z62 - x2a_e6).divided_by_i();
    tDComplex z9 = (z91 - z1_9).divided_by_i();
    tDComplex zb = (zb3 - z3_b).divided_by_i();

    this_FFT_plan->DComplex[0] = (y08 + y0_8);
    this_FFT_plan->DComplex[1] = (y1 + z1);
    this_FFT_plan->DComplex[2] = (y2 + z2);
    this_FFT_plan->DComplex[3] = (y3 + z3);
    this_FFT_plan->DComplex[4] = (y4 + z4);
    this_FFT_plan->DComplex[5] = (yb - zb);
    this_FFT_plan->DComplex[6] = (y6 + z6);
    this_FFT_plan->DComplex[7] = (y9 - z9);
    this_FFT_plan->DComplex[8] = (y08 - y0_8);
    this_FFT_plan->DComplex[9] = (y9 + z9);
    this_FFT_plan->DComplex[10] = (y6 - z6);
    this_FFT_plan->DComplex[11] = (yb + zb);
    this_FFT_plan->DComplex[12] = (y4 - z4);
    this_FFT_plan->DComplex[13] = (y3 - z3);
    this_FFT_plan->DComplex[14] = (y2 - z2);
    this_FFT_plan->DComplex[15] = (y1 - z1);
}


void Radix_16RS_DIT(FFT_Proc_Rec* this_FFT_plan)
{
    this_FFT_plan->BlockBitLen += 4;

    tDComplex x08  = (this_FFT_plan->DComplex[0] + this_FFT_plan->DComplex[8]);
    tDComplex x0_8 = (this_FFT_plan->DComplex[0] - this_FFT_plan->DComplex[8]);

    tDComplex x1f  = (this_FFT_plan->DComplex[1] + this_FFT_plan->DComplex[15]);
    tDComplex x1_f = (this_FFT_plan->DComplex[1] - this_FFT_plan->DComplex[15]);

    tDComplex x2e  = (this_FFT_plan->DComplex[2] + this_FFT_plan->DComplex[14]);
    tDComplex x2_e = (this_FFT_plan->DComplex[2] - this_FFT_plan->DComplex[14]);

    tDComplex x3d  = (this_FFT_plan->DComplex[3] + this_FFT_plan->DComplex[13]);
    tDComplex x3_d = (this_FFT_plan->DComplex[3] - this_FFT_plan->DComplex[13]);

    tDComplex x4c  = (this_FFT_plan->DComplex[4] + this_FFT_plan->DComplex[12]);
    tDComplex x4_c = (this_FFT_plan->DComplex[4] - this_FFT_plan->DComplex[12]);

    tDComplex xb5  = (this_FFT_plan->DComplex[11] + this_FFT_plan->DComplex[5]);
    tDComplex xb_5 = (this_FFT_plan->DComplex[11] - this_FFT_plan->DComplex[5]);

    tDComplex x6a  = (this_FFT_plan->DComplex[6] + this_FFT_plan->DComplex[10]);
    tDComplex x6_a = (this_FFT_plan->DComplex[6] - this_FFT_plan->DComplex[10]);

    tDComplex x97  = (this_FFT_plan->DComplex[9] + this_FFT_plan->DComplex[7]);
    tDComplex x9_7 = (this_FFT_plan->DComplex[9] - this_FFT_plan->DComplex[7]);

    tDComplex x1f97  = (x1f + x97);
    tDComplex x1f_97 = (x1f - x97);

    tDComplex x3db5  = (x3d + xb5);
    tDComplex x3d_b5 = (x3d - xb5);

    tDComplex x2e6a  =  (x2e + x6a);
    tDComplex x2e_6a =  (x2e - x6a);

    tDComplex y91_b3 = (x2e_6a * cos_pi_over_32[8]);

    tDComplex y91 = (x0_8 + y91_b3);
    tDComplex yb3 = (x0_8 - y91_b3);

    tDComplex y3_b = ((x1f_97 * cos_pi_over_32[12]) - (x3d_b5 * cos_pi_over_32[4]));
    tDComplex y1_9 = ((x1f_97 * cos_pi_over_32[4]) + (x3d_b5 * cos_pi_over_32[12]));

    tDComplex y1 = (y91 + y1_9);
    tDComplex y9 = (y91 - y1_9);

    tDComplex y3 = (yb3 + y3_b);
    tDComplex yb = (yb3 - y3_b);

    tDComplex y62 = (x08 - x4c);

    tDComplex y2_6 = ((x1f97 - x3db5)* cos_pi_over_32[8]);

    tDComplex y2 = (y62 + y2_6);
    tDComplex y6 = (y62 - y2_6);

    tDComplex y048 = (x08 + x4c);

    tDComplex y4  = (y048 - x2e6a);
    tDComplex y08 = (y048 + x2e6a);

    tDComplex y0_8 = (x1f97 + x3db5);

//    tDComplex y0 = (y08 + y0_8);
//    tDComplex y8 = (y08 - y0_8);

    tDComplex x19_f7 = (x1_f + x9_7);
    tDComplex x17_f9 = (x1_f - x9_7);

    tDComplex x3b_d5 = (x3_d + xb_5);
    tDComplex x35_db = (x3_d - xb_5);

    tDComplex x26_ea = (x2_e + x6_a);
    tDComplex x2a_e6 = (x2_e - x6_a);

    tDComplex zb931 = (x26_ea * cos_pi_over_32[8]);

    tDComplex zb3 = (zb931 - x4_c);
    tDComplex z91 = (zb931 + x4_c);

    tDComplex z3_b = ((x17_f9 * cos_pi_over_32[4]) - (x35_db * cos_pi_over_32[12]));
    tDComplex z1_9 = ((x17_f9 * cos_pi_over_32[12]) + (x35_db * cos_pi_over_32[4]));

    tDComplex z62 = ((x19_f7 + x3b_d5) * cos_pi_over_32[8]);

    tDComplex z1 = (z91 + z1_9).multiplied_by_i();
    tDComplex z2 = (z62 + x2a_e6).multiplied_by_i();
    tDComplex z3 = (zb3 + z3_b).multiplied_by_i();
    tDComplex z4 = (x19_f7 - x3b_d5).multiplied_by_i();
    tDComplex z6 = (z62 - x2a_e6).multiplied_by_i();
    tDComplex z9 = (z91 - z1_9).multiplied_by_i();
    tDComplex zb = (zb3 - z3_b).multiplied_by_i();

    this_FFT_plan->DComplex[0] = (y08 + y0_8);
    this_FFT_plan->DComplex[1] = (y1 + z1);
    this_FFT_plan->DComplex[2] = (y2 + z2);
    this_FFT_plan->DComplex[3] = (y3 + z3);
    this_FFT_plan->DComplex[4] = (y4 + z4);
    this_FFT_plan->DComplex[5] = (yb - zb);
    this_FFT_plan->DComplex[6] = (y6 + z6);
    this_FFT_plan->DComplex[7] = (y9 - z9);
    this_FFT_plan->DComplex[8] = (y08 - y0_8);
    this_FFT_plan->DComplex[9] = (y9 + z9);
    this_FFT_plan->DComplex[10] = (y6 - z6);
    this_FFT_plan->DComplex[11] = (yb + zb);
    this_FFT_plan->DComplex[12] = (y4 - z4);
    this_FFT_plan->DComplex[13] = (y3 - z3);
    this_FFT_plan->DComplex[14] = (y2 - z2);
    this_FFT_plan->DComplex[15] = (y1 - z1);
}


void Radix_16F_DIT(FFT_Proc_Rec* this_FFT_plan)
{
    int32_t DataStride = 1 << this_FFT_plan->BlockBitLen;
    int32_t DataStride_shl_1 = DataStride << 1;
    int32_t DataStride_shl_2 = DataStride << 2;
    int32_t DataStride_shl_3 = DataStride << 3;

    this_FFT_plan->BlockBitLen += 4;

    int32_t BlockLength = 1 << this_FFT_plan->BlockBitLen;

    if (this_FFT_plan->BlockBitLen > nFFT_MAX_FFT_BIT_LENGTH)
        throw(-1);

    tDComplex* this_A_Arr = A_Arr[this_FFT_plan->BlockBitLen];

    for (int32_t i = 0; i < this_FFT_plan->FFT->length; i += BlockLength)
    {
        int32_t D0 = i;
        int32_t D1 = D0 + DataStride;
        int32_t D2 = D0 + DataStride_shl_1;
        int32_t D3 = D1 + DataStride_shl_1;
        int32_t D4 = D0 + DataStride_shl_2;
        int32_t D5 = D1 + DataStride_shl_2;
        int32_t D6 = D2 + DataStride_shl_2;
        int32_t D7 = D3 + DataStride_shl_2;
        int32_t D8 = D0 + DataStride_shl_3;
        int32_t D9 = D1 + DataStride_shl_3;
        int32_t Da = D2 + DataStride_shl_3;
        int32_t Db = D3 + DataStride_shl_3;
        int32_t Dc = D4 + DataStride_shl_3;
        int32_t Dd = D5 + DataStride_shl_3;
        int32_t De = D6 + DataStride_shl_3;
        int32_t Df = D7 + DataStride_shl_3;

        tDComplex V0 = this_FFT_plan->DComplex[D0];
        tDComplex V1 = this_FFT_plan->DComplex[D8];
        tDComplex V2 = this_FFT_plan->DComplex[D4];
        tDComplex V3 = this_FFT_plan->DComplex[Dc];
        tDComplex V4 = this_FFT_plan->DComplex[D2];
        tDComplex V5 = this_FFT_plan->DComplex[Da];
        tDComplex V6 = this_FFT_plan->DComplex[D6];
        tDComplex V7 = this_FFT_plan->DComplex[De];
        tDComplex V8 = this_FFT_plan->DComplex[D1];
        tDComplex V9 = this_FFT_plan->DComplex[D9];
        tDComplex Va = this_FFT_plan->DComplex[D5];
        tDComplex Vb = this_FFT_plan->DComplex[Dd];
        tDComplex Vc = this_FFT_plan->DComplex[D3];
        tDComplex Vd = this_FFT_plan->DComplex[Db];
        tDComplex Ve = this_FFT_plan->DComplex[D7];
        tDComplex Vf = this_FFT_plan->DComplex[Df];

        tDComplex x08  = (V0 + V8);
        tDComplex x0_8 = (V0 - V8);

        tDComplex x1f  = (V1 + Vf);
        tDComplex x1_f = (V1 - Vf);

        tDComplex x2e  = (V2 + Ve);
        tDComplex x2_e = (V2 - Ve);

        tDComplex x3d  = (V3 + Vd);
        tDComplex x3_d = (V3 - Vd);

        tDComplex x4c  = (V4 + Vc);
        tDComplex x4_c = (V4 - Vc);

        tDComplex xb5  = (Vb + V5);
        tDComplex xb_5 = (Vb - V5);

        tDComplex x6a  = (V6 + Va);
        tDComplex x6_a = (V6 - Va);

        tDComplex x97  = (V9 + V7);
        tDComplex x9_7 = (V9 - V7);

        tDComplex x1f97  = (x1f + x97);
        tDComplex x1f_97 = (x1f - x97);

        tDComplex x3db5  = (x3d + xb5);
        tDComplex x3d_b5 = (x3d - xb5);

        tDComplex x2e6a  =  (x2e + x6a);
        tDComplex x2e_6a =  (x2e - x6a);

        tDComplex y91_b3 = (x2e_6a * cos_pi_over_32[8]);

        tDComplex y91 = (x0_8 + y91_b3);
        tDComplex yb3 = (x0_8 - y91_b3);

        tDComplex y3_b = ((x1f_97 * cos_pi_over_32[12]) - (x3d_b5 * cos_pi_over_32[4]));
        tDComplex y1_9 = ((x1f_97 * cos_pi_over_32[4]) + (x3d_b5 * cos_pi_over_32[12]));

        tDComplex y1 = (y91 + y1_9);
        tDComplex y9 = (y91 - y1_9);

        tDComplex y3 = (yb3 + y3_b);
        tDComplex yb = (yb3 - y3_b);

        tDComplex y62 = (x08 - x4c);

        tDComplex y2_6 = ((x1f97 - x3db5)* cos_pi_over_32[8]);

        tDComplex y2 = (y62 + y2_6);
        tDComplex y6 = (y62 - y2_6);

        tDComplex y048 = (x08 + x4c);

        tDComplex y4  = (y048 - x2e6a);
        tDComplex y08 = (y048 + x2e6a);

        tDComplex y0_8 = (x1f97 + x3db5);

    //    tDComplex y0 = (y08 + y0_8);
    //    tDComplex y8 = (y08 - y0_8);

        tDComplex x19_f7 = (x1_f + x9_7);
        tDComplex x17_f9 = (x1_f - x9_7);

        tDComplex x3b_d5 = (x3_d + xb_5);
        tDComplex x35_db = (x3_d - xb_5);

        tDComplex x26_ea = (x2_e + x6_a);
        tDComplex x2a_e6 = (x2_e - x6_a);

        tDComplex zb931 = (x26_ea * cos_pi_over_32[8]);

        tDComplex zb3 = (zb931 - x4_c);
        tDComplex z91 = (zb931 + x4_c);

        tDComplex z3_b = ((x17_f9 * cos_pi_over_32[4]) - (x35_db * cos_pi_over_32[12]));
        tDComplex z1_9 = ((x17_f9 * cos_pi_over_32[12]) + (x35_db * cos_pi_over_32[4]));

        tDComplex z62 = ((x19_f7 + x3b_d5) * cos_pi_over_32[8]);

        tDComplex z1 = (z91 + z1_9).divided_by_i();
        tDComplex z2 = (z62 + x2a_e6).divided_by_i();
        tDComplex z3 = (zb3 + z3_b).divided_by_i();
        tDComplex z4 = (x19_f7 - x3b_d5).divided_by_i();
        tDComplex z6 = (z62 - x2a_e6).divided_by_i();
        tDComplex z9 = (z91 - z1_9).divided_by_i();
        tDComplex zb = (zb3 - z3_b).divided_by_i();

        this_FFT_plan->DComplex[D0++] = (y08 + y0_8);
        this_FFT_plan->DComplex[D1++] = (y1 + z1);
        this_FFT_plan->DComplex[D2++] = (y2 + z2);
        this_FFT_plan->DComplex[D3++] = (y3 + z3);
        this_FFT_plan->DComplex[D4++] = (y4 + z4);
        this_FFT_plan->DComplex[D5++] = (yb - zb);
        this_FFT_plan->DComplex[D6++] = (y6 + z6);
        this_FFT_plan->DComplex[D7++] = (y9 - z9);
        this_FFT_plan->DComplex[D8++] = (y08 - y0_8);
        this_FFT_plan->DComplex[D9++] = (y9 + z9);
        this_FFT_plan->DComplex[Da++] = (y6 - z6);
        this_FFT_plan->DComplex[Db++] = (yb + zb);
        this_FFT_plan->DComplex[Dc++] = (y4 - z4);
        this_FFT_plan->DComplex[Dd++] = (y3 - z3);
        this_FFT_plan->DComplex[De++] = (y2 - z2);
        this_FFT_plan->DComplex[Df++] = (y1 - z1);

        for (int32_t n = 1; n < DataStride; n++)
        {
            int32_t P=n;

            V0 = this_FFT_plan->DComplex[D0];
            V1 = this_FFT_plan->DComplex[D8] * this_A_Arr[P];
            V2 = this_FFT_plan->DComplex[D4] * this_A_Arr[P+=n];
            V3 = this_FFT_plan->DComplex[Dc] * this_A_Arr[P+=n];
            V4 = this_FFT_plan->DComplex[D2] * this_A_Arr[P+=n];
            V5 = this_FFT_plan->DComplex[Da] * this_A_Arr[P+=n];
            V6 = this_FFT_plan->DComplex[D6] * this_A_Arr[P+=n];
            V7 = this_FFT_plan->DComplex[De] * this_A_Arr[P+=n];
            V8 = this_FFT_plan->DComplex[D1] * this_A_Arr[P+=n];
            V9 = this_FFT_plan->DComplex[D9] * this_A_Arr[P+=n];
            Va = this_FFT_plan->DComplex[D5] * this_A_Arr[P+=n];
            Vb = this_FFT_plan->DComplex[Dd] * this_A_Arr[P+=n];
            Vc = this_FFT_plan->DComplex[D3] * this_A_Arr[P+=n];
            Vd = this_FFT_plan->DComplex[Db] * this_A_Arr[P+=n];
            Ve = this_FFT_plan->DComplex[D7] * this_A_Arr[P+=n];
            Vf = this_FFT_plan->DComplex[Df] * this_A_Arr[P+=n];

            tDComplex x08  = (V0 + V8);
            tDComplex x0_8 = (V0 - V8);

            tDComplex x1f  = (V1 + Vf);
            tDComplex x1_f = (V1 - Vf);

            tDComplex x2e  = (V2 + Ve);
            tDComplex x2_e = (V2 - Ve);

            tDComplex x3d  = (V3 + Vd);
            tDComplex x3_d = (V3 - Vd);

            tDComplex x4c  = (V4 + Vc);
            tDComplex x4_c = (V4 - Vc);

            tDComplex xb5  = (Vb + V5);
            tDComplex xb_5 = (Vb - V5);

            tDComplex x6a  = (V6 + Va);
            tDComplex x6_a = (V6 - Va);

            tDComplex x97  = (V9 + V7);
            tDComplex x9_7 = (V9 - V7);

            tDComplex x1f97  = (x1f + x97);
            tDComplex x1f_97 = (x1f - x97);

            tDComplex x3db5  = (x3d + xb5);
            tDComplex x3d_b5 = (x3d - xb5);

            tDComplex x2e6a  =  (x2e + x6a);
            tDComplex x2e_6a =  (x2e - x6a);

            tDComplex y91_b3 = (x2e_6a * cos_pi_over_32[8]);

            tDComplex y91 = (x0_8 + y91_b3);
            tDComplex yb3 = (x0_8 - y91_b3);

            tDComplex y3_b = ((x1f_97 * cos_pi_over_32[12]) - (x3d_b5 * cos_pi_over_32[4]));
            tDComplex y1_9 = ((x1f_97 * cos_pi_over_32[4]) + (x3d_b5 * cos_pi_over_32[12]));

            tDComplex y1 = (y91 + y1_9);
            tDComplex y9 = (y91 - y1_9);

            tDComplex y3 = (yb3 + y3_b);
            tDComplex yb = (yb3 - y3_b);

            tDComplex y62 = (x08 - x4c);

            tDComplex y2_6 = ((x1f97 - x3db5)* cos_pi_over_32[8]);

            tDComplex y2 = (y62 + y2_6);
            tDComplex y6 = (y62 - y2_6);

            tDComplex y048 = (x08 + x4c);

            tDComplex y4  = (y048 - x2e6a);
            tDComplex y08 = (y048 + x2e6a);

            tDComplex y0_8 = (x1f97 + x3db5);

        //    tDComplex y0 = (y08 + y0_8);
        //    tDComplex y8 = (y08 - y0_8);

            tDComplex x19_f7 = (x1_f + x9_7);
            tDComplex x17_f9 = (x1_f - x9_7);

            tDComplex x3b_d5 = (x3_d + xb_5);
            tDComplex x35_db = (x3_d - xb_5);

            tDComplex x26_ea = (x2_e + x6_a);
            tDComplex x2a_e6 = (x2_e - x6_a);

            tDComplex zb931 = (x26_ea * cos_pi_over_32[8]);

            tDComplex zb3 = (zb931 - x4_c);
            tDComplex z91 = (zb931 + x4_c);

            tDComplex z3_b = ((x17_f9 * cos_pi_over_32[4]) - (x35_db * cos_pi_over_32[12]));
            tDComplex z1_9 = ((x17_f9 * cos_pi_over_32[12]) + (x35_db * cos_pi_over_32[4]));

            tDComplex z1 = (z91 + z1_9).divided_by_i();
            tDComplex z2 = (z62 + x2a_e6).divided_by_i();
            tDComplex z3 = (zb3 + z3_b).divided_by_i();
            tDComplex z4 = (x19_f7 - x3b_d5).divided_by_i();
            tDComplex z6 = (z62 - x2a_e6).divided_by_i();
            tDComplex z9 = (z91 - z1_9).divided_by_i();
            tDComplex zb = (zb3 - z3_b).divided_by_i();

            this_FFT_plan->DComplex[D0++] = (y08 + y0_8);
            this_FFT_plan->DComplex[D1++] = (y1 + z1);
            this_FFT_plan->DComplex[D2++] = (y2 + z2);
            this_FFT_plan->DComplex[D3++] = (y3 + z3);
            this_FFT_plan->DComplex[D4++] = (y4 + z4);
            this_FFT_plan->DComplex[D5++] = (yb - zb);
            this_FFT_plan->DComplex[D6++] = (y6 + z6);
            this_FFT_plan->DComplex[D7++] = (y9 - z9);
            this_FFT_plan->DComplex[D8++] = (y08 - y0_8);
            this_FFT_plan->DComplex[D9++] = (y9 + z9);
            this_FFT_plan->DComplex[Da++] = (y6 - z6);
            this_FFT_plan->DComplex[Db++] = (yb + zb);
            this_FFT_plan->DComplex[Dc++] = (y4 - z4);
            this_FFT_plan->DComplex[Dd++] = (y3 - z3);
            this_FFT_plan->DComplex[De++] = (y2 - z2);
            this_FFT_plan->DComplex[Df++] = (y1 - z1);
        }
    }
}


void Radix_16R_DIT(FFT_Proc_Rec* this_FFT_plan)
{
    int32_t DataStride = 1 << this_FFT_plan->BlockBitLen;
    int32_t DataStride_shl_1 = DataStride << 1;
    int32_t DataStride_shl_2 = DataStride << 2;
    int32_t DataStride_shl_3 = DataStride << 3;

    this_FFT_plan->BlockBitLen += 4;

    int32_t BlockLength = 1 << this_FFT_plan->BlockBitLen;

    if (this_FFT_plan->BlockBitLen > nFFT_MAX_FFT_BIT_LENGTH)
        throw(-1);

    tDComplex* this_A_Arr = A_Arr_Conj[this_FFT_plan->BlockBitLen];

    for (int32_t i = 0; i < this_FFT_plan->FFT->length; i += BlockLength)
    {
        int32_t D0 = i;
        int32_t D1 = D0 + DataStride;
        int32_t D2 = D0 + DataStride_shl_1;
        int32_t D3 = D1 + DataStride_shl_1;
        int32_t D4 = D0 + DataStride_shl_2;
        int32_t D5 = D1 + DataStride_shl_2;
        int32_t D6 = D2 + DataStride_shl_2;
        int32_t D7 = D3 + DataStride_shl_2;
        int32_t D8 = D0 + DataStride_shl_3;
        int32_t D9 = D1 + DataStride_shl_3;
        int32_t Da = D2 + DataStride_shl_3;
        int32_t Db = D3 + DataStride_shl_3;
        int32_t Dc = D4 + DataStride_shl_3;
        int32_t Dd = D5 + DataStride_shl_3;
        int32_t De = D6 + DataStride_shl_3;
        int32_t Df = D7 + DataStride_shl_3;

        tDComplex V0 = this_FFT_plan->DComplex[D0];
        tDComplex V1 = this_FFT_plan->DComplex[D8];
        tDComplex V2 = this_FFT_plan->DComplex[D4];
        tDComplex V3 = this_FFT_plan->DComplex[Dc];
        tDComplex V4 = this_FFT_plan->DComplex[D2];
        tDComplex V5 = this_FFT_plan->DComplex[Da];
        tDComplex V6 = this_FFT_plan->DComplex[D6];
        tDComplex V7 = this_FFT_plan->DComplex[De];
        tDComplex V8 = this_FFT_plan->DComplex[D1];
        tDComplex V9 = this_FFT_plan->DComplex[D9];
        tDComplex Va = this_FFT_plan->DComplex[D5];
        tDComplex Vb = this_FFT_plan->DComplex[Dd];
        tDComplex Vc = this_FFT_plan->DComplex[D3];
        tDComplex Vd = this_FFT_plan->DComplex[Db];
        tDComplex Ve = this_FFT_plan->DComplex[D7];
        tDComplex Vf = this_FFT_plan->DComplex[Df];

        tDComplex x08  = (V0 + V8);
        tDComplex x0_8 = (V0 - V8);

        tDComplex x1f  = (V1 + Vf);
        tDComplex x1_f = (V1 - Vf);

        tDComplex x2e  = (V2 + Ve);
        tDComplex x2_e = (V2 - Ve);

        tDComplex x3d  = (V3 + Vd);
        tDComplex x3_d = (V3 - Vd);

        tDComplex x4c  = (V4 + Vc);
        tDComplex x4_c = (V4 - Vc);

        tDComplex xb5  = (Vb + V5);
        tDComplex xb_5 = (Vb - V5);

        tDComplex x6a  = (V6 + Va);
        tDComplex x6_a = (V6 - Va);

        tDComplex x97  = (V9 + V7);
        tDComplex x9_7 = (V9 - V7);

        tDComplex x1f97  = (x1f + x97);
        tDComplex x1f_97 = (x1f - x97);

        tDComplex x3db5  = (x3d + xb5);
        tDComplex x3d_b5 = (x3d - xb5);

        tDComplex x2e6a  =  (x2e + x6a);
        tDComplex x2e_6a =  (x2e - x6a);

        tDComplex y91_b3 = (x2e_6a * cos_pi_over_32[8]);

        tDComplex y91 = (x0_8 + y91_b3);
        tDComplex yb3 = (x0_8 - y91_b3);

        tDComplex y3_b = ((x1f_97 * cos_pi_over_32[12]) - (x3d_b5 * cos_pi_over_32[4]));
        tDComplex y1_9 = ((x1f_97 * cos_pi_over_32[4]) + (x3d_b5 * cos_pi_over_32[12]));

        tDComplex y1 = (y91 + y1_9);
        tDComplex y9 = (y91 - y1_9);

        tDComplex y3 = (yb3 + y3_b);
        tDComplex yb = (yb3 - y3_b);

        tDComplex y62 = (x08 - x4c);

        tDComplex y2_6 = ((x1f97 - x3db5)* cos_pi_over_32[8]);

        tDComplex y2 = (y62 + y2_6);
        tDComplex y6 = (y62 - y2_6);

        tDComplex y048 = (x08 + x4c);

        tDComplex y4  = (y048 - x2e6a);
        tDComplex y08 = (y048 + x2e6a);

        tDComplex y0_8 = (x1f97 + x3db5);

    //    tDComplex y0 = (y08 + y0_8);
    //    tDComplex y8 = (y08 - y0_8);

        tDComplex x19_f7 = (x1_f + x9_7);
        tDComplex x17_f9 = (x1_f - x9_7);

        tDComplex x3b_d5 = (x3_d + xb_5);
        tDComplex x35_db = (x3_d - xb_5);

        tDComplex x26_ea = (x2_e + x6_a);
        tDComplex x2a_e6 = (x2_e - x6_a);

        tDComplex zb931 = (x26_ea * cos_pi_over_32[8]);

        tDComplex zb3 = (zb931 - x4_c);
        tDComplex z91 = (zb931 + x4_c);

        tDComplex z3_b = ((x17_f9 * cos_pi_over_32[4]) - (x35_db * cos_pi_over_32[12]));
        tDComplex z1_9 = ((x17_f9 * cos_pi_over_32[12]) + (x35_db * cos_pi_over_32[4]));

        tDComplex z62 = ((x19_f7 + x3b_d5) * cos_pi_over_32[8]);

        tDComplex z1 = (z91 + z1_9).multiplied_by_i();
        tDComplex z2 = (z62 + x2a_e6).multiplied_by_i();
        tDComplex z3 = (zb3 + z3_b).multiplied_by_i();
        tDComplex z4 = (x19_f7 - x3b_d5).multiplied_by_i();
        tDComplex z6 = (z62 - x2a_e6).multiplied_by_i();
        tDComplex z9 = (z91 - z1_9).multiplied_by_i();
        tDComplex zb = (zb3 - z3_b).multiplied_by_i();

        this_FFT_plan->DComplex[D0++] = (y08 + y0_8);
        this_FFT_plan->DComplex[D1++] = (y1 + z1);
        this_FFT_plan->DComplex[D2++] = (y2 + z2);
        this_FFT_plan->DComplex[D3++] = (y3 + z3);
        this_FFT_plan->DComplex[D4++] = (y4 + z4);
        this_FFT_plan->DComplex[D5++] = (yb - zb);
        this_FFT_plan->DComplex[D6++] = (y6 + z6);
        this_FFT_plan->DComplex[D7++] = (y9 - z9);
        this_FFT_plan->DComplex[D8++] = (y08 - y0_8);
        this_FFT_plan->DComplex[D9++] = (y9 + z9);
        this_FFT_plan->DComplex[Da++] = (y6 - z6);
        this_FFT_plan->DComplex[Db++] = (yb + zb);
        this_FFT_plan->DComplex[Dc++] = (y4 - z4);
        this_FFT_plan->DComplex[Dd++] = (y3 - z3);
        this_FFT_plan->DComplex[De++] = (y2 - z2);
        this_FFT_plan->DComplex[Df++] = (y1 - z1);

        for (int32_t n = 1; n < DataStride; n++)
        {
            int32_t P=n;

            V0 = this_FFT_plan->DComplex[D0];
            V1 = this_FFT_plan->DComplex[D8] * this_A_Arr[P];
            V2 = this_FFT_plan->DComplex[D4] * this_A_Arr[P+=n];
            V3 = this_FFT_plan->DComplex[Dc] * this_A_Arr[P+=n];
            V4 = this_FFT_plan->DComplex[D2] * this_A_Arr[P+=n];
            V5 = this_FFT_plan->DComplex[Da] * this_A_Arr[P+=n];
            V6 = this_FFT_plan->DComplex[D6] * this_A_Arr[P+=n];
            V7 = this_FFT_plan->DComplex[De] * this_A_Arr[P+=n];
            V8 = this_FFT_plan->DComplex[D1] * this_A_Arr[P+=n];
            V9 = this_FFT_plan->DComplex[D9] * this_A_Arr[P+=n];
            Va = this_FFT_plan->DComplex[D5] * this_A_Arr[P+=n];
            Vb = this_FFT_plan->DComplex[Dd] * this_A_Arr[P+=n];
            Vc = this_FFT_plan->DComplex[D3] * this_A_Arr[P+=n];
            Vd = this_FFT_plan->DComplex[Db] * this_A_Arr[P+=n];
            Ve = this_FFT_plan->DComplex[D7] * this_A_Arr[P+=n];
            Vf = this_FFT_plan->DComplex[Df] * this_A_Arr[P+=n];

            tDComplex x08  = (V0 + V8);
            tDComplex x0_8 = (V0 - V8);

            tDComplex x1f  = (V1 + Vf);
            tDComplex x1_f = (V1 - Vf);

            tDComplex x2e  = (V2 + Ve);
            tDComplex x2_e = (V2 - Ve);

            tDComplex x3d  = (V3 + Vd);
            tDComplex x3_d = (V3 - Vd);

            tDComplex x4c  = (V4 + Vc);
            tDComplex x4_c = (V4 - Vc);

            tDComplex xb5  = (Vb + V5);
            tDComplex xb_5 = (Vb - V5);

            tDComplex x6a  = (V6 + Va);
            tDComplex x6_a = (V6 - Va);

            tDComplex x97  = (V9 + V7);
            tDComplex x9_7 = (V9 - V7);

            tDComplex x1f97  = (x1f + x97);
            tDComplex x1f_97 = (x1f - x97);

            tDComplex x3db5  = (x3d + xb5);
            tDComplex x3d_b5 = (x3d - xb5);

            tDComplex x2e6a  =  (x2e + x6a);
            tDComplex x2e_6a =  (x2e - x6a);

            tDComplex y91_b3 = (x2e_6a * cos_pi_over_32[8]);

            tDComplex y91 = (x0_8 + y91_b3);
            tDComplex yb3 = (x0_8 - y91_b3);

            tDComplex y3_b = ((x1f_97 * cos_pi_over_32[12]) - (x3d_b5 * cos_pi_over_32[4]));
            tDComplex y1_9 = ((x1f_97 * cos_pi_over_32[4]) + (x3d_b5 * cos_pi_over_32[12]));

            tDComplex y1 = (y91 + y1_9);
            tDComplex y9 = (y91 - y1_9);

            tDComplex y3 = (yb3 + y3_b);
            tDComplex yb = (yb3 - y3_b);

            tDComplex y62 = (x08 - x4c);

            tDComplex y2_6 = ((x1f97 - x3db5)* cos_pi_over_32[8]);

            tDComplex y2 = (y62 + y2_6);
            tDComplex y6 = (y62 - y2_6);

            tDComplex y048 = (x08 + x4c);

            tDComplex y4  = (y048 - x2e6a);
            tDComplex y08 = (y048 + x2e6a);

            tDComplex y0_8 = (x1f97 + x3db5);

        //    tDComplex y0 = (y08 + y0_8);
        //    tDComplex y8 = (y08 - y0_8);

            tDComplex x19_f7 = (x1_f + x9_7);
            tDComplex x17_f9 = (x1_f - x9_7);

            tDComplex x3b_d5 = (x3_d + xb_5);
            tDComplex x35_db = (x3_d - xb_5);

            tDComplex x26_ea = (x2_e + x6_a);
            tDComplex x2a_e6 = (x2_e - x6_a);

            tDComplex zb931 = (x26_ea * cos_pi_over_32[8]);

            tDComplex zb3 = (zb931 - x4_c);
            tDComplex z91 = (zb931 + x4_c);

            tDComplex z3_b = ((x17_f9 * cos_pi_over_32[4]) - (x35_db * cos_pi_over_32[12]));
            tDComplex z1_9 = ((x17_f9 * cos_pi_over_32[12]) + (x35_db * cos_pi_over_32[4]));

            tDComplex z62 = ((x19_f7 + x3b_d5) * cos_pi_over_32[8]);

            tDComplex z1 = (z91 + z1_9).multiplied_by_i();
            tDComplex z2 = (z62 + x2a_e6).multiplied_by_i();
            tDComplex z3 = (zb3 + z3_b).multiplied_by_i();
            tDComplex z4 = (x19_f7 - x3b_d5).multiplied_by_i();
            tDComplex z6 = (z62 - x2a_e6).multiplied_by_i();
            tDComplex z9 = (z91 - z1_9).multiplied_by_i();
            tDComplex zb = (zb3 - z3_b).multiplied_by_i();

            this_FFT_plan->DComplex[D0++] = (y08 + y0_8);
            this_FFT_plan->DComplex[D1++] = (y1 + z1);
            this_FFT_plan->DComplex[D2++] = (y2 + z2);
            this_FFT_plan->DComplex[D3++] = (y3 + z3);
            this_FFT_plan->DComplex[D4++] = (y4 + z4);
            this_FFT_plan->DComplex[D5++] = (yb - zb);
            this_FFT_plan->DComplex[D6++] = (y6 + z6);
            this_FFT_plan->DComplex[D7++] = (y9 - z9);
            this_FFT_plan->DComplex[D8++] = (y08 - y0_8);
            this_FFT_plan->DComplex[D9++] = (y9 + z9);
            this_FFT_plan->DComplex[Da++] = (y6 - z6);
            this_FFT_plan->DComplex[Db++] = (yb + zb);
            this_FFT_plan->DComplex[Dc++] = (y4 - z4);
            this_FFT_plan->DComplex[Dd++] = (y3 - z3);
            this_FFT_plan->DComplex[De++] = (y2 - z2);
            this_FFT_plan->DComplex[Df++] = (y1 - z1);
        }
    }
}


#ifdef FACTOR_2
void FFT_DIT_01(FFT_Proc_Rec* this_FFT)     { Radix_02FS_DIT(this_FFT); }
void FFT_DIT_02(FFT_Proc_Rec* this_FFT)     { Radix_04FS_DIT(this_FFT); }
void FFT_DIT_03(FFT_Proc_Rec* this_FFT)     { Radix_08FS_DIT(this_FFT); }
void FFT_DIT_04(FFT_Proc_Rec* this_FFT)     { Radix_16FS_DIT(this_FFT); }
void FFT_DIT_05(FFT_Proc_Rec* this_FFT)     { Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); }
void FFT_DIT_06(FFT_Proc_Rec* this_FFT)     { Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); }
void FFT_DIT_07(FFT_Proc_Rec* this_FFT)     { Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); }
void FFT_DIT_08(FFT_Proc_Rec* this_FFT)     { Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); }
void FFT_DIT_09(FFT_Proc_Rec* this_FFT)     { Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); }
void FFT_DIT_10(FFT_Proc_Rec* this_FFT)     { Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); }
void FFT_DIT_11(FFT_Proc_Rec* this_FFT)     { Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); }
void FFT_DIT_12(FFT_Proc_Rec* this_FFT)     { Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); }
void FFT_DIT_13(FFT_Proc_Rec* this_FFT)     { Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); Radix_02F_DIT(this_FFT); }
void FFT_DIT_XX(FFT_Proc_Rec* this_FFT)
{
    while (this_FFT->BlockBitLen<this_FFT->NumberOfBitsNeeded)
        Radix_02F_DIT(this_FFT);
}


void IFFT_DIT_01(FFT_Proc_Rec* this_FFT)    { Radix_02R_DIT(this_FFT); }
void IFFT_DIT_02(FFT_Proc_Rec* this_FFT)    { Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); }
void IFFT_DIT_03(FFT_Proc_Rec* this_FFT)    { Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); }
void IFFT_DIT_04(FFT_Proc_Rec* this_FFT)    { Radix_02R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); }
void IFFT_DIT_05(FFT_Proc_Rec* this_FFT)    { Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); }
void IFFT_DIT_06(FFT_Proc_Rec* this_FFT)    { Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); }
void IFFT_DIT_07(FFT_Proc_Rec* this_FFT)    { Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); }
void IFFT_DIT_08(FFT_Proc_Rec* this_FFT)    { Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); }
void IFFT_DIT_09(FFT_Proc_Rec* this_FFT)    { Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); }
void IFFT_DIT_10(FFT_Proc_Rec* this_FFT)    { Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); }
void IFFT_DIT_11(FFT_Proc_Rec* this_FFT)    { Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); }
void IFFT_DIT_12(FFT_Proc_Rec* this_FFT)    { Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); }
void IFFT_DIT_13(FFT_Proc_Rec* this_FFT)    { Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); Radix_02R_DIT(this_FFT); }
void IFFT_DIT_XX(FFT_Proc_Rec* this_FFT)
{
    while (this_FFT->BlockBitLen<(this_FFT->NumberOfBitsNeeded-1))
        Radix_02R_DIT(this_FFT);
}

#endif


#ifdef FACTOR_4_2
void FFT_DIT_01(FFT_Proc_Rec* this_FFT)     { Radix_02FS_DIT(this_FFT); }
void FFT_DIT_02(FFT_Proc_Rec* this_FFT)     { Radix_04FS_DIT(this_FFT); }
void FFT_DIT_03(FFT_Proc_Rec* this_FFT)     { Radix_08FS_DIT(this_FFT); }
void FFT_DIT_04(FFT_Proc_Rec* this_FFT)     { Radix_16FS_DIT(this_FFT); }
void FFT_DIT_05(FFT_Proc_Rec* this_FFT)     { Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_02F_DIT(this_FFT); }
void FFT_DIT_06(FFT_Proc_Rec* this_FFT)     { Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); }
void FFT_DIT_07(FFT_Proc_Rec* this_FFT)     { Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_02F_DIT(this_FFT); }
void FFT_DIT_08(FFT_Proc_Rec* this_FFT)     { Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); }
void FFT_DIT_09(FFT_Proc_Rec* this_FFT)     { Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_02F_DIT(this_FFT); }
void FFT_DIT_10(FFT_Proc_Rec* this_FFT)     { Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); }
void FFT_DIT_11(FFT_Proc_Rec* this_FFT)     { Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_02F_DIT(this_FFT); }
void FFT_DIT_12(FFT_Proc_Rec* this_FFT)     { Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); }
void FFT_DIT_13(FFT_Proc_Rec* this_FFT)     { Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_02F_DIT(this_FFT); }

void IFFT_DIT_01(FFT_Proc_Rec* this_FFT)    { Radix_02RS_DIT(this_FFT); }
void IFFT_DIT_02(FFT_Proc_Rec* this_FFT)    { Radix_04RS_DIT(this_FFT); }
void IFFT_DIT_03(FFT_Proc_Rec* this_FFT)    { Radix_08RS_DIT(this_FFT); }
void IFFT_DIT_04(FFT_Proc_Rec* this_FFT)    { Radix_16RS_DIT(this_FFT); }
void IFFT_DIT_05(FFT_Proc_Rec* this_FFT)    { Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_02R_DIT(this_FFT); }
void IFFT_DIT_06(FFT_Proc_Rec* this_FFT)    { Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); }
void IFFT_DIT_07(FFT_Proc_Rec* this_FFT)    { Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_02R_DIT(this_FFT); }
void IFFT_DIT_08(FFT_Proc_Rec* this_FFT)    { Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); }
void IFFT_DIT_09(FFT_Proc_Rec* this_FFT)    { Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_02R_DIT(this_FFT); }
void IFFT_DIT_10(FFT_Proc_Rec* this_FFT)    { Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); }
void IFFT_DIT_11(FFT_Proc_Rec* this_FFT)    { Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_02R_DIT(this_FFT); }
void IFFT_DIT_12(FFT_Proc_Rec* this_FFT)    { Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); }
void IFFT_DIT_13(FFT_Proc_Rec* this_FFT)    { Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_02R_DIT(this_FFT); }

void FFT_DIT_XX(FFT_Proc_Rec* this_FFT)
{
    while ((this_FFT->NumberOfBitsNeeded-this_FFT->BlockBitLen)>1)
        Radix_04F_DIT(this_FFT);

    if ((this_FFT->NumberOfBitsNeeded-this_FFT->BlockBitLen)>0)
        Radix_02F_DIT(this_FFT);
}

void IFFT_DIT_XX(FFT_Proc_Rec* this_FFT)
{
    while ((this_FFT->NumberOfBitsNeeded-this_FFT->BlockBitLen)>1)
    {
        Radix_04R_DIT(this_FFT);
    }

    if ((this_FFT->NumberOfBitsNeeded-this_FFT->BlockBitLen)>0)
    {
        Radix_02R_DIT(this_FFT);
    }
}
#endif


#ifdef FACTOR_4_8
void FFT_DIT_01(FFT_Proc_Rec* this_FFT)     { Radix_02FS_DIT(this_FFT); }
void FFT_DIT_02(FFT_Proc_Rec* this_FFT)     { Radix_04FS_DIT(this_FFT); }
void FFT_DIT_03(FFT_Proc_Rec* this_FFT)     { Radix_08FS_DIT(this_FFT); }
void FFT_DIT_04(FFT_Proc_Rec* this_FFT)     { Radix_16FS_DIT(this_FFT); }
void FFT_DIT_05(FFT_Proc_Rec* this_FFT)     { Radix_04F_DIT(this_FFT); Radix_08F_DIT(this_FFT); }
void FFT_DIT_06(FFT_Proc_Rec* this_FFT)     { Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); }
void FFT_DIT_07(FFT_Proc_Rec* this_FFT)     { Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_08F_DIT(this_FFT); }
void FFT_DIT_08(FFT_Proc_Rec* this_FFT)     { Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); }
void FFT_DIT_09(FFT_Proc_Rec* this_FFT)     { Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_08F_DIT(this_FFT); }
void FFT_DIT_10(FFT_Proc_Rec* this_FFT)     { Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); }
void FFT_DIT_11(FFT_Proc_Rec* this_FFT)     { Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_08F_DIT(this_FFT); }
void FFT_DIT_12(FFT_Proc_Rec* this_FFT)     { Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); }
void FFT_DIT_13(FFT_Proc_Rec* this_FFT)     { Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_08F_DIT(this_FFT); }

void IFFT_DIT_01(FFT_Proc_Rec* this_FFT)    { Radix_02RS_DIT(this_FFT); }
void IFFT_DIT_02(FFT_Proc_Rec* this_FFT)    { Radix_04RS_DIT(this_FFT); }
void IFFT_DIT_03(FFT_Proc_Rec* this_FFT)    { Radix_08RS_DIT(this_FFT); }
void IFFT_DIT_04(FFT_Proc_Rec* this_FFT)    { Radix_16RS_DIT(this_FFT); }
void IFFT_DIT_05(FFT_Proc_Rec* this_FFT)    { Radix_04R_DIT(this_FFT); Radix_08R_DIT(this_FFT); }
void IFFT_DIT_06(FFT_Proc_Rec* this_FFT)    { Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); }
void IFFT_DIT_07(FFT_Proc_Rec* this_FFT)    { Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_08R_DIT(this_FFT); }
void IFFT_DIT_08(FFT_Proc_Rec* this_FFT)    { Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); }
void IFFT_DIT_09(FFT_Proc_Rec* this_FFT)    { Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_08R_DIT(this_FFT); }
void IFFT_DIT_10(FFT_Proc_Rec* this_FFT)    { Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); }
void IFFT_DIT_11(FFT_Proc_Rec* this_FFT)    { Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_08R_DIT(this_FFT); }
void IFFT_DIT_12(FFT_Proc_Rec* this_FFT)    { Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); }
void IFFT_DIT_13(FFT_Proc_Rec* this_FFT)    { Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_08R_DIT(this_FFT); }

void FFT_DIT_XX(FFT_Proc_Rec* this_FFT)
{
    while ((this_FFT->NumberOfBitsNeeded-this_FFT->BlockBitLen)>1)
    {
        Radix_04F_DIT(this_FFT);
    }
    if ((this_FFT->NumberOfBitsNeeded-this_FFT->BlockBitLen)>0)
    {
        Radix_02F_DIT(this_FFT);
    }
}

void IFFT_DIT_XX(FFT_Proc_Rec* this_FFT)
{
    while ((this_FFT->NumberOfBitsNeeded-this_FFT->BlockBitLen)>1)
    {
        Radix_04R_DIT(this_FFT);
    }
    if ((this_FFT->NumberOfBitsNeeded-this_FFT->BlockBitLen)>0)
    {
        Radix_02R_DIT(this_FFT);
    }
}

#endif


#ifdef FACTOR_8_4
void FFT_DIT_01(FFT_Proc_Rec* this_FFT)     { Radix_02FS_DIT(this_FFT); }
void FFT_DIT_02(FFT_Proc_Rec* this_FFT)     { Radix_04FS_DIT(this_FFT); }
void FFT_DIT_03(FFT_Proc_Rec* this_FFT)     { Radix_08FS_DIT(this_FFT); }
void FFT_DIT_04(FFT_Proc_Rec* this_FFT)     { Radix_16FS_DIT(this_FFT); }
void FFT_DIT_05(FFT_Proc_Rec* this_FFT)     { Radix_08F_DIT(this_FFT); Radix_04F_DIT(this_FFT); }
void FFT_DIT_06(FFT_Proc_Rec* this_FFT)     { Radix_08F_DIT(this_FFT); Radix_08F_DIT(this_FFT); }
void FFT_DIT_07(FFT_Proc_Rec* this_FFT)     { Radix_08F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); }
void FFT_DIT_08(FFT_Proc_Rec* this_FFT)     { Radix_08F_DIT(this_FFT); Radix_08F_DIT(this_FFT); Radix_04F_DIT(this_FFT);}
void FFT_DIT_09(FFT_Proc_Rec* this_FFT)     { Radix_08F_DIT(this_FFT); Radix_08F_DIT(this_FFT); Radix_08F_DIT(this_FFT);}
void FFT_DIT_10(FFT_Proc_Rec* this_FFT)     { Radix_08F_DIT(this_FFT); Radix_08F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); }
void FFT_DIT_11(FFT_Proc_Rec* this_FFT)     { Radix_08F_DIT(this_FFT); Radix_08F_DIT(this_FFT); Radix_08F_DIT(this_FFT); Radix_04F_DIT(this_FFT); }
void FFT_DIT_12(FFT_Proc_Rec* this_FFT)     { Radix_08F_DIT(this_FFT); Radix_08F_DIT(this_FFT); Radix_08F_DIT(this_FFT); Radix_08F_DIT(this_FFT); }
void FFT_DIT_13(FFT_Proc_Rec* this_FFT)     { Radix_08F_DIT(this_FFT); Radix_08F_DIT(this_FFT); Radix_08F_DIT(this_FFT); Radix_04F_DIT(this_FFT); Radix_04F_DIT(this_FFT); }

void IFFT_DIT_01(FFT_Proc_Rec* this_FFT)    { Radix_02RS_DIT(this_FFT); }
void IFFT_DIT_02(FFT_Proc_Rec* this_FFT)    { Radix_04RS_DIT(this_FFT); }
void IFFT_DIT_03(FFT_Proc_Rec* this_FFT)    { Radix_08RS_DIT(this_FFT); }
void IFFT_DIT_04(FFT_Proc_Rec* this_FFT)    { Radix_16RS_DIT(this_FFT); }
void IFFT_DIT_05(FFT_Proc_Rec* this_FFT)    { Radix_08R_DIT(this_FFT); Radix_04R_DIT(this_FFT); }
void IFFT_DIT_06(FFT_Proc_Rec* this_FFT)    { Radix_08R_DIT(this_FFT); Radix_08R_DIT(this_FFT); }
void IFFT_DIT_07(FFT_Proc_Rec* this_FFT)    { Radix_08R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); }
void IFFT_DIT_08(FFT_Proc_Rec* this_FFT)    { Radix_08R_DIT(this_FFT); Radix_08R_DIT(this_FFT); Radix_04R_DIT(this_FFT);}
void IFFT_DIT_09(FFT_Proc_Rec* this_FFT)    { Radix_08R_DIT(this_FFT); Radix_08R_DIT(this_FFT); Radix_08R_DIT(this_FFT);}
void IFFT_DIT_10(FFT_Proc_Rec* this_FFT)    { Radix_08R_DIT(this_FFT); Radix_08R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); }
void IFFT_DIT_11(FFT_Proc_Rec* this_FFT)    { Radix_08R_DIT(this_FFT); Radix_08R_DIT(this_FFT); Radix_08R_DIT(this_FFT); Radix_04R_DIT(this_FFT); }
void IFFT_DIT_12(FFT_Proc_Rec* this_FFT)    { Radix_08R_DIT(this_FFT); Radix_08R_DIT(this_FFT); Radix_08R_DIT(this_FFT); Radix_08R_DIT(this_FFT); }
void IFFT_DIT_13(FFT_Proc_Rec* this_FFT)    { Radix_08R_DIT(this_FFT); Radix_08R_DIT(this_FFT); Radix_08R_DIT(this_FFT); Radix_04R_DIT(this_FFT); Radix_04R_DIT(this_FFT); }

void FFT_DIT_XX(FFT_Proc_Rec* this_FFT)
{
    while ((this_FFT->NumberOfBitsNeeded-this_FFT->BlockBitLen)>2)
    {
        Radix_08F_DIT(this_FFT);
    }
    if ((this_FFT->NumberOfBitsNeeded-this_FFT->BlockBitLen)>1)
    {
        Radix_04F_DIT(this_FFT);
    }
    if ((this_FFT->NumberOfBitsNeeded-this_FFT->BlockBitLen)>0)
    {
        Radix_02F_DIT(this_FFT);
    }
}

void IFFT_DIT_XX(FFT_Proc_Rec* this_FFT)
{
    while ((this_FFT->NumberOfBitsNeeded-this_FFT->BlockBitLen)>2)
    {
        Radix_08R_DIT(this_FFT);
    }
    if ((this_FFT->NumberOfBitsNeeded-this_FFT->BlockBitLen)>1)
    {
        Radix_04R_DIT(this_FFT);
    }
    if ((this_FFT->NumberOfBitsNeeded-this_FFT->BlockBitLen)>0)
    {
        Radix_02R_DIT(this_FFT);
    }
}

#endif


#ifdef FACTOR_16
void FFT_DIT_01(FFT_Proc_Rec* this_FFT)     { Radix_02FS_DIT(this_FFT); }
void FFT_DIT_02(FFT_Proc_Rec* this_FFT)     { Radix_04FS_DIT(this_FFT); }
void FFT_DIT_03(FFT_Proc_Rec* this_FFT)     { Radix_08FS_DIT(this_FFT); }
void FFT_DIT_04(FFT_Proc_Rec* this_FFT)     { Radix_16FS_DIT(this_FFT); }
void FFT_DIT_05(FFT_Proc_Rec* this_FFT)     { Radix_08F_DIT(this_FFT); Radix_04F_DIT(this_FFT); }
void FFT_DIT_06(FFT_Proc_Rec* this_FFT)     { Radix_08F_DIT(this_FFT); Radix_08F_DIT(this_FFT); }
void FFT_DIT_07(FFT_Proc_Rec* this_FFT)     { Radix_16F_DIT(this_FFT); Radix_08F_DIT(this_FFT); }
void FFT_DIT_08(FFT_Proc_Rec* this_FFT)     { Radix_16F_DIT(this_FFT); Radix_16F_DIT(this_FFT); }
void FFT_DIT_09(FFT_Proc_Rec* this_FFT)     { Radix_16F_DIT(this_FFT); Radix_08F_DIT(this_FFT); Radix_04F_DIT(this_FFT); }
void FFT_DIT_10(FFT_Proc_Rec* this_FFT)     { Radix_16F_DIT(this_FFT); Radix_08F_DIT(this_FFT); Radix_08F_DIT(this_FFT); }
void FFT_DIT_11(FFT_Proc_Rec* this_FFT)     { Radix_16F_DIT(this_FFT); Radix_16F_DIT(this_FFT); Radix_08F_DIT(this_FFT); }
void FFT_DIT_12(FFT_Proc_Rec* this_FFT)     { Radix_16F_DIT(this_FFT); Radix_16F_DIT(this_FFT); Radix_16F_DIT(this_FFT); }
void FFT_DIT_13(FFT_Proc_Rec* this_FFT)     { Radix_16F_DIT(this_FFT); Radix_16F_DIT(this_FFT); Radix_08F_DIT(this_FFT); Radix_04F_DIT(this_FFT); }

void IFFT_DIT_01(FFT_Proc_Rec* this_FFT)    { Radix_02R_DIT(this_FFT); }
void IFFT_DIT_02(FFT_Proc_Rec* this_FFT)    { Radix_04R_DIT(this_FFT); }
void IFFT_DIT_03(FFT_Proc_Rec* this_FFT)    { Radix_08R_DIT(this_FFT); }
void IFFT_DIT_04(FFT_Proc_Rec* this_FFT)    { Radix_16R_DIT(this_FFT); }
void IFFT_DIT_05(FFT_Proc_Rec* this_FFT)    { Radix_08R_DIT(this_FFT); Radix_04R_DIT(this_FFT); }
void IFFT_DIT_06(FFT_Proc_Rec* this_FFT)    { Radix_08R_DIT(this_FFT); Radix_08R_DIT(this_FFT); }
void IFFT_DIT_07(FFT_Proc_Rec* this_FFT)    { Radix_16R_DIT(this_FFT); Radix_08R_DIT(this_FFT); }
void IFFT_DIT_08(FFT_Proc_Rec* this_FFT)    { Radix_16R_DIT(this_FFT); Radix_16R_DIT(this_FFT); }
void IFFT_DIT_09(FFT_Proc_Rec* this_FFT)    { Radix_16R_DIT(this_FFT); Radix_08R_DIT(this_FFT); Radix_04R_DIT(this_FFT); }
void IFFT_DIT_10(FFT_Proc_Rec* this_FFT)    { Radix_16R_DIT(this_FFT); Radix_08R_DIT(this_FFT); Radix_08R_DIT(this_FFT); }
void IFFT_DIT_11(FFT_Proc_Rec* this_FFT)    { Radix_16R_DIT(this_FFT); Radix_16R_DIT(this_FFT); Radix_08R_DIT(this_FFT); }
void IFFT_DIT_12(FFT_Proc_Rec* this_FFT)    { Radix_16R_DIT(this_FFT); Radix_16R_DIT(this_FFT); Radix_16R_DIT(this_FFT); }
void IFFT_DIT_13(FFT_Proc_Rec* this_FFT)    { Radix_16R_DIT(this_FFT); Radix_16R_DIT(this_FFT); Radix_08R_DIT(this_FFT); Radix_04R_DIT(this_FFT); }

void FFT_DIT_XX(FFT_Proc_Rec* this_FFT)
{
    while ((this_FFT->NumberOfBitsNeeded-this_FFT->BlockBitLen)>3)
    {
        Radix_16F_DIT(this_FFT);
    }
    if ((this_FFT->NumberOfBitsNeeded-this_FFT->BlockBitLen)>2)
    {
        Radix_08F_DIT(this_FFT);
    }
    if ((this_FFT->NumberOfBitsNeeded-this_FFT->BlockBitLen)>1)
    {
        Radix_04F_DIT(this_FFT);
    }
    if ((this_FFT->NumberOfBitsNeeded-this_FFT->BlockBitLen)>0)
    {
        Radix_02F_DIT(this_FFT);
    }
}

void IFFT_DIT_XX(FFT_Proc_Rec* this_FFT)
{
    while ((this_FFT->NumberOfBitsNeeded-this_FFT->BlockBitLen)>2)
    {
        Radix_16R_DIT(this_FFT);
    }
    if ((this_FFT->NumberOfBitsNeeded-this_FFT->BlockBitLen)>2)
    {
        Radix_08R_DIT(this_FFT);
    }
    if ((this_FFT->NumberOfBitsNeeded-this_FFT->BlockBitLen)>1)
    {
        Radix_04R_DIT(this_FFT);
    }
    if ((this_FFT->NumberOfBitsNeeded-this_FFT->BlockBitLen)>0)
    {
        Radix_02R_DIT(this_FFT);
    }
}

#endif


//void (*  FFT_DIT[32+1])(FFT_Proc_Rec*) = {nullptr,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX};
//void (* IFFT_DIT[32+1])(FFT_Proc_Rec*) = {nullptr, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX};

void (*  FFT_DIT[32+1])(FFT_Proc_Rec*) = {nullptr,  FFT_DIT_01,  FFT_DIT_02,  FFT_DIT_03,  FFT_DIT_04,  FFT_DIT_05,  FFT_DIT_06,  FFT_DIT_07,  FFT_DIT_08,  FFT_DIT_09,  FFT_DIT_10,  FFT_DIT_11,  FFT_DIT_12,  FFT_DIT_13,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX,  FFT_DIT_XX};
void (* IFFT_DIT[32+1])(FFT_Proc_Rec*) = {nullptr, IFFT_DIT_01, IFFT_DIT_02, IFFT_DIT_03, IFFT_DIT_04, IFFT_DIT_05, IFFT_DIT_06, IFFT_DIT_07, IFFT_DIT_08, IFFT_DIT_09, IFFT_DIT_10, IFFT_DIT_11, IFFT_DIT_12, IFFT_DIT_13, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX, IFFT_DIT_XX};


void FFT_DIT_Complex(FFT_Proc_Rec* this_FFT_plan)
{
    if (((this_FFT_plan->NumberOfBitsNeeded <= nFFT_MAX_FFT_BIT_LENGTH) && (this_FFT_plan->NumberOfBitsNeeded > 0)) && (this_FFT_plan->DComplex != nullptr))
    {
        this_FFT_plan->FFT = &FFT_PreCalc_Data_Rec[this_FFT_plan->NumberOfBitsNeeded];
        this_FFT_plan->BlockBitLen = 0;
    }
    else
        throw(-1);

    if (this_FFT_plan->NumberOfBitsNeeded>4)
    {
        int32_t sp_i, sp_j;
        for (sp_i = 1; sp_i < this_FFT_plan->FFT->length; sp_i++)
        {
            sp_j = RevBits[this_FFT_plan->NumberOfBitsNeeded][sp_i];

            if (sp_i < sp_j)
                swap(this_FFT_plan->DComplex[sp_i],this_FFT_plan->DComplex[sp_j]);
        }
    }

    FFT_DIT[this_FFT_plan->NumberOfBitsNeeded](this_FFT_plan);
}


void IFFT_DIT_Complex(FFT_Proc_Rec* this_FFT_plan)
{
    if (((this_FFT_plan->NumberOfBitsNeeded<=nFFT_MAX_FFT_BIT_LENGTH) && (this_FFT_plan->NumberOfBitsNeeded>0)) && (this_FFT_plan->DComplex != nullptr))
    {
        this_FFT_plan->FFT = &FFT_PreCalc_Data_Rec[this_FFT_plan->NumberOfBitsNeeded];
        this_FFT_plan->BlockBitLen = 0;
    }
    else
        throw(-1);

    if (this_FFT_plan->NumberOfBitsNeeded>4)
    {
        int32_t sp_i, sp_j;
        for (sp_i = 1; sp_i < this_FFT_plan->FFT->length; sp_i++)
        {
            sp_j = RevBits[this_FFT_plan->NumberOfBitsNeeded][sp_i];

            if (sp_i < sp_j)
                swap(this_FFT_plan->DComplex[sp_i],this_FFT_plan->DComplex[sp_j]);
        }
    }

    IFFT_DIT[this_FFT_plan->NumberOfBitsNeeded](this_FFT_plan);
}


void FFT_DIT_Real(FFT_Proc_Rec* this_FFT_plan)
{

    if (((this_FFT_plan->NumberOfBitsNeeded<=nFFT_MAX_FFT_BIT_LENGTH) && (this_FFT_plan->NumberOfBitsNeeded>1)) && (this_FFT_plan->DComplex != nullptr))
    {
        this_FFT_plan->NumberOfBitsNeeded--;
        this_FFT_plan->FFT = &FFT_PreCalc_Data_Rec[this_FFT_plan->NumberOfBitsNeeded];
        this_FFT_plan->BlockBitLen = 0;
    }
    else
        throw(-1);

    tDComplex* this_A_Arr = A_Arr[this_FFT_plan->NumberOfBitsNeeded+1];

    if (this_FFT_plan->NumberOfBitsNeeded>4)
    {
        int32_t sp_i, sp_j;
        for (sp_i = 1; sp_i < this_FFT_plan->FFT->length; sp_i++)
        {
            sp_j = RevBits[this_FFT_plan->NumberOfBitsNeeded][sp_i];

            if (sp_i < sp_j)
                swap(this_FFT_plan->DComplex[sp_i],this_FFT_plan->DComplex[sp_j]);
        }
    }

    FFT_DIT[this_FFT_plan->NumberOfBitsNeeded](this_FFT_plan);

    //========================================================================
    // Resolve HalfN Real FFT Output to N Complex FFT Output.
    //========================================================================

    int32_t D0 = 1;
    int32_t D1 = this_FFT_plan->FFT->length - 1;
    int32_t D2 = this_FFT_plan->FFT->length + 1;
    int32_t D3 = this_FFT_plan->FFT->length + D1;

    for (int32_t J = 1; J <= this_FFT_plan->FFT->length_half; J++)
    {
        tDComplex X0 = (this_FFT_plan->DComplex[D0] + this_FFT_plan->DComplex[D1]) * 0.5;
        tDComplex X1 = (this_FFT_plan->DComplex[D0] - this_FFT_plan->DComplex[D1]) * 0.5;

        tDComplex V = DComplex(X0.Im, -X1.Re) * this_A_Arr[D0];

        tDComplex W = DComplex(X0.Re, X1.Im);

        this_FFT_plan->DComplex[D0] = (W + V);

        this_FFT_plan->DComplex[D2] = (W - V);

        this_FFT_plan->DComplex[D3--] = this_FFT_plan->DComplex[D0++].conj();

        this_FFT_plan->DComplex[D1--] = this_FFT_plan->DComplex[D2++].conj();
    }

    this_FFT_plan->DComplex[this_FFT_plan->FFT->length] = DComplex(this_FFT_plan->DComplex[0].Re - this_FFT_plan->DComplex[0].Im, 0.);

    this_FFT_plan->DComplex[0] = DComplex(this_FFT_plan->DComplex[0].Re + this_FFT_plan->DComplex[0].Im, 0.);

    this_FFT_plan->NumberOfBitsNeeded++;
    this_FFT_plan->FFT = &FFT_PreCalc_Data_Rec[this_FFT_plan->NumberOfBitsNeeded];
}


void nFFT_Init(int32_t desired_max_fft_bit_length)
{
    nFFT_MAX_FFT_BIT_LENGTH = desired_max_fft_bit_length;

    if (nFFT_MAX_FFT_BIT_LENGTH > 32)
        throw(-1);

    for (int32_t nf_i = 0; nf_i < 32; nf_i++)
    {
        cos_pi_over_32[nf_i] = std::cos(nf_i * OneOver[32] * Pi);
    }

    uint32_t nFFT_MAX_FFT_LENGTH = 1 << nFFT_MAX_FFT_BIT_LENGTH;

    double delta = (-TwoPi) / nFFT_MAX_FFT_LENGTH;

    //=========================================================================================================================
    // Create twiddle factor lookup arrays.
    //=========================================================================================================================
    for (int32_t nf_i = 0; nf_i <= nFFT_MAX_FFT_BIT_LENGTH; nf_i++)
    {
        A_Arr[nf_i] = new tDComplex[sizeof(tDComplex)<< nf_i];
        A_Arr_Conj[nf_i] = new tDComplex[sizeof(tDComplex)<< nf_i];
        RevBits[nf_i] = new int32_t[sizeof(int32_t) << nf_i];
    }
    //=========================================================================================================================
    // Calculate twiddle factor values for MAX_FFT_LENGTH array.
    //=========================================================================================================================
    for (uint32_t nf_j = 0; nf_j < nFFT_MAX_FFT_LENGTH; nf_j++)
    {
        A_Arr[nFFT_MAX_FFT_BIT_LENGTH][nf_j] = complex_exp(nf_j*delta);
        A_Arr_Conj[nFFT_MAX_FFT_BIT_LENGTH][nf_j] = A_Arr[nFFT_MAX_FFT_BIT_LENGTH][nf_j].conj();
    }
    //=========================================================================================================================
    // Backfill shorter twiddle factor arrays quickly using shifted indices.
    //=========================================================================================================================
    for (int32_t nf_i = 0; nf_i < nFFT_MAX_FFT_BIT_LENGTH; nf_i++)
    {
        for (int32_t nf_j = 0; nf_j < (1 << nf_i); nf_j++)
        {
            A_Arr[nf_i][nf_j] = A_Arr[nFFT_MAX_FFT_BIT_LENGTH][nf_j << (nFFT_MAX_FFT_BIT_LENGTH - nf_i)];
            A_Arr_Conj[nf_i][nf_j] = A_Arr[nf_i][nf_j].conj();
        }

        DATA32 sp_i, sp_j;
        for (sp_i.Integer = 0; sp_i.Integer < (1 << nf_i); sp_i.Integer++)
        {
            sp_j.Bytes[3] = BitReversedLookupTable[sp_i.Bytes[0]];
            sp_j.Bytes[2] = BitReversedLookupTable[sp_i.Bytes[1]];
            sp_j.Bytes[1] = BitReversedLookupTable[sp_i.Bytes[2]];
            sp_j.Bytes[0] = BitReversedLookupTable[sp_i.Bytes[3]];
            sp_j.Cardinal >>= (32-nf_i);
            RevBits[nf_i][sp_i.Integer] = sp_j.Integer;
        }

    }
    //=========================================================================================================================
}

void nFFT_Cleanup()
{
    for (int32_t nf_i = 0; nf_i <= nFFT_MAX_FFT_BIT_LENGTH; nf_i++)
    {
        if (A_Arr[nf_i] != nullptr)
        {
            delete[] A_Arr[nf_i];
        }

        if (A_Arr_Conj[nf_i] != nullptr)
        {
            delete[] A_Arr_Conj[nf_i];
        }

        if (RevBits[nf_i] != nullptr)
        {
            delete[] RevBits[nf_i];
        }
    }
}
