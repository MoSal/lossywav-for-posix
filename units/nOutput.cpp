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

#include <iostream>
#include <cmath>

#if defined(_MSC_VER)
#include <intrin.h>
#pragma intrinsic(_BitScanForward)
#pragma intrinsic(_BitScanReverse)
#endif

#include "nOutput.h"
#include "nCore.h"
#include "nMaths.h"
#include "nSpreading.h"
#include "fftw_interface.h"
#include "nSGNS.h"
#include "nParameter.h"

// Global
uint8_t* bit_removal_history;

const char hyphen_string[256] = "---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------";
const char bits_filled[256]   = "OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO";
const char bits_empty[256]    = "...............................................................................................................................................................................................................................................................";
const char space_string[256]  = "                                                                                                                                                                                                                                                               ";

int64_t DATA_block_lsb [MAX_CHANNELS][34];
int64_t DATA_block_msb [MAX_CHANNELS][34];
int64_t DATA_sample_lsb [MAX_CHANNELS][34];
int64_t DATA_sample_msb [MAX_CHANNELS][34];

int64_t BTRD_block_lsb [MAX_CHANNELS][34];
int64_t BTRD_block_msb [MAX_CHANNELS][34];
int64_t BTRD_sample_lsb [MAX_CHANNELS][34];
int64_t BTRD_sample_msb [MAX_CHANNELS][34];

int64_t CORR_block_lsb [MAX_CHANNELS][34];
int64_t CORR_block_msb [MAX_CHANNELS][34];
int64_t CORR_sample_lsb [MAX_CHANNELS][34];
int64_t CORR_sample_msb [MAX_CHANNELS][34];

std::string head_bar;
std::string top_bar;
std::string mid_bar;
int Display_Width;
int bar_length;
std::string Titles [16];
std::string Header;


void Make_Bars(int32_t MB_Width, int32_t MB_Item, int32_t MB_Value, int32_t MB_Count, int32_t MB_Summary = 0)
{
    int32_t mb_i;

    if (MB_Summary > 0)
    {
        bar_length = ((MB_Width - (MB_Item + 2) - (MB_Value + 1) * MB_Count - (MB_Summary + 1)) / MB_Count) - 1;
    }
    else
    {
        bar_length = ((MB_Width - (MB_Item + 2) - (MB_Value + 1) * MB_Count) / MB_Count) - 1;
    }

    top_bar = "+" + std::string(hyphen_string, MB_Item) + '+';
    mid_bar = top_bar;
    head_bar = top_bar + std::string(hyphen_string, Display_Width - top_bar.length() - 1) + "+";

    if (bar_length > 0)
        for (mb_i = 1; mb_i <= MB_Count;  ++mb_i)
        {
            top_bar += std::string(hyphen_string, MB_Value + 1 + bar_length) + "+";
            mid_bar += std::string(hyphen_string, MB_Value) + "+" + std::string(hyphen_string, bar_length) + "+";
        }
    else
        for (mb_i = 1; mb_i <= MB_Count;  ++mb_i)
        {
            top_bar += std::string(hyphen_string, MB_Value) + "+";
            mid_bar += std::string(hyphen_string, MB_Value) + "+";
        }

    if (MB_Summary > 0)
    {
        top_bar += std::string(hyphen_string, MB_Summary) + "+";
        mid_bar += std::string(hyphen_string, MB_Summary) + "+";
    }

    if (bar_length > 0)
    {
        for (mb_i = 0; mb_i < MB_Count;  ++mb_i)
        {
            Titles[mb_i] += std::string(space_string, bar_length + 2 + MB_Value - Titles[mb_i].length());
        }
        Header += std::string(space_string, MB_Count * (bar_length + 2 + MB_Value) - Header.length());
    }
}


void remove_bits_detailed_output(std::ostream& ToOutput)
{
    uint64_t rb_i, rb_l;
    int32_t rb_j, rb_k;
    uint64_t ss_k;
    double btr_row_tot;
    int32_t channel;
    std::string time_string;
    ToOutput << std::endl << "Detailed bits-to-remove data per channel per codec-block.\n";
    ss_k = (Display_Width - (8 + 2) - (4 + 1)) / (3 * Global.Channels);
    Titles[0] = "";
    Titles[1] = "";
    Titles[2] = "";
    Header = "";
    Make_Bars(Display_Width, 8, (3 * Global.Channels - 1), ss_k, 4);
    ToOutput << top_bar << std::endl << "|  Time  |";

    for (rb_i = 0; rb_i < ss_k; rb_i ++)
        for (channel = 0; channel < Global.Channels; channel ++)
        {
            ToOutput << '#' << channel;

            if (channel == Global.Channels - 1)
            {
                ToOutput << '|';
            }
            else
            {
                ToOutput << ',';
            }
        }

    ToOutput << "Tot |" << std::endl << mid_bar << std::endl;
    ss_k = ss_k * Global.Channels;
    btr_row_tot = 0;
    rb_j = 0;

    for (rb_i = 0; rb_i < Global.blocks_processed; rb_i ++)
        for (channel = 0; channel < Global.Channels; channel ++)
        {
            rb_l = rb_j % ss_k;

            if (rb_l == 0)
            {
                time_string_make(time_string, double(rb_i) * Global.Codec_Block.Size / Global.sample_rate);
                btr_row_tot = 0;
                ToOutput << '|' << time_string << '|';
            }

            rb_k = bit_removal_history[rb_j];
            btr_row_tot = btr_row_tot + rb_k;
            ToOutput << std::setw(2) << rb_k;

            if (channel == Global.Channels - 1)
            {
                ToOutput << '|';
            }
            else
            {
                ToOutput << ',';
            }

            if (rb_l == ss_k - 1)
            {
                ToOutput << std::fixed << std::setw(4) << std::setprecision(0) << btr_row_tot << '|' << std::endl;
            }

            ++ rb_j;
        }

    rb_l = rb_j % ss_k;

    if ((rb_l < ss_k - 1) && (rb_l > 0))
    {
        for (rb_i = 0; rb_i <= (ss_k - rb_l - 1) / Global.Channels; rb_i ++)
            for (channel = 0; channel < Global.Channels; channel ++)
                if (channel < Global.Channels - 1)
                {
                    ToOutput << "   ";
                }
                else
                {
                    ToOutput << "  |";
                }

        ToOutput << std::fixed << std::setw(4) << std::setprecision(0) << btr_row_tot << '|' << std::endl;
    }

    ToOutput << mid_bar << std::endl;
}


#if defined(__GNUC__)
    static inline uint32_t clz(uint32_t x)
    {
        return __builtin_clz(x); // requires -march=486 on g++ command line on PC
    }
    static inline uint32_t ctz(uint32_t x)
    {
        return __builtin_ctz(x); // requires -march=486 on g++ command line on PC
    }
#elif defined(_MSC_VER)
    static inline uint32_t clz(uint32_t x)
    {
        uint32_t r = 0;
        _BitScanForward(&r, x);
        return r;
    }
    static inline uint32_t ctz(uint32_t x)
    {
        uint32_t long r = 0;
        _BitScanReverse(&r, x);
        return r;
    }
#else
    static inline uint32_t popcnt(uint32_t x)
    {
        x -= ((x >> 1) & 0x55555555);
        x = (((x >> 2) & 0x33333333) + (x & 0x33333333));
        x = (((x >> 4) + x) & 0x0f0f0f0f);
        x += (x >> 8);
        x += (x >> 16);
        return x & 0x0000003f;
    }
    static inline uint32_t clz(uint32_t x)
    {
        x |= (x >> 1);
        x |= (x >> 2);
        x |= (x >> 4);
        x |= (x >> 8);
        x |= (x >> 16);
        return 32 - popcnt(x);
    }
    static inline uint32_t ctz(uint32_t x)
    {
        return popcnt((x & -int32_t(x)) - 1);
    }
#endif


inline int32_t Bit_Scan_Left(int32_t BSF_Int)
{
    // __asm {    bsf BSF_Int, Eax }
    return int(ctz(uint32_t(BSF_Int)));
}

inline int32_t Bit_Scan_Right(int32_t BSR_Int)
{
    // __asm { bsr BSR_Int, Eax    }
    return int(clz(uint32_t(BSR_Int)));
}


void lsb_analysis()
{
    for (int32_t sa_j = 0; sa_j < Global.Channels; sa_j ++)
    {
        int32_t DATA_block_value = 0;

        for (int32_t sa_i = 0; sa_i < AudioData.Size.This; sa_i ++)
        {
            int32_t temp_val = fabs(AudioData.WAVEPTR[THIS_CODEC_BLOCK][sa_j][sa_i].Integers[0]);

            if (temp_val == 0)
            {
                ++ DATA_sample_lsb[sa_j][0];
                ++ DATA_sample_msb[sa_j][0];
            }
            else
            {
                ++ DATA_sample_lsb[sa_j][Bit_Scan_Left(temp_val) + 1];
                ++ DATA_sample_msb[sa_j][32 - Bit_Scan_Right(temp_val) + 1];
            }

            DATA_block_value = DATA_block_value | temp_val;
        }

        if (DATA_block_value == 0)
        {
            ++ DATA_block_lsb[sa_j][0];
            ++ DATA_block_msb[sa_j][0];
        }
        else
        {
            ++ DATA_block_lsb[sa_j][Bit_Scan_Left(DATA_block_value) + 1];
            ++ DATA_block_msb[sa_j][32 - Bit_Scan_Right(DATA_block_value) + 1];
        }

        int32_t BTRD_block_value = 0;

        for (int32_t sa_i = 0; sa_i < AudioData.Size.This; sa_i ++)
        {
            int32_t temp_val = fabs(AudioData.BTRDPTR[THIS_CODEC_BLOCK][sa_j][sa_i].Integers[0]);

            if (temp_val == 0)
            {
                ++ BTRD_sample_lsb[sa_j][0];
                ++ BTRD_sample_msb[sa_j][0];
            }
            else
            {
                ++ BTRD_sample_lsb[sa_j][Bit_Scan_Left(temp_val) + 1];
                ++ BTRD_sample_msb[sa_j][32 - Bit_Scan_Right(temp_val) + 1];
            }

            BTRD_block_value = BTRD_block_value | temp_val;
        }

        if (BTRD_block_value == 0)
        {
            ++ BTRD_block_lsb[sa_j][0];
            ++ BTRD_block_msb[sa_j][0];
        }
        else
        {
            ++ BTRD_block_lsb[sa_j][Bit_Scan_Left(BTRD_block_value) + 1];
            ++ BTRD_block_msb[sa_j][32 - Bit_Scan_Right(BTRD_block_value) + 1];
        }

        int32_t CORR_block_value = 0;

        for (int32_t sa_i = 0; sa_i < AudioData.Size.This; sa_i ++)
        {
            int32_t temp_val = fabs(AudioData.CORRPTR[THIS_CODEC_BLOCK][sa_j][sa_i].Integers[0]);

            if (temp_val == 0)
            {
                ++ CORR_sample_lsb[sa_j][0];
                ++ CORR_sample_msb[sa_j][0];
            }
            else
            {
                ++ CORR_sample_lsb[sa_j][Bit_Scan_Left(temp_val) + 1];
                ++ CORR_sample_msb[sa_j][32 - Bit_Scan_Right(temp_val) + 1];
            }

            CORR_block_value = CORR_block_value | temp_val;
        }

        if (CORR_block_value == 0)
        {
            ++ CORR_block_lsb[sa_j][0];
            ++ CORR_block_msb[sa_j][0];
        }
        else
        {
            ++ CORR_block_lsb[sa_j][Bit_Scan_Left(CORR_block_value) + 1];
            ++ CORR_block_msb[sa_j][32 - Bit_Scan_Right(CORR_block_value) + 1];
        }
    }
}


void Write_bitdist(std::ostream& ToOutput)
{
    int32_t nd_i, nd_j;
    float this_bit_removed_percentage;
    float this_bit_lost_percentage;
    int32_t bits_removed_filled;
    int32_t bits_lost_filled;

    if (!parameters.output.perchannel)
    {
        Titles[0] = "";
        Titles[1] = "";
        Titles[2] = "";
        Header = "| Distribution of bits removed from each codec-block";
        Make_Bars(Display_Width, 3, 6, 1);
        ToOutput << std::endl << "Bits removed distribution." << std::endl << head_bar << std::endl << "|Bit" << Header << '|' << std::endl << mid_bar << std::endl;

        for (nd_i = 0; nd_i < Global.bits_per_sample; nd_i ++)
        {
            this_bit_removed_percentage = 0;

            for (nd_j = 0; nd_j < Global.Channels; nd_j ++)
            {
                this_bit_removed_percentage += float(Stats.bits_removed[nd_j][nd_i] * Global.blocks_processed_recip * OneOver[Global.Channels]);
            }

            bits_removed_filled = nRoundEvenInt32(this_bit_removed_percentage * bar_length);
            ToOutput << '|' << std::setw(3) << nd_i;
            ToOutput << std::setw(5) << std::fixed << std::setprecision(1) << this_bit_removed_percentage * 100 << "%|";
            ToOutput << std::string(bits_filled, bits_removed_filled)
                     << std::string(bits_empty, bar_length - bits_removed_filled)
                     << '|' << std::endl;
        }

        ToOutput << mid_bar << std::endl;

        Titles[0] = "";
        Titles[1] = "";
        Titles[2] = "";
        Header = "| Distribution of bits lost from each codec-block";
        Make_Bars(Display_Width, 3, 6, 1);
        ToOutput << std::endl << "Bits lost distribution." << std::endl << head_bar << std::endl << "|Bit" << Header << '|' << std::endl << mid_bar << std::endl;

        for (nd_i = 0; nd_i < Global.bits_per_sample; nd_i ++)
        {
            this_bit_lost_percentage = 0;

            for (nd_j = 0; nd_j < Global.Channels; nd_j ++)
            {
                this_bit_lost_percentage += float(Stats.bits_lost[nd_j][nd_i] * Global.blocks_processed_recip * OneOver[Global.Channels]);
            }

            bits_lost_filled = nRoundEvenInt32(this_bit_lost_percentage * bar_length);
            ToOutput << '|' << std::setw(3) << nd_i;
            ToOutput << std::setw(5) << std::fixed << std::setprecision(1) << this_bit_lost_percentage * 100 << "%|";
            ToOutput << std::string(bits_filled, bits_lost_filled)
                     << std::string(bits_empty, bar_length - bits_lost_filled)
                     << '|' << std::endl;
        }

        ToOutput << mid_bar << std::endl;
    }
    else
    {
        for (nd_i = 0; nd_i < Global.Channels; nd_i ++)
        {
            Titles[nd_i] = std::string("| Ch. #") + char('0' + nd_i);
        }

        Make_Bars(Display_Width, 3, 6, Global.Channels);
        ToOutput << std::endl << "Distribution of bits removed from each codec-block / channel." << std::endl << top_bar << std::endl << "|Bit";

        for (nd_i = 0; nd_i < Global.Channels; nd_i ++)
        {
            ToOutput << Titles[nd_i];
        }

        ToOutput << '|' << std::endl << mid_bar << std::endl;

        for (nd_i = 0; nd_i < Global.bits_per_sample; nd_i ++)
        {
            ToOutput << '|' << std::setw(3) << nd_i << '|';

            for (nd_j = 0; nd_j < Global.Channels; nd_j ++)
            {
                this_bit_removed_percentage = float(Stats.bits_removed[nd_j][nd_i] * Global.blocks_processed_recip);
                bits_removed_filled = nRoundEvenInt32(this_bit_removed_percentage * bar_length);
                ToOutput << std::setw(5) << std::fixed << std::setprecision(1) << this_bit_removed_percentage * 100 << "%|";
                ToOutput << std::string(bits_filled, bits_removed_filled)
                         << std::string(bits_empty, bar_length - bits_removed_filled) << '|';
            }

            ToOutput << std::endl;
        }

        ToOutput << mid_bar << std::endl;

        for (nd_i = 0; nd_i < Global.Channels; nd_i ++)
        {
            Titles[nd_i] = std::string("| Ch. #") + char('0' + nd_i);
        }

        Make_Bars(Display_Width, 3, 6, Global.Channels);
        ToOutput << std::endl << "Distribution of bits lost from each codec-block / channel." << std::endl << top_bar << std::endl << "|Bit";

        for (nd_i = 0; nd_i < Global.Channels; nd_i ++)
        {
            ToOutput << Titles[nd_i];
        }

        ToOutput << '|' << std::endl << mid_bar << std::endl;

        for (nd_i = 0; nd_i < Global.bits_per_sample; nd_i ++)
        {
            ToOutput << '|' << std::setw(3) << nd_i << '|';

            for (nd_j = 0; nd_j < Global.Channels; nd_j ++)
            {
                this_bit_lost_percentage = float(Stats.bits_lost[nd_j][nd_i] * Global.blocks_processed_recip);
                bits_lost_filled = nRoundEvenInt32(this_bit_lost_percentage * bar_length);
                ToOutput << std::setw(5) << std::fixed << std::setprecision(1) << this_bit_lost_percentage * 100 << "%|";
                ToOutput << std::string(bits_filled, bits_lost_filled)
                         << std::string(bits_empty, bar_length - bits_lost_filled) << '|';
            }

            ToOutput << std::endl;
        }

        ToOutput << mid_bar << std::endl;
    }
}


void Write_blockdist(std::ostream& ToOutput)
{
    int32_t nd_i, nd_j;
    int32_t bits_filled_wave;
    int32_t bits_filled_btrd;
    int32_t bits_filled_corr;
    double this_DATA_percentage;
    double this_BTRD_percentage;
    double this_CORR_percentage;
    Titles[0] = "|Input lsb distribution";
    Titles[1] = "|Lossy lsb distribution";
    Titles[2] = "|LWCDF lsb distribution";
    Header = "";
    Make_Bars(Display_Width, 3, 6, 3);
    ToOutput << std::endl << "Codec-block least significant bit (lsb) distribution." << std::endl << top_bar
             << std::endl << "|Bit" << Titles[0] << Titles[1] << Titles[2] << '|' << std::endl << mid_bar << std::endl;

    for (nd_i = 0; nd_i <= Global.bits_per_sample; nd_i ++)
    {
        this_DATA_percentage = 0;
        this_BTRD_percentage = 0;
        this_CORR_percentage = 0;

        for (nd_j = 0; nd_j < Global.Channels; nd_j ++)
        {
            this_DATA_percentage += DATA_block_lsb[nd_j][nd_i] * Global.blocks_processed_recip * OneOver[Global.Channels];
            this_BTRD_percentage += BTRD_block_lsb[nd_j][nd_i] * Global.blocks_processed_recip * OneOver[Global.Channels];
            this_CORR_percentage += CORR_block_lsb[nd_j][nd_i] * Global.blocks_processed_recip * OneOver[Global.Channels];
        }

        bits_filled_wave = nRoundEvenInt32(std::min(1.0,this_DATA_percentage) * bar_length);
        bits_filled_btrd = nRoundEvenInt32(std::min(1.0,this_BTRD_percentage) * bar_length);
        bits_filled_corr = nRoundEvenInt32(std::min(1.0,this_CORR_percentage) * bar_length);

        if (nd_i == 0)
        {
            ToOutput << "|NUL";
        }
        else
        {
            ToOutput << '|' << std::setw(3) << (nd_i - 1);
        }

        ToOutput << '|' << std::setw(5) << NumToStr(this_DATA_percentage * 100, 1) << "%|";
        ToOutput << std::string(bits_filled, bits_filled_wave) << std::string(bits_empty, bar_length - bits_filled_wave) << '|';
        ToOutput << std::setw(5) << NumToStr(this_BTRD_percentage * 100, 1) << "%|";
        ToOutput << std::string(bits_filled, bits_filled_btrd) << std::string(bits_empty, bar_length - bits_filled_btrd) << '|';
        ToOutput << std::setw(5) << NumToStr(this_CORR_percentage * 100, 1) << "%|";
        ToOutput << std::string(bits_filled, bits_filled_corr) << std::string(bits_empty, bar_length - bits_filled_corr) << '|' << std::endl;
    }

    Titles[0] = "|Input msb distribution";
    Titles[1] = "|Lossy msb distribution";
    Titles[2] = "|LWCDF msb distribution";
    Header = "";

    Make_Bars(Display_Width, 3, 6, 3);
    ToOutput << mid_bar << std::endl << std::endl << "Codec-block most significant bit (msb) distribution." << std::endl
             << top_bar << std::endl << "|Bit" << Titles[0] << Titles[1] << Titles[2] << '|' << std::endl << mid_bar << std::endl;

    for (nd_i = 0; nd_i <= Global.bits_per_sample; nd_i ++)
    {
        this_DATA_percentage = 0;
        this_BTRD_percentage = 0;
        this_CORR_percentage = 0;

        for (nd_j = 0; nd_j < Global.Channels; nd_j ++)
        {
            this_DATA_percentage += DATA_block_msb[nd_j][nd_i] * Global.blocks_processed_recip * OneOver[Global.Channels];
            this_BTRD_percentage += BTRD_block_msb[nd_j][nd_i] * Global.blocks_processed_recip * OneOver[Global.Channels];
            this_CORR_percentage += CORR_block_msb[nd_j][nd_i] * Global.blocks_processed_recip * OneOver[Global.Channels];
        }

        bits_filled_wave = nRoundEvenInt32(std::min(1.0,this_DATA_percentage) * bar_length);
        bits_filled_btrd = nRoundEvenInt32(std::min(1.0,this_BTRD_percentage) * bar_length);
        bits_filled_corr = nRoundEvenInt32(std::min(1.0,this_CORR_percentage) * bar_length);

        if (nd_i == 0)
        {
            ToOutput << "|NUL";
        }
        else
        {
            ToOutput << '|' << std::setw(3) << (nd_i - 1);
        }

        ToOutput << '|' << std::setw(5) << NumToStr(this_DATA_percentage * 100, 1) << "%|";
        ToOutput << std::string(bits_filled, bits_filled_wave) << std::string(bits_empty, bar_length - bits_filled_wave) << '|';
        ToOutput << std::setw(5) << NumToStr(this_BTRD_percentage * 100, 1) << "%|";
        ToOutput << std::string(bits_filled, bits_filled_btrd) << std::string(bits_empty, bar_length - bits_filled_btrd) << '|';
        ToOutput << std::setw(5) << NumToStr(this_CORR_percentage * 100, 1) << "%|";
        ToOutput << std::string(bits_filled, bits_filled_corr) << std::string(bits_empty, bar_length - bits_filled_corr) << '|' << std::endl;
    }

    ToOutput << mid_bar << std::endl;
}


void Write_SampleDist(std::ostream& ToOutput)
{
    int32_t nd_i, nd_j;
    int32_t bits_filled_wave;
    int32_t bits_filled_btrd;
    int32_t bits_filled_corr;
    double this_DATA_percentage;
    double this_BTRD_percentage;
    double this_CORR_percentage;

    Titles[0] = "|Input lsb distribution";
    Titles[1] = "|Lossy lsb distribution";
    Titles[2] = "|LWCDF lsb distribution";
    Header = "";

    Make_Bars(Display_Width, 3, 6, 3);
    ToOutput << std::endl << "Sample least significant bit (lsb) distribution." << std::endl << top_bar
             << std::endl << "|Bit" << Titles[0] << Titles[1] << Titles[2] << '|' << std::endl << mid_bar << std::endl;

    double total_samples_processed_recip = OneOver[Global.Channels] / Global.samples_processed;

    for (nd_i = 0; nd_i <= Global.bits_per_sample; nd_i ++)
    {
        this_DATA_percentage = 0;
        this_BTRD_percentage = 0;
        this_CORR_percentage = 0;

        for (nd_j = 0; nd_j < Global.Channels; nd_j ++)
        {
            this_DATA_percentage += DATA_sample_lsb[nd_j][nd_i] * total_samples_processed_recip;
            this_BTRD_percentage += BTRD_sample_lsb[nd_j][nd_i] * total_samples_processed_recip;
            this_CORR_percentage += CORR_sample_lsb[nd_j][nd_i] * total_samples_processed_recip;
        }

        bits_filled_wave = nRoundEvenInt32(std::min(1.0,this_DATA_percentage) * bar_length);
        bits_filled_btrd = nRoundEvenInt32(std::min(1.0,this_BTRD_percentage) * bar_length);
        bits_filled_corr = nRoundEvenInt32(std::min(1.0,this_CORR_percentage) * bar_length);

        if (nd_i == 0)
        {
            ToOutput << "|NUL";
        }
        else
        {
            ToOutput << '|' << std::setw(3) << (nd_i - 1);
        }

        ToOutput << '|' << std::setw(5) << NumToStr(this_DATA_percentage * 100, 1) << "%|";
        ToOutput << std::string(bits_filled, bits_filled_wave) << std::string(bits_empty, bar_length - bits_filled_wave) << '|';
        ToOutput << std::setw(5) << NumToStr(this_BTRD_percentage * 100, 1) << "%|";
        ToOutput << std::string(bits_filled, bits_filled_btrd) << std::string(bits_empty, bar_length - bits_filled_btrd) << '|';
        ToOutput << std::setw(5) << NumToStr(this_CORR_percentage * 100, 1) << "%|";
        ToOutput << std::string(bits_filled, bits_filled_corr) << std::string(bits_empty, bar_length - bits_filled_corr) << '|' << std::endl;
    }

    Titles[0] = "|Input msb distribution";
    Titles[1] = "|Lossy msb distribution";
    Titles[2] = "|LWCDF msb distribution";
    Header = "";

    Make_Bars(Display_Width, 3, 6, 3);
    ToOutput << mid_bar << std::endl << std::endl << "Sample most significant bit (msb) distribution." << std::endl
             << top_bar << std::endl << "|Bit" << Titles[0] << Titles[1] << Titles[2] << '|' << std::endl
             << mid_bar << std::endl;

    for (nd_i = 0; nd_i <= Global.bits_per_sample; nd_i ++)
    {
        this_DATA_percentage = 0;
        this_BTRD_percentage = 0;
        this_CORR_percentage = 0;

        for (nd_j = 0; nd_j < Global.Channels; nd_j ++)
        {
            this_DATA_percentage += DATA_sample_msb[nd_j][nd_i] * total_samples_processed_recip;
            this_BTRD_percentage += BTRD_sample_msb[nd_j][nd_i] * total_samples_processed_recip;
            this_CORR_percentage += CORR_sample_msb[nd_j][nd_i] * total_samples_processed_recip;
        }

        bits_filled_wave = nRoundEvenInt32(std::min(1.0,this_DATA_percentage) * bar_length);
        bits_filled_btrd = nRoundEvenInt32(std::min(1.0,this_BTRD_percentage) * bar_length);
        bits_filled_corr = nRoundEvenInt32(std::min(1.0,this_CORR_percentage) * bar_length);

        if (nd_i == 0)
        {
            ToOutput << "|NUL";
        }
        else
        {
            ToOutput << '|' << std::setw(3) << (nd_i - 1);
        }

        ToOutput << '|' << std::setw(5) << NumToStr(this_DATA_percentage * 100, 1) << "%|";
        ToOutput << std::string(bits_filled, bits_filled_wave) << std::string(bits_empty, bar_length - bits_filled_wave) << '|';
        ToOutput << std::setw(5) << NumToStr(this_BTRD_percentage * 100, 1) << "%|";
        ToOutput << std::string(bits_filled, bits_filled_btrd) << std::string(bits_empty, bar_length - bits_filled_btrd) << '|';
        ToOutput << std::setw(5) << NumToStr(this_CORR_percentage * 100, 1) << "%|";
        ToOutput << std::string(bits_filled, bits_filled_corr) << std::string(bits_empty, bar_length - bits_filled_corr) << '|' << std::endl;
    }

    ToOutput << mid_bar << std::endl;
}


void Write_Histogram(std::ostream& ToOutput)
{
    int32_t nd_i;
    int32_t total_samples;
    int32_t temp_val;

    int32_t DATA_Max;
    int32_t BTRD_Max;
    int32_t CORR_Max;

    double DATA_val;
    double BTRD_val;
    double CORR_val;

    int32_t DATA_len;
    int32_t BTRD_len;
    int32_t CORR_len;

    std::string DATA_bar;
    std::string BTRD_bar;
    std::string CORR_bar;

    Titles[0] = "|Input Value Histogram";
    Titles[1] = "|Lossy Value Histogram";
    Titles[2] = "|LWCDF Value Histogram";
    Header = "";

    Make_Bars(Display_Width, 3, 8, 3);

    ToOutput << std::endl << "Sample Value Histogram [log scale]." << std::endl << top_bar
             << std::endl << "|Bin" << Titles[0] << Titles[1] << Titles[2] << '|' << std::endl << mid_bar << std::endl;

    total_samples = 0;
    DATA_Max = 0;
    BTRD_Max = 0;
    CORR_Max = 0;

    for (nd_i = 0; nd_i <= history.Histogram_Length; nd_i ++)
    {
        temp_val = history.Histogram_DATA[nd_i];
        total_samples = total_samples + temp_val;

        if (temp_val > DATA_Max)
        {
            DATA_Max = temp_val;
        }

        temp_val = history.Histogram_BTRD[nd_i];

        if (temp_val > BTRD_Max)
        {
            BTRD_Max = temp_val;
        }

        temp_val = history.Histogram_CORR[nd_i];

        if (temp_val > CORR_Max)
        {
            CORR_Max = temp_val;
        }
    }

    for (nd_i = 0; nd_i <= history.Histogram_Length; nd_i ++)
    {
        temp_val = history.Histogram_DATA[nd_i];
        DATA_val = double(temp_val) / total_samples * 100;
        DATA_len = int(nlog2(temp_val) / nlog2(DATA_Max) * bar_length);
        DATA_bar = std::string(bits_filled, DATA_len) + std::string(bits_empty, bar_length - DATA_len);

        temp_val = history.Histogram_BTRD[nd_i];
        BTRD_val = double(temp_val) / total_samples * 100;
        BTRD_len = int(nlog2(temp_val) / nlog2(BTRD_Max) * bar_length);
        BTRD_bar = std::string(bits_filled, BTRD_len) + std::string(bits_empty, bar_length - BTRD_len);

        temp_val = history.Histogram_CORR[nd_i];
        CORR_val = double(temp_val) / total_samples * 100;
        CORR_len = int(nlog2(temp_val) / nlog2(CORR_Max) * bar_length);
        CORR_bar = std::string(bits_filled, CORR_len) + std::string(bits_empty, bar_length - CORR_len);
        ToOutput << '|' << std::setw(3) << (nd_i - history.Histogram_Offset)
                 << '|' << std::fixed << std::setprecision(3) << std::setw(7) << DATA_val << "%|" << DATA_bar
                 << '|' << std::fixed << std::setprecision(3) << std::setw(7) << BTRD_val << "%|" << BTRD_bar
                 << '|' << std::fixed << std::setprecision(3) << std::setw(7) << CORR_val << "%|" << CORR_bar
                 << '|' << std::endl;
    }

    ToOutput << mid_bar << std::endl;
}


void Write_FreqDist(std::ostream& ToOutput, int32_t this_analysis_number)
{
    int32_t nt_i, nt_j;
    double nt_x, nt_m, nt_r;
    int32_t fr_l, fr_h, fr_m;

    std::string inp_DATA_string;
    std::string Out_BTRD_String;
    std::string Out_CORR_String;

    double nt_inp_DATA_val;
    double nt_Out_BTRD_val;
    double nt_Out_CORR_val;

    int32_t nt_inp_DATA_len;
    int32_t nt_Out_BTRD_len;
    int32_t nt_Out_CORR_len;

    char BinChar;

    if ((!parameters.output.postanalyse) && (Global.Channels == 2))
    {
        nt_x = Global.blocks_processed_recip;
    }
    else
    {
        nt_x = OneOver[Global.Channels] * Global.blocks_processed_recip;
    }

    nt_m = -log10_2x20 * (Global.bits_per_sample + 1.5 + nlog2(1.5) * 0.50f);
    nt_r = 1.0 / nt_m;

    fr_l = std::max(1, nRoundEvenInt32(double(Global.lower_freq_limit) / Global.sample_rate * settings.analysis[this_analysis_number].FFT.length));
    fr_h = nRoundEvenInt32(double(Global.upper_freq_limit) / Global.sample_rate * settings.analysis[this_analysis_number].FFT.length);
    fr_m = nRoundEvenInt32(20000.0 / Global.sample_rate * settings.analysis[this_analysis_number].FFT.length);


    if ((settings.analysis[this_analysis_number].active) && ((this_analysis_number == SHORT_ANALYSIS) || (parameters.output.longdist)))
    {
        ToOutput << std::endl << "Frequency Analysis of audio data.\n";

        if (parameters.output.postanalyse)
        {
            Titles[0] = std::string("| Input (dBFS)");
            Titles[1] = std::string("| Lossy (dBFS)");
            Titles[2] = std::string("| LWCDF (dBFS)");
            Header = "";

            Make_Bars(Display_Width, 4, 7, 3, 8);

            ToOutput << top_bar << std::endl << "| Bin" << Titles[0] << Titles[1] << Titles[2] << "| Delta  |" << std::endl << mid_bar << std::endl;

            for (nt_i = 0; nt_i <=  (settings.analysis[this_analysis_number].FFT.length/2); nt_i ++)
            {
                nt_inp_DATA_val = 0;
                nt_Out_BTRD_val = 0;
                nt_Out_CORR_val = 0;

                for (nt_j = 0; nt_j < Global.Channels; nt_j ++)
                {
                    nt_inp_DATA_val = nt_inp_DATA_val + results.WAVE[this_analysis_number][nt_j].History[nt_i];
                    nt_Out_BTRD_val = nt_Out_BTRD_val + results.BTRD[this_analysis_number][nt_j].History[nt_i];
                    nt_Out_CORR_val = nt_Out_CORR_val + results.CORR[this_analysis_number][nt_j].History[nt_i];
                }

                nt_inp_DATA_val = std::max(nt_m, log10_2x20 * (nlog2(nt_inp_DATA_val * nt_x) * 0.50f - settings.analysis[this_analysis_number].FFT.bit_length + 3 - Global.bits_per_sample));
                nt_Out_BTRD_val = std::max(nt_m, log10_2x20 * (nlog2(nt_Out_BTRD_val * nt_x) * 0.50f - settings.analysis[this_analysis_number].FFT.bit_length + 3 - Global.bits_per_sample));
                nt_Out_CORR_val = std::max(nt_m, log10_2x20 * (nlog2(nt_Out_CORR_val * nt_x) * 0.50f - settings.analysis[this_analysis_number].FFT.bit_length + 3 - Global.bits_per_sample));

                nt_inp_DATA_len = bar_length - std::min(bar_length, std::max(0, nRoundEvenInt32(bar_length * nt_inp_DATA_val * nt_r)));
                nt_Out_BTRD_len = bar_length - std::min(bar_length, std::max(0, nRoundEvenInt32(bar_length * nt_Out_BTRD_val * nt_r)));
                nt_Out_CORR_len = bar_length - std::min(bar_length, std::max(0, nRoundEvenInt32(bar_length * nt_Out_CORR_val * nt_r)));

                inp_DATA_string = std::string(bits_filled, nt_inp_DATA_len) + std::string(bits_empty, bar_length - nt_inp_DATA_len);
                Out_BTRD_String = std::string(bits_filled, nt_Out_BTRD_len) + std::string(bits_empty, bar_length - nt_Out_BTRD_len);
                Out_CORR_String = std::string(bits_filled, nt_Out_CORR_len) + std::string(bits_empty, bar_length - nt_Out_CORR_len);

                if (nt_i == fr_l)
                {
                    BinChar = 'L';
                }
                else if (nt_i == fr_h)
                {
                    BinChar = 'U';
                }
                else if (nt_i == fr_m)
                {
                    BinChar = '*';
                }
                else
                {
                    BinChar = ' ';
                }

                if (nt_i == 0)
                {
                    ToOutput << "| DC ";
                }
                else if (nt_i == PowersOf.TwoInt64[settings.analysis[this_analysis_number].FFT.bit_length - 1])
                {
                    ToOutput << "|FS/2";
                }
                else
                {
                    ToOutput << '|' << BinChar << std::setw(3) << nt_i;
                }

                ToOutput << '|' << std::setw(7) << std::fixed << std::setprecision(2) << nt_inp_DATA_val << '|' << inp_DATA_string << '|'
                                 << std::setw(7) << std::fixed << std::setprecision(2) << nt_Out_BTRD_val << '|' << Out_BTRD_String << '|'
                                 << std::setw(7) << std::fixed << std::setprecision(2) << nt_Out_CORR_val << '|' << Out_CORR_String << '|'
                                 << BinChar << std::setw(7) << std::fixed << std::setprecision(3) << nround10(nt_Out_BTRD_val - nt_inp_DATA_val, 3) << '|' << std::endl;
            }

            ToOutput << mid_bar << std::endl;
        }
        else if (Global.Channels != 2)
        {
            Titles[0] = "| Input Average (dBFS)";
            Titles[1] = "";
            Titles[2] = "";
            Header = "";

            Make_Bars(Display_Width, 4, 7, 1);
            ToOutput << top_bar << std::endl << "| Bin" << Titles[0] << std::string("|") << std::endl << mid_bar << std::endl;

            for (nt_i = 0; nt_i <= settings.analysis[this_analysis_number].FFT.length/2; nt_i ++)
            {
                nt_inp_DATA_val = 0;

                for (nt_j = 0; nt_j < Global.Channels; nt_j ++)
                {
                    nt_inp_DATA_val = nt_inp_DATA_val + results.WAVE[this_analysis_number][nt_j].History[nt_i];
                }

                nt_inp_DATA_val = std::max(nt_m, log10_2x20 * (nlog2(nt_inp_DATA_val * nt_x) * 0.50f - settings.analysis[this_analysis_number].FFT.bit_length + 3 - Global.bits_per_sample));
                nt_inp_DATA_len = bar_length - std::min(bar_length, std::max(0, nRoundEvenInt32(bar_length * nt_inp_DATA_val * nt_r)));
                inp_DATA_string = std::string(bits_filled, nt_inp_DATA_len) + std::string(bits_empty, bar_length - nt_inp_DATA_len);

                if (nt_i == fr_l)
                {
                    BinChar = 'L';
                }
                else if (nt_i == fr_h)
                {
                    BinChar = 'U';
                }
                else if (nt_i == fr_m)
                {
                    BinChar = '*';
                }
                else
                {
                    BinChar = ' ';
                }

                if (nt_i == 0)
                {
                    ToOutput << "| DC ";
                }
                else if (nt_i == PowersOf.TwoInt64[settings.analysis[this_analysis_number].FFT.bit_length - 1])
                {
                    ToOutput << "|FS/2";
                }
                else
                {
                    ToOutput << '|' << BinChar << std::setw(3) << nt_i;
                }

                ToOutput << '|' << std::setw(7) << std::fixed << std::setprecision(2) << nt_inp_DATA_val << '|' << inp_DATA_string << '|' << std::endl;
            }

            ToOutput << mid_bar << std::endl;
        }
        else
        {
            Titles[0] = "| Input Average (dBFS) Channel #0";
            Titles[1] = "| Input Average (dBFS) Channel #1";
            Titles[2] = "";
            Header = "";

            Make_Bars(Display_Width, 4, 7, 2);
            ToOutput << top_bar << std::endl << "| Bin" << Titles[0] << Titles[1] << std::string("|") << std::endl << mid_bar << std::endl;

            for (nt_i = 0; nt_i <= settings.analysis[this_analysis_number].FFT.length/2; nt_i ++)
            {
                nt_inp_DATA_val = results.WAVE[this_analysis_number][0].History[nt_i];
                nt_Out_BTRD_val = results.WAVE[this_analysis_number][1].History[nt_i];
                nt_inp_DATA_val = std::max(nt_m, log10_2x20 * (nlog2(nt_inp_DATA_val * nt_x) * 0.50f - settings.analysis[this_analysis_number].FFT.bit_length + 3 - Global.bits_per_sample));
                nt_Out_BTRD_val = std::max(nt_m, log10_2x20 * (nlog2(nt_Out_BTRD_val * nt_x) * 0.50f - settings.analysis[this_analysis_number].FFT.bit_length + 3 - Global.bits_per_sample));
                nt_inp_DATA_len = bar_length - std::min(bar_length, std::max(0, nRoundEvenInt32(bar_length * nt_inp_DATA_val * nt_r)));
                nt_Out_BTRD_len = bar_length - std::min(bar_length, std::max(0, nRoundEvenInt32(bar_length * nt_Out_BTRD_val * nt_r)));
                inp_DATA_string = std::string(bits_filled, nt_inp_DATA_len) + std::string(bits_empty, bar_length - nt_inp_DATA_len);
                Out_BTRD_String = std::string(bits_filled, nt_Out_BTRD_len) + std::string(bits_empty, bar_length - nt_Out_BTRD_len);

                if (nt_i == fr_l)
                {
                    BinChar = 'L';
                }
                else if (nt_i == fr_h)
                {
                    BinChar = 'U';
                }
                else if (nt_i == fr_m)
                {
                    BinChar = '*';
                }
                else
                {
                    BinChar = ' ';
                }

                if (nt_i == 0)
                {
                    ToOutput << "| DC ";
                }
                else if (nt_i == PowersOf.TwoInt64[settings.analysis[this_analysis_number].FFT.bit_length - 1])
                {
                    ToOutput << "|FS/2";
                }
                else
                {
                    ToOutput << '|' << BinChar << std::setw(3) << nt_i;
                }

                ToOutput << '|' << std::setw(7) << std::fixed << std::setprecision(2) << nt_inp_DATA_val << '|' << inp_DATA_string << '|'
                                << std::setw(7) << std::fixed << std::setprecision(2) << nt_Out_BTRD_val << '|' << Out_BTRD_String << '|' << std::endl;
            }

            ToOutput << mid_bar << std::endl;
        }

        ToOutput << "Legend: L = Lower Frequency Calculation Limit (" << std::setw(6) << NumToStr(Global.lower_freq_limit * OneOver[1000], 3) << "kHz);" << std::endl
                 << "        U = Upper Frequency Calculation Limit (" << std::setw(6) << NumToStr(Global.upper_freq_limit * OneOver[1000], 3) << "kHz);" << std::endl
                 << "        * = 20.000kHz\n";
    }
}

void Write_Min_Bin_Dist(std::ostream& ToOutput)
{
    int32_t nt_i, nt_j;
    int32_t Bar_Full_Len;
    double Old_Percent;
    double New_Percent;
    double Alt_Percent;
    double Tot_Old_Percent;
    double Tot_New_Percent;
    double Tot_Alt_Percent;
    std::string Old_Bar_Str;
    std::string New_Bar_Str;
    std::string Alt_Bar_Str;
    uint64_t Local_Max;
    double Percent_Used [9];
    int32_t fr_l, fr_h, fr_m;
    char BinChar;

    if (parameters.output.spread == 2)
    {
        Titles[0] = "| Old Minimum Value Used";
        Titles[1] = "| New Minimum Value Used";
        Titles[2] = "| Minimum Value Used";
        Header = "";
        Make_Bars(Display_Width, 4, 7, 3);

        for (nt_j = 1; nt_j <= PRECALC_ANALYSES; nt_j ++)
        {
            if (settings.analysis[nt_j].active)
            {
                Local_Max = 0;
                Tot_Old_Percent = 0;
                Tot_New_Percent = 0;
                Tot_Alt_Percent = 0;

                for (nt_i = 0; nt_i <= spreading.Bins.Upper[Current.Analysis.bits[nt_j]] + 1; nt_i ++)
                {
                    Local_Max = std::max(Local_Max, process.Old_Min_Used_History[nt_j][nt_i] + process.New_Min_Used_History[nt_j][nt_i]);
                }

                fr_l = spreading.Bins.Lower[Current.Analysis.bits[nt_j]];
                fr_h = spreading.Bins.Upper[Current.Analysis.bits[nt_j]];
                fr_m = nRoundEvenInt32(20000.0 / Global.sample_rate * Current.Analysis.length[nt_j]);
                ToOutput << std::endl << "Spreading Algorithm Results : Minima, FFT Length = " << NumToStr(Current.Analysis.length[nt_j])
                         << std::endl << top_bar << std::endl << "|Bin " << Titles[0] << Titles[1] << Titles[2] << std::string("|") << std::endl << mid_bar << std::endl;

                for (nt_i = 0; nt_i <= spreading.Bins.Upper[Current.Analysis.bits[nt_j]] + 1; nt_i ++)
                {
                    if (nt_i == fr_l)
                    {
                        BinChar = 'L';
                    }
                    else if (nt_i == fr_h)
                    {
                        BinChar = 'U';
                    }
                    else if (nt_i == fr_m)
                    {
                        BinChar = '*';
                    }
                    else
                    {
                        BinChar = ' ';
                    }

                    if (nt_i == 0)
                    {
                        ToOutput << "| DC ";
                    }
                    else if (nt_i == Current.Analysis.length[nt_j])
                    {
                        ToOutput << "|FS/2";
                    }
                    else
                    {
                        ToOutput << '|' << BinChar << std::setw(3) << nt_i;
                    }

                    Old_Percent = double(process.Old_Min_Used_History[nt_j][nt_i]) / process.Analyses_Completed[nt_j] * 100;
                    New_Percent = double(process.New_Min_Used_History[nt_j][nt_i]) / process.Analyses_Completed[nt_j] * 100;
                    Alt_Percent = Old_Percent + New_Percent;
                    Bar_Full_Len = nRoundEvenInt32(double(process.Old_Min_Used_History[nt_j][nt_i]) / Local_Max * bar_length);
                    Old_Bar_Str = std::string(bits_filled, Bar_Full_Len) + std::string(bits_empty, bar_length - Bar_Full_Len);
                    Bar_Full_Len = nRoundEvenInt32(double(process.New_Min_Used_History[nt_j][nt_i]) / Local_Max * bar_length);
                    New_Bar_Str = std::string(bits_filled, Bar_Full_Len) + std::string(bits_empty, bar_length - Bar_Full_Len);
                    Bar_Full_Len = nRoundEvenInt32(double(process.Old_Min_Used_History[nt_j][nt_i] + process.New_Min_Used_History[nt_j][nt_i]) / Local_Max * bar_length);
                    Alt_Bar_Str = std::string(bits_filled, Bar_Full_Len) + std::string(bits_empty, bar_length - Bar_Full_Len);
                    Tot_Old_Percent = Tot_Old_Percent + Old_Percent;
                    Tot_New_Percent = Tot_New_Percent + New_Percent;
                    Tot_Alt_Percent = Tot_Alt_Percent + Alt_Percent;
                    ToOutput << '|' << std::fixed << std::setw(6) << std::setprecision(2) << Old_Percent << "%|" << Old_Bar_Str
                             << '|' << std::fixed << std::setw(6) << std::setprecision(2) << New_Percent << "%|" << New_Bar_Str
                             << '|' << std::fixed << std::setw(6) << std::setprecision(2) << Alt_Percent << "%|" << Alt_Bar_Str << '|' << std::endl;
                }

                ToOutput << mid_bar << std::endl;
                Bar_Full_Len = nRoundEvenInt32(Tot_Old_Percent * bar_length * OneOver[100]);
                Old_Bar_Str = std::string(bits_filled, Bar_Full_Len) + std::string(bits_empty, bar_length - Bar_Full_Len);
                Bar_Full_Len = nRoundEvenInt32(Tot_New_Percent * bar_length * OneOver[100]);
                New_Bar_Str = std::string(bits_filled, Bar_Full_Len) + std::string(bits_empty, bar_length - Bar_Full_Len);
                Bar_Full_Len = nRoundEvenInt32(Tot_Alt_Percent * bar_length * OneOver[100]);
                Alt_Bar_Str = std::string(bits_filled, Bar_Full_Len) + std::string(bits_empty, bar_length - Bar_Full_Len);
                ToOutput << "|Tot.|" << std::fixed << std::setw(6) << std::setprecision(2) << Tot_Old_Percent << "%|" << Old_Bar_Str
                         <<      '|' << std::fixed << std::setw(6) << std::setprecision(2) << Tot_New_Percent << "%|" << New_Bar_Str
                         <<      '|' << std::fixed << std::setw(6) << std::setprecision(2) << Tot_Alt_Percent << "%|" << Alt_Bar_Str << '|' << std::endl;
                ToOutput << mid_bar << std::endl;
            }
        }
    }

    Titles[0] = "| Minimum used per codec-block";
    Titles[1] = "";
    Titles[2] = "";
    Header = "";
    Make_Bars(Display_Width, 4, 7, 1);
    ToOutput << std::endl << "Minimum used to determine bits-to-remove, by FFT length:" << std::endl << top_bar;
    ToOutput << std::endl << "|FFT " << Titles[0] << std::string("|") << std::endl << mid_bar << std::endl;

    for (nt_j = 1; nt_j <= PRECALC_ANALYSES; nt_j ++)
        if (settings.analysis[nt_j].active)
        {
            Percent_Used[nt_j] = 0;

            for (nt_i = 0; nt_i < Global.Channels; nt_i ++)
            {
                Percent_Used[nt_j] = Percent_Used[nt_j] + results.minima[nt_i][nt_j];
            }

            Percent_Used[nt_j] = Percent_Used[nt_j] * Global.blocks_processed_recip * OneOver[Global.Channels] * 100;
            Bar_Full_Len = nRoundEvenInt32(Percent_Used[nt_j] * bar_length * OneOver[100]);
            Old_Bar_Str = std::string(bits_filled, Bar_Full_Len) + std::string(bits_empty, bar_length - Bar_Full_Len);
            ToOutput << '|' << std::setw(4) << Current.Analysis.length[nt_j]
                     << '|' << std::setw(6) << std::fixed << std::setprecision(2) << Percent_Used[nt_j]
                     << "%|" << Old_Bar_Str << '|' << std::endl;
        }

    nt_j = 7;
    Percent_Used[nt_j] = 0;

    for (nt_i = 0; nt_i < Global.Channels; nt_i ++)
    {
        Percent_Used[nt_j] += results.minima[nt_i][nt_j];
    }

    Percent_Used[nt_j] = Percent_Used[nt_j] * Global.blocks_processed_recip * OneOver[Global.Channels] * 100;
    Bar_Full_Len = nRoundEvenInt32(Percent_Used[nt_j] * bar_length * OneOver[100]);
    Old_Bar_Str = std::string(bits_filled, Bar_Full_Len) + std::string(bits_empty, bar_length - Bar_Full_Len);
    ToOutput << "|Stat|" << std::setw(6) << std::fixed << std::setprecision(2) << Percent_Used[nt_j] << "%|" << Old_Bar_Str << '|' << std::endl;
    nt_j = 8;
    Percent_Used[nt_j] = 0;

    for (nt_i = 0; nt_i < Global.Channels; nt_i ++)
    {
        Percent_Used[nt_j] = Percent_Used[nt_j] + results.minima[nt_i][nt_j];
    }

    Percent_Used[nt_j] = Percent_Used[nt_j] * Global.blocks_processed_recip * OneOver[Global.Channels] * 100;
    Bar_Full_Len = nRoundEvenInt32(Percent_Used[nt_j] * bar_length * OneOver[100]);
    Old_Bar_Str = std::string(bits_filled, Bar_Full_Len) + std::string(bits_empty, bar_length - Bar_Full_Len);
    ToOutput << "|Dyn.|" << std::setw(6) << std::fixed << std::setprecision(2) << Percent_Used[nt_j] << "%|" << Old_Bar_Str << '|' << std::endl;
    ToOutput << mid_bar << std::endl;

    Titles[0] = "| Old Minimum Value used";
    Titles[1] = "| New Minimum Value used";
    Titles[2] = "| Average Value used";
    Header = "";

    Make_Bars(Display_Width, 4, 7, 3);
    ToOutput << std::endl << "Spreading Algorithm Results : Minima, per analysis" << std::endl << top_bar;
    ToOutput << std::endl << "|FFT " << Titles[0] << Titles[1] << Titles[2] << "|" << std::endl << mid_bar << std::endl;

    for (nt_j = 1; nt_j <= PRECALC_ANALYSES; nt_j ++)
        if (settings.analysis[nt_j].active)
        {
            Old_Percent = double(process.Old_Min_Used[nt_j]) / process.Analyses_Completed[nt_j] * 100;
            New_Percent = double(process.New_Min_Used[nt_j]) / process.Analyses_Completed[nt_j] * 100;
            Alt_Percent = double(process.Alt_Ave_Used[nt_j]) / process.Analyses_Completed[nt_j] * 100;
            Bar_Full_Len = nRoundEvenInt32(Old_Percent * bar_length * OneOver[100]);
            Old_Bar_Str = std::string(bits_filled, Bar_Full_Len) + std::string(bits_empty, bar_length - Bar_Full_Len);
            Bar_Full_Len = nRoundEvenInt32(New_Percent * bar_length * OneOver[100]);
            New_Bar_Str = std::string(bits_filled, Bar_Full_Len) + std::string(bits_empty, bar_length - Bar_Full_Len);
            Bar_Full_Len = nRoundEvenInt32(Alt_Percent * bar_length * OneOver[100]);
            Alt_Bar_Str = std::string(bits_filled, Bar_Full_Len) + std::string(bits_empty, bar_length - Bar_Full_Len);
            ToOutput << '|' << std::setw(4) << Current.Analysis.length[nt_j]
                     << '|' << std::fixed << std::setw(6) << std::setprecision(2) << Old_Percent << "%|" << Old_Bar_Str
                     << '|' << std::fixed << std::setw(6) << std::setprecision(2) << New_Percent << "%|" << New_Bar_Str
                     << '|' << std::fixed << std::setw(6) << std::setprecision(2) << Alt_Percent << "%|" << Alt_Bar_Str
                     << '|' << std::endl;
        }

    ToOutput << mid_bar << std::endl;
    Titles[0] = "| BTR > Static-Max-Bits-To-Remove";
    Titles[1] = "| BTR > Dynamic-Max-Bits-To-Remove";
    Titles[2] = "";
    Header = "";

    Make_Bars(Display_Width, 4, 7, 2);
    ToOutput << std::endl << "Minimum_Bits_To_Keep Results (BTR limited to relevant maximum)" << std::endl << top_bar
             << std::endl << "|FFT " << Titles[0] << Titles[1] << '|' << std::endl << mid_bar << std::endl;

    for (nt_j = 1; nt_j <= PRECALC_ANALYSES; nt_j ++)
        if (settings.analysis[nt_j].active)
        {
            Old_Percent = double(process.Over_Static[nt_j]) / process.Analyses_Completed[nt_j] * 100;
            New_Percent = double(process.Over_Dynamic[nt_j]) / process.Analyses_Completed[nt_j] * 100;
            Bar_Full_Len = nRoundEvenInt32(Old_Percent * bar_length * OneOver[100]);
            Old_Bar_Str = std::string(bits_filled, Bar_Full_Len) + std::string(bits_empty, bar_length - Bar_Full_Len);
            Bar_Full_Len = nRoundEvenInt32(New_Percent * bar_length * OneOver[100]);
            New_Bar_Str = std::string(bits_filled, Bar_Full_Len) + std::string(bits_empty, bar_length - Bar_Full_Len);
            ToOutput << '|' << std::setw(4) << Current.Analysis.length[nt_j]
                     << '|' << std::fixed << std::setw(6) << std::setprecision(2) << Old_Percent << "%|" << Old_Bar_Str
                     << '|' << std::fixed << std::setw(6) << std::setprecision(2) << New_Percent << "%|" << New_Bar_Str
                     << '|' << std::endl;
        }

    ToOutput << mid_bar << std::endl;
}


void Process_Output()
{
    if ((parameters.output.verbosity) && (!parameters.output.silent) && ((Global.last_print == 0) || (Global.blocks_processed % Global.output_blocks == 0)))
    {
        gettimer();
        time_string_make(strings.Elapsed, timer.Elapsed);
        std::cerr << '\r' << "Progress  :";

        if (Global.WAVE_size == MAX_WAVE_SIZE)
        {
            size_string_make(strings.Size, 1.0 * Global.samples_processed * Global.Channels * Global.bytes_per_sample);
            std::cerr << ' ' << strings.Size << ';';
        }
        else
        {
            std::cerr << std::fixed << std::setw(6) << std::setprecision(2) << (double(Global.samples_processed) / Global.Total_Samples * 100) << "%;";
        }

        std::cerr << std::fixed << std::setw(7) << std::setprecision(4) << (OneOver[Global.Channels] * Stats.total_bits_removed * Global.blocks_processed_recip) << " bits;"
                  << std::fixed << std::setw(6) << std::setprecision(2) << Global.processing_rate << "x;";
        std::cerr << ' ' << strings.Elapsed;

        if ((Global.WAVE_size != MAX_WAVE_SIZE) && (Global.samples_processed > 0))
        {
            time_string_make(strings.estimate, std::min(86399.0, double(Global.Total_Samples) / Global.samples_processed * timer.Elapsed));
            std::cerr << '/' << strings.estimate << ' ';
        }

        Global.last_print = Global.blocks_processed;
    }
}

std::string SChar(int32_t SVal)
{
    if (SVal != 1)
        return "s";
    else
        return "";
}

void Write_Results(std::ostream& ToOutput)
{
    double bits_lost;

    ToOutput << "Results   : " << std::fixed << std::setprecision(4) << (OneOver[Global.Channels] * Stats.total_bits_removed * Global.blocks_processed_recip)
             << " bits; " << std::fixed << std::setw(0) << std::setprecision(2) << Global.processing_rate << "x; " << strings.Elapsed;

    if (FFTW_Initialised())
    {
        ToOutput << "; [F]";
    }
    else
    {
        ToOutput << "; [I]";
    }

    ToOutput << std::endl;

    bits_lost = nround10(Stats.total_bits_lost * OneOver[Global.Channels] * Global.blocks_processed_recip, 4);

    ToOutput << "Feedback  :";

    bool output_written = false;

    if (bits_lost > 0.0)
    {
        ToOutput << std::setw(7) << std::fixed << std::setprecision(4) << bits_lost/*: 7: 4*/ << " bits lost;";
        output_written = true;
    }

    if (Stats.Count.eclips > 0)
    {
        ToOutput << " Extant clip" << SChar(Stats.Count.eclips) << ": " << Stats.Count.eclips << ";";
        output_written = true;
    }

    if (Stats.Count.sclips > 0)
    {
        ToOutput << " Scaling clip" << SChar(Stats.Count.sclips) << ": " << Stats.Count.sclips << ";";
        output_written = true;
    }

    if (Stats.Count.rclips > 0)
    {
        ToOutput << " Rounding clip" << SChar(Stats.Count.rclips) << ": " << Stats.Count.rclips << ";";
        output_written = true;
    }

    if (parameters.feedback.active)
    {
        if (Stats.Incidence.round > 0)
        {
            ToOutput << " Rounding noise: " << Stats.Incidence.round << ";";
            output_written = true;
        }

        if (Stats.Incidence.noise > 0)
        {
            ToOutput << " Noise shaping noise: " << Stats.Incidence.noise << ";";
            output_written = true;
        }

        if (Stats.Count.aclips > 0)
        {
            ToOutput << " Noise shaping clip" << SChar(Stats.Count.aclips) << ": " << Stats.Count.aclips << ";";
            output_written = true;
        }

        if (Stats.Count.xclips > 0)
        {
            ToOutput << " Reactive Noise shaping" << SChar(Stats.Count.xclips) << ": " << Stats.Count.xclips << ";";
            output_written = true;
        }

        if ((parameters.shaping.active) && (Stats.Skipped_Filters > 0))
        {
            ToOutput << " ANS Filters Skipped : " << Stats.Skipped_Filters << ";";
            output_written = true;
        }
    }

    if (parameters.feedback.verbose)
    {
        ToOutput << std::endl << "  Clipping  : Incidence: Extant: " << std::setw(6) << Stats.Incidence.eclip << "; Scaling: " << std::setw(6) << Stats.Incidence.sclip << ";";
        ToOutput << " Rounding: " << Stats.Incidence.rclip << "; Shaping: " << std::setw(6) << Stats.Incidence.aclip << ";" << "; Reactive: " << std::setw(6) << Stats.Incidence.xclip << ";";
        ToOutput << std::endl << "  Clipping  : Count    : Extant: " << std::setw(6) << Stats.Count.eclips << "; Scaling: " << std::setw(6) << Stats.Count.sclips << ";";
        ToOutput << " Rounding: " << Stats.Count.rclips << "; Shaping: " << std::setw(6) << Stats.Count.aclips << ";" << "; Reactive: " << std::setw(6) << Stats.Count.xclips << ";";
        ToOutput << std::endl << "  Exceedence: Rounding : " << std::setw(6) << Stats.Incidence.round << "; Shaping: " << std::setw(6) << Stats.Incidence.noise << ";";
        ToOutput << " ANS Filters Skipped : " << Stats.Skipped_Filters << ";";
        output_written = true;
    }

    if (!output_written)
    {
        ToOutput << " None.";
    }

    ToOutput << std::endl;

    if (parameters.output.freqdist)
    {
        for (int32_t this_analysis_number = 1; this_analysis_number < (PRECALC_ANALYSES + 1); this_analysis_number++)
        {
            if (settings.analysis[this_analysis_number].active)
                Write_FreqDist(ToOutput, this_analysis_number);
        }
    }

    if (parameters.output.bitdist)
    {
        Write_bitdist(ToOutput);
    }

    if (parameters.output.blockdist)
    {
        Write_blockdist(ToOutput);
    }

    if (parameters.output.sampledist)
    {
        Write_SampleDist(ToOutput);
    }

    if (parameters.output.detail)
    {
        remove_bits_detailed_output(ToOutput);
    }

    if (parameters.output.spread != -1)
    {
        Write_Min_Bin_Dist(ToOutput);
    }

    if (parameters.output.histogram)
    {
        Write_Histogram(ToOutput);
    }

    flush(ToOutput);
}


void write_cleanup()
{
    gettimer();
    time_string_make(strings.Elapsed, timer.Elapsed);

    if (!parameters.output.silent)
    {
        if (parameters.output.verbosity)
        {
            std::cerr << "\r                                                                               \r";
            Write_Results(std::cerr);
        }
        else
        {
            std::cerr << std::fixed << std::setw(7) << std::setprecision(4) << (OneOver[Global.Channels] * Stats.total_bits_removed * Global.blocks_processed_recip) << ';'
                      << std::fixed << std::setw(6) << std::setprecision(2) << Global.processing_rate << "x; " << strings.Elapsed;

            if (FFTW_Initialised())
            {
                std::cerr << "; [F]";
            }
            else
            {
                std::cerr << "; [I]";
            }

            if ((parameters.shaping.active) && (Stats.Skipped_Filters > 0))
            {
                std::cerr << "; F" << Stats.Skipped_Filters;
            }

            std::cerr << std::endl;
        }
    }

    if (parameters.output.writetolog)
    {
        time_string_make(strings.time, double(Global.samples_processed) / Global.sample_rate);
        size_string_make(strings.Size, 1.0 * Global.samples_processed * Global.Channels * Global.bytes_per_sample);

        open_log_file();

        LogOutput << version_string << strings.version_short << lossyWAVHead1
                  << "Processed : " << strings.datestamp << std::endl
                  << "Settings  : " << strings.parameter << std::endl
                  << "Filename  : " << WAVFilePrintName() << std::endl
                  << "File Info : " << std::fixed << std::setw(0) << std::setprecision(2) << (Global.sample_rate / 1000.0) << "kHz; "
                  << Global.Channels << " channel; " << Global.bits_per_sample << " bit, "
                  << strings.time << ", "
                  << strings.Size << std::endl;

        Write_Results(LogOutput);

        LogOutput << std::endl;

        close_log_file();
    }
}


void nOutput_Init()
{
    int32_t nd_i, nd_j;

    for (nd_j = 0; nd_j < MAX_CHANNELS; nd_j ++)
    {
        for (nd_i = 0; nd_i <= 33; nd_i ++)
        {
            DATA_block_lsb[nd_j][nd_i] = 0;
            DATA_block_msb[nd_j][nd_i] = 0;
            DATA_sample_lsb[nd_j][nd_i] = 0;
            DATA_sample_msb[nd_j][nd_i] = 0;
            BTRD_block_lsb[nd_j][nd_i] = 0;
            BTRD_block_msb[nd_j][nd_i] = 0;
            BTRD_sample_lsb[nd_j][nd_i] = 0;
            BTRD_sample_msb[nd_j][nd_i] = 0;
            CORR_block_lsb[nd_j][nd_i] = 0;
            CORR_block_msb[nd_j][nd_i] = 0;
            CORR_sample_lsb[nd_j][nd_i] = 0;
            CORR_sample_msb[nd_j][nd_i] = 0;
        }
    }

    if (parameters.output.width == -1)
    {
        parameters.output.width = 79;
    }
    Display_Width = parameters.output.width;
}
