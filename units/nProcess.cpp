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
#include "nMaths.h"
#include "nSpreading.h"
#include "fftw_interface.h"
#include "nFFT.h"
#include "nFillFFT.h"
#include "nRemoveBits.h"
#include "nOutput.h"

#include "nProcess.h"

void Fill_Last_from_SGNS()
{
    for (int32_t sa_i = 0; sa_i <= Global.FFT.length_half; ++sa_i)
    {
        double temp_result = SGNSFFT.Last_Results[Global.Channel][sa_i];
        FFT_last_unity_result[sa_i] = temp_result;
        FFT_unity_result[sa_i] += temp_result;
    }
}


void Fill_Last_from_History()
{
    for (int32_t sa_i = 0; sa_i <= Global.FFT.length_half;  ++sa_i)
    {
        double temp_result = history.Last_Results[Global.Channel][sa_i];
        FFT_last_unity_result[sa_i] = temp_result;
        FFT_unity_result[sa_i] += temp_result;
    }
}


void Fill_History_from_Last(FFT_Proc_Rec* this_FFT_plan)
{
    for (int32_t sa_i = 0; sa_i <= this_FFT_plan->FFT->length_half; sa_i++)
    {
        double temp_result = FFT_last_unity_result[sa_i];
        history.Last_Results[Global.Channel][sa_i] = temp_result;

        double divisor = OneOver[this_FFT_plan->Task.analyses_performed-1];

        SGNSFFT.FFT_results_short_root[sa_i] = divisor * (FFT_root_result[sa_i] - FFT_last_root_result[sa_i]);

        temp_result = divisor * (FFT_unity_result[sa_i] - temp_result);
        history.WAVE_results[Global.Channel][sa_i] += temp_result;
        SGNSFFT.FFT_results_short[sa_i] = temp_result;
    }
}


void Fill_SGNSFFT_from_Last(FFT_Proc_Rec* this_FFT_plan)
{
    for (int32_t sa_i = 0; sa_i <= this_FFT_plan->FFT->length_half; ++sa_i)
    {
        double temp_result = FFT_last_unity_result[sa_i];
        SGNSFFT.Last_Results[Global.Channel][sa_i] = temp_result;
        SGNSFFT.FFT_results_long[sa_i] = FFT_unity_result[sa_i] - temp_result;
        SGNSFFT.FFT_results_long_root[sa_i] = FFT_root_result[sa_i];
    }
}


void Fill_Last_with_Zero()
{
    for (int32_t sa_i = 0; sa_i <= Global.FFT.length_half; ++sa_i)
    {
        FFT_last_unity_result[sa_i] = 0;
        FFT_last_root_result[sa_i] = 0;
        FFT_skew_result[sa_i] = 0;
    }
}


void Zero_FFT_unity_results()
{
    for (int32_t sa_i = 0; sa_i <= Global.FFT.length_half; ++sa_i)
    {
        FFT_unity_result[sa_i] = 0;
        FFT_root_result[sa_i] = 0;
    }
}

void Post_Process_FFT_Array_to_Last(FFT_Proc_Rec* this_FFT_plan)
{
    for (int32_t sa_i = 0; sa_i <= this_FFT_plan->FFT->length_half; ++sa_i)
    {
        double sa_x = magsquared(this_FFT_plan->DComplex[sa_i]);
        FFT_last_unity_result[sa_i] = sa_x;
        FFT_unity_result[sa_i] += sa_x;
    }
}


void Post_Process_FFT_Results(FFT_Proc_Rec* this_FFT_plan)
{
    for (int32_t sa_i = 0; sa_i <= this_FFT_plan->FFT->length_half; ++sa_i)
    {
        double sc_y = magnitude(this_FFT_plan->DComplex[sa_i]);

        FFT_last_root_result[sa_i] = sc_y;
        FFT_root_result[sa_i] += sc_y;

        FFT_skew_result[sa_i] = sc_y * Skewing_Gain[sa_i << Global.FFT.bit_shift_from_max];

        sc_y *= sc_y;
        FFT_last_unity_result[sa_i] = sc_y;
        FFT_unity_result[sa_i] += sc_y;
    }
}


void Post_Process_FFT_Results_No_Raw(FFT_Proc_Rec* this_FFT_plan)
{
    for (int32_t sc_i = 0; sc_i <= spreading.Bins.Upper[this_FFT_plan->FFT->bit_length] + 1; ++sc_i)
    {
        double sc_y = magsquared(this_FFT_plan->DComplex[sc_i]);
        FFT_last_unity_result[sc_i] = sc_y;
        FFT_unity_result[sc_i] += sc_y;
        FFT_skew_result[sc_i] = nroot(sc_y) * Skewing_Gain[sc_i << this_FFT_plan->FFT->bit_shift_from_max];
    }
}


void Add_to_Histogram()
{
    for (int32_t sa_j = 0; sa_j < Global.Channels; ++sa_j)
        for (int32_t sa_i = 0; sa_i < AudioData.Size.This; ++sa_i)
        {
            history.Histogram_DATA[nRoundEvenInt32(AudioData.WAVEPTR[THIS_CODEC_BLOCK][sa_j][sa_i].Integers[0] * history.Histogram_Multiplier) + history.Histogram_Offset]++;
            history.Histogram_BTRD[nRoundEvenInt32(AudioData.BTRDPTR[THIS_CODEC_BLOCK][sa_j][sa_i].Integers[0] * history.Histogram_Multiplier) + history.Histogram_Offset]++;
            history.Histogram_CORR[nRoundEvenInt32(AudioData.CORRPTR[THIS_CODEC_BLOCK][sa_j][sa_i].Integers[0] * history.Histogram_Multiplier) + history.Histogram_Offset]++;
        }
}


void Process_Data_to_LongDist(FFT_Proc_Rec* this_FFT_plan)
{
    Zero_FFT_unity_results();

    this_FFT_plan->Task.analyses_performed = 0;

    for (int32_t this_analysis_block_number = 0; this_analysis_block_number < process.analysis_blocks[Global.FFT.bit_length]; ++this_analysis_block_number)
    {
        int32_t this_block_start = floor(process.actual_analysis_blocks_start[this_FFT_plan->FFT->bit_length] + this_analysis_block_number * process.FFT_underlap_length[this_FFT_plan->FFT->bit_length]);
        this_FFT_plan->Task.block_start = std::min(std::max(this_block_start, process.limits.minstart),process.limits.maxend-this_FFT_plan->FFT->length);
        this_FFT_plan->Task.analyses_performed++;

        if (this_FFT_plan->Task.Fill_FFT_Proc(this_FFT_plan) != 0)
        {
            if (FFTW_Initialised())
                FFTW.Execute_R2C_New_Array(FFTW.Plans[this_FFT_plan->FFT->bit_length],&this_FFT_plan->DReal[0],&this_FFT_plan->DReal[0]);
            else
                FFT_DIT_Real(this_FFT_plan);

            Post_Process_FFT_Array_to_Last(this_FFT_plan);
        }
    }

    for (int32_t sa_i = 0; sa_i <= this_FFT_plan->FFT->length_half; ++sa_i)
        this_FFT_plan->Task.Dist_Array[sa_i] += (FFT_unity_result[sa_i] * OneOver[this_FFT_plan->Task.analyses_performed]);
}


void Process_This_Codec_Block()
{
    int32_t this_channel;
    int32_t this_analysis_number;
    int32_t this_analysis_block_number;
    int32_t codec_block_dependent_bits_to_remove;
    int32_t local_channels;
    double Spreading_result;
    double bits_removed_this_codec_block;
    int32_t Spreading_Used;

    process.limits.minstart = -(AudioData.Size.Prev+AudioData.Size.Last);
    process.limits.maxend = (AudioData.Size.This+AudioData.Size.Next);

    codec_block_dependent_bits_to_remove = Global.bits_per_sample;
    bits_removed_this_codec_block = 0;

    if (parameters.midside && (Global.Channels == 2))
        local_channels = 4;
    else
        local_channels = Global.Channels;

    for (this_channel = 0; this_channel < local_channels; this_channel++)
    {
        Global.Channel = this_channel;
        process.Channel_Data[Global.Channel].maximum_bits_to_remove = settings.static_maximum_bits_to_remove;
        process.Channel_Data[Global.Channel].min_FFT_result.btr = settings.static_maximum_bits_to_remove;
        process.Channel_Data[Global.Channel].min_FFT_result.analysis = 7;

        process.dynamic_maximum_bits_to_remove = int(std::max(0.0, floor(AudioData.Channel_Log2_RMS[Global.Channel] - settings.dynamic_minimum_bits_to_keep)));

        if (process.dynamic_maximum_bits_to_remove < process.Channel_Data[Global.Channel].maximum_bits_to_remove)
        {
            process.Channel_Data[Global.Channel].maximum_bits_to_remove = process.dynamic_maximum_bits_to_remove;
            process.Channel_Data[Global.Channel].min_FFT_result.btr = process.dynamic_maximum_bits_to_remove;
            process.Channel_Data[Global.Channel].min_FFT_result.analysis = 8;
        }

        for (this_analysis_number = 1; this_analysis_number <= PRECALC_ANALYSES; ++this_analysis_number)
        {
            Global.analysis.number = this_analysis_number;
            Global.FFT = FFT_PreCalc_Data_Rec[Global.analysis.bits[Global.analysis.number]];

            FFT_Proc_Rec this_FFT_plan;
            this_FFT_plan.NumberOfBitsNeeded = Global.analysis.bits[Global.analysis.number];
            this_FFT_plan.FFT = &FFT_PreCalc_Data_Rec[this_FFT_plan.NumberOfBitsNeeded];
            this_FFT_plan.FFT_Array = &FFT_Array;

            results.this_FFT_result.btr = 99;
            results.this_FFT_result.spreading = Max_dB;
            results.this_FFT_result.analysis = Global.analysis.number;

            if (((settings.FFT_analysis_switches & Global.FFT.length) > 0) && (Global.FFT.length <= (Global.Codec_Block.Size << 1)))
            {
                Zero_FFT_unity_results();
                this_FFT_plan.Task.analyses_performed = 0;

                for (this_analysis_block_number = 0; this_analysis_block_number <= process.analysis_blocks[Global.FFT.bit_length]; ++this_analysis_block_number)
                {
                    int32_t this_block_start = floor(process.actual_analysis_blocks_start[this_FFT_plan.FFT->bit_length] + this_analysis_block_number * process.FFT_underlap_length[this_FFT_plan.FFT->bit_length]);
                    this_FFT_plan.Task.block_start = std::min(std::max(this_block_start, process.limits.minstart),process.limits.maxend-this_FFT_plan.FFT->length);

                    if ((this_analysis_block_number == 0) && (this_FFT_plan.Task.block_start == results.saved_FFT_results[Global.Channel][Global.analysis.number].start))
                    {
                        if ((Global.FFT.length & settings.raw_result_short) > 0)
                            Fill_Last_from_History();

                        if ((Global.FFT.length & settings.raw_result_sgns) > 0)
                            Fill_Last_from_SGNS();

                        spreading.old_minimum = process.FFT_spreading[Global.analysis.number][Global.Channel].old_minimum;
                        spreading.new_minimum = process.FFT_spreading[Global.analysis.number][Global.Channel].new_minimum;
                        spreading.alt_average = process.FFT_spreading[Global.analysis.number][Global.Channel].alt_average;
                    }
                    else
                        if (FillFFT_Input_From_WAVE(&this_FFT_plan) == 0)
                        {
                            if ((Global.FFT.length & (settings.raw_result_short | settings.raw_result_sgns)) > 0)
                                Fill_Last_with_Zero();

                            spreading.old_minimum = Max_dB;
                            spreading.new_minimum = Max_dB;
                            spreading.alt_average = Max_dB;
                        }
                        else
                        {
                            if (FFTW_Initialised())
                                FFTW.Execute_R2C_New_Array(FFTW.Plans[this_FFT_plan.FFT->bit_length],&this_FFT_plan.DReal[0],&this_FFT_plan.DReal[0]);
                            else
                                FFT_DIT_Real(&this_FFT_plan);

                            process.Post_Analysis[Global.analysis.number](&this_FFT_plan);

                            Spreading_Function();
                        }

                    this_FFT_plan.Task.analyses_performed++;
                    Spreading_result = spreading.alt_average;

                    Spreading_Used = 1;

                    if (spreading.new_minimum < Spreading_result)
                    {
                        Spreading_result = spreading.new_minimum;
                        Spreading_Used = 2;
                    }

                    if (spreading.old_minimum < Spreading_result)
                    {
                        Spreading_result = spreading.old_minimum;
                        Spreading_Used = 3;
                    }

                    Spreading_result -= HannWindowRMS;
                    Spreading_result -= Spreading_result * (Spreading_result < 0);
                    results.this_FFT_result.spreading = Spreading_result;
                    int32_t this_spreading_index = int(std::min((THRESHOLD_INDEX_SPREAD_RANGE - 1.0d), Spreading_result) * THRESHOLD_INDEX_SPREAD);
                    results.this_FFT_result.btr = spreading.threshold_index[this_spreading_index];
                    results.this_FFT_result.start = this_FFT_plan.Task.block_start - Global.Codec_Block.Size;
                    results.this_FFT_result.analysis = Global.analysis.number;

                    if (parameters.output.spread != -1)
                    {
                        switch (Spreading_Used)
                        {
                            case 1:
                                ++ process.Alt_Ave_Used[Global.analysis.number];
                                break;

                            case 2:
                                ++ process.New_Min_Used_History[Global.analysis.number][process.new_min_bin];
                                ++ process.New_Min_Used[Global.analysis.number];
                                break;

                            case 3:
                                ++ process.Old_Min_Used_History[Global.analysis.number][process.old_min_bin];
                                ++ process.Old_Min_Used[Global.analysis.number];
                                break;

                            default:
                                lossyWAVError("Invalid spreading result",0x99);
                        }

                        if (results.this_FFT_result.btr > settings.static_maximum_bits_to_remove)
                            ++ process.Over_Static[Global.analysis.number];

                        if (results.this_FFT_result.btr > process.dynamic_maximum_bits_to_remove)
                            ++ process.Over_Dynamic[Global.analysis.number];
                    }

                    if ((results.this_FFT_result.btr < process.Channel_Data[Global.Channel].min_FFT_result.btr) || (process.Channel_Data[Global.Channel].min_FFT_result.btr == -1))
                        process.Channel_Data[Global.Channel].min_FFT_result = results.this_FFT_result;
                }

                process.Analyses_Completed[Global.analysis.number] += this_FFT_plan.Task.analyses_performed;

                if ((Global.FFT.length & settings.raw_result_short) > 0)
                    Fill_History_from_Last(&this_FFT_plan);

                if ((Global.FFT.length & settings.raw_result_sgns) > 0)
                    Fill_SGNSFFT_from_Last(&this_FFT_plan);

                results.saved_FFT_results[Global.Channel][Global.analysis.number] = results.this_FFT_result;
                process.FFT_spreading[Global.analysis.number][Global.Channel].old_minimum = spreading.old_minimum;
                process.FFT_spreading[Global.analysis.number][Global.Channel].new_minimum = spreading.new_minimum;
                process.FFT_spreading[Global.analysis.number][Global.Channel].alt_average = spreading.alt_average;
            }
        }


        if ((parameters.shaping.active) && (!parameters.shaping.fixed))
        {
            Make_Filter(Global.Channel);
        }


        if (process.Channel_Data[Global.Channel].min_FFT_result.btr < 0)
        {
            process.Channel_Data[Global.Channel].min_FFT_result.btr = 0;
        }


        process.Channel_Data[Global.Channel].calc_bits_to_remove = process.Channel_Data[Global.Channel].min_FFT_result.btr;
        process.Channel_Data[Global.Channel].bits_to_remove = process.Channel_Data[Global.Channel].calc_bits_to_remove;

        if (Global.Channel < Global.Channels)
        {
            Remove_Bits();

            codec_block_dependent_bits_to_remove = std::min(codec_block_dependent_bits_to_remove, process.Channel_Data[Global.Channel].calc_bits_to_remove);
        }
    }

    for (this_channel = 0; this_channel < Global.Channels; ++this_channel)
    {
        Global.Channel = this_channel;
        ++ results.minima[Global.Channel][process.Channel_Data[Global.Channel].min_FFT_result.analysis];

        if (parameters.linkchannels || parameters.midside)
        {
            if (process.Channel_Data[Global.Channel].bits_removed != codec_block_dependent_bits_to_remove)
            {
                if (Global.Channel < Global.Channels)
                {
                    process.Channel_Data[Global.Channel].bits_to_remove = codec_block_dependent_bits_to_remove;
                    Remove_Bits();
                }
            }
        }

        Stats.Incidence.eclip += process.Channel_Data[Global.Channel].Total.eclip;
        Stats.Incidence.sclip += process.Channel_Data[Global.Channel].Total.sclip;
        Stats.Incidence.rclip += process.Channel_Data[Global.Channel].Total.rclip;
        Stats.Incidence.aclip += process.Channel_Data[Global.Channel].Total.aclip;
        Stats.Incidence.round += process.Channel_Data[Global.Channel].Total.round;
        Stats.Incidence.noise += process.Channel_Data[Global.Channel].Total.noise;

        Stats.Count.eclips += process.Channel_Data[Global.Channel].Count.eclips;
        Stats.Count.sclips += process.Channel_Data[Global.Channel].Count.sclips;
        Stats.Count.rclips += process.Channel_Data[Global.Channel].Count.rclips;
        Stats.Count.aclips += process.Channel_Data[Global.Channel].Count.aclips;

        Stats.total_bits_removed += process.Channel_Data[Global.Channel].bits_removed;
        Stats.total_bits_lost += process.Channel_Data[Global.Channel].bits_lost;

        ++ Stats.bits_removed[Global.Channel][process.Channel_Data[Global.Channel].bits_removed];
        ++ Stats.bits_lost[Global.Channel][process.Channel_Data[Global.Channel].bits_lost];

        if (parameters.output.detail)
        {
            bit_removal_history[(Global.blocks_processed - 1) * Global.Channels + Global.Channel] = process.Channel_Data[Global.Channel].bits_removed;
        }

        bits_removed_this_codec_block += OneOver[Global.Channels] * process.Channel_Data[Global.Channel].bits_removed;
    }

    //==========================================================================
    // Post analyse bit removed / correction audio data using uint16_t FFT.
    //==========================================================================
    if (parameters.output.postanalyse)
    {
        Global.FFT = FFT_PreCalc_Data_Rec[history.FFT.bit_length];

        FFT_Proc_Rec this_FFT_plan;
        this_FFT_plan.NumberOfBitsNeeded = history.FFT.bit_length;
        this_FFT_plan.FFT_Array = &FFT_Array;
        this_FFT_plan.FFT = &FFT_PreCalc_Data_Rec[history.FFT.bit_length];

        for (this_channel = 0; this_channel < Global.Channels; ++this_channel)
        {
            Global.Channel = this_channel;

            //======================================================================
            // Create spectrum for BTRD-Codec-Block
            //======================================================================
            this_FFT_plan.Task.Fill_FFT_Proc = FillFFT_Input_From_BTRD;
            this_FFT_plan.Task.Dist_Array = &history.BTRD_results[Global.Channel][0];
            Process_Data_to_LongDist(&this_FFT_plan);

            //======================================================================
            // Create spectrum for CORR-Codec-Block
            //======================================================================
            this_FFT_plan.Task.Fill_FFT_Proc = FillFFT_Input_From_CORR;
            this_FFT_plan.Task.Dist_Array = &history.CORR_results[Global.Channel][0];
            Process_Data_to_LongDist(&this_FFT_plan);
        }
    }

    //==========================================================================
    // Post analyse DATA/BTRD/CORR audio data using Codec_Block_Length FFT.
    //==========================================================================
    if (parameters.output.longdist)
    {
        Global.FFT = FFT_PreCalc_Data_Rec[LongDist.FFT.bit_length];

        FFT_Proc_Rec this_FFT_plan;
        this_FFT_plan.NumberOfBitsNeeded = LongDist.FFT.bit_length;
        this_FFT_plan.FFT_Array = &FFT_Array;
        this_FFT_plan.FFT = &FFT_PreCalc_Data_Rec[LongDist.FFT.bit_length];

        for (this_channel = 0; this_channel < Global.Channels; ++this_channel)
        {
            Global.Channel = this_channel;

            //======================================================================
            // Create spectrum for DATA-Codec-Block
            //======================================================================
            this_FFT_plan.Task.Fill_FFT_Proc = FillFFT_Input_From_WAVE;
            this_FFT_plan.Task.Dist_Array = &LongDist.WAVE_results[Global.Channel][0];
            Process_Data_to_LongDist(&this_FFT_plan);

            //======================================================================
            // Create spectrum for BTRD-Codec-Block
            //======================================================================
            this_FFT_plan.Task.Fill_FFT_Proc = FillFFT_Input_From_BTRD;
            this_FFT_plan.Task.Dist_Array = &LongDist.BTRD_results[Global.Channel][0];
            Process_Data_to_LongDist(&this_FFT_plan);

            //======================================================================
            // Create spectrum for CORR-Codec-Block
            //======================================================================
            this_FFT_plan.Task.Fill_FFT_Proc = FillFFT_Input_From_CORR;
            this_FFT_plan.Task.Dist_Array = &LongDist.CORR_results[Global.Channel][0];
            Process_Data_to_LongDist(&this_FFT_plan);
        }
    }

    //==========================================================================
    // Create DATA, BTRD and CORR sample value histogram
    //==========================================================================
    if (parameters.output.histogram)
        Add_to_Histogram();

    Process_Output();

    //==========================================================================
    // Create histogram of LSBs and MSBs on a sample / block basis.
    //==========================================================================
    if (parameters.output.blockdist || parameters.output.sampledist)
        lsb_analysis();
}


void nProcess_Init()
{
    int32_t total_overlap_length;

    for (int32_t sa_j = 1; sa_j <= PRECALC_ANALYSES; ++sa_j)
    {
        for (int32_t sa_i = 0; sa_i < int(MAX_FFT_LENGTH_HALF); ++sa_i)
        {
            process.Old_Min_Used_History[sa_j][sa_i] = 0;
            process.New_Min_Used_History[sa_j][sa_i] = 0;
        }

        process.Old_Min_Used[sa_j] = 0;
        process.New_Min_Used[sa_j] = 0;
        process.Alt_Ave_Used[sa_j] = 0;
        process.Over_Dynamic[sa_j] = 0;
        process.Over_Static[sa_j] = 0;
        process.Post_Analysis[sa_j] = Post_Process_FFT_Results_No_Raw;
    }

    Global.feedback.eclip = 0;
    Global.feedback.sclip = 0;
    Global.feedback.aclip = 0;
    Global.feedback.noise = 0;
    Global.feedback.rclip = 0;
    Global.feedback.round = 0;

    if (settings.raw_result_short != 0)
        process.Post_Analysis[3] = Post_Process_FFT_Results;

    if (settings.raw_result_sgns != 0)
        process.Post_Analysis[7] = Post_Process_FFT_Results;

    for (int32_t sa_i = 1; sa_i <= MAX_FFT_BIT_LENGTH; ++sa_i)
    {
        Global.FFT = FFT_PreCalc_Data_Rec[sa_i];

        if (FFTW_Initialised())
            if ((((settings.FFT_analysis_switches & Global.FFT.length) > 0) || (parameters.output.longdist && (Global.FFT.bit_length == LongDist.FFT.bit_length))) || parameters.noisecalc)
            {
                FFTW.Plans[Global.FFT.bit_length] = FFTW.Plan_DFT_r2c_1d(Global.FFT.length, &FFT_Array.DReal[0], &FFT_Array.DReal[0], 0x69);

                if (parameters.noisecalc)
                {
                    FFTW.Plans_Inv[Global.FFT.bit_length] = FFTW.Plan_DFT_1d(Global.FFT.length, &FFT_Array.DReal[0], &FFT_Array.DReal[0], FFTW_BACKWARD, 0x69);
                }
            }
        //========================================================================
        // Calculate actual analysis_time values for FFT_lengths
        //========================================================================
        process.end_overlap_length[Global.FFT.bit_length] = std::min(Global.FFT.length / 2, Global.Codec_Block.Size);
        process.actual_analysis_blocks_start[Global.FFT.bit_length] = -process.end_overlap_length[Global.FFT.bit_length];

        total_overlap_length = Global.Codec_Block.Size + (process.end_overlap_length[Global.FFT.bit_length] << 1) - Global.FFT.length;
        process.FFT_underlap_length[Global.FFT.bit_length] = std::max(1, Global.FFT.length / std::max(2, parameters.fft.underlap));

        process.analysis_blocks[Global.FFT.bit_length] = std::max(0, nRoundEvenInt32(total_overlap_length / process.FFT_underlap_length[Global.FFT.bit_length]));

        process.FFT_underlap_length[Global.FFT.bit_length] = total_overlap_length / std::max(1, process.analysis_blocks[Global.FFT.bit_length]);
    }
}
