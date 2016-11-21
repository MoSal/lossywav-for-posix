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
#include "nMaths.h"
#include "nSpreading.h"
#include "fftw_interface.h"
#include "nFFT.h"
#include "nFillFFT.h"
#include "nRemoveBits.h"
#include "nOutput.h"
#include "nProcess.h"


void Add_to_Unity(Results_Type* this_result)
{
    for (int32_t sa_i = 0; sa_i <= Current.Analysis.upper_process_bin[Current.Analysis.number]; ++sa_i)
    {
        this_result->Unity[sa_i] += this_result->LastUnity[sa_i];
        this_result->Root[sa_i] += this_result->LastRoot[sa_i];
    }
}


void Add_to_History(FFT_Proc_Rec* this_FFT_plan, Results_Type* this_result)
{
    double divisor_a = OneOver[this_FFT_plan->Task.analyses_performed - 1];

    for (int32_t sa_i = 0; sa_i <= Current.Analysis.upper_process_bin[Current.Analysis.number]; sa_i++)
    {
        double unity_diff = (this_result->Unity[sa_i] - this_result->LastUnity[sa_i]);

        double temp_result = divisor_a * unity_diff;

        this_result->History[sa_i] += temp_result;

        this_result->SGNSUnity[sa_i] = temp_result;

        double root_diff = (this_result->Root[sa_i] - this_result->LastRoot[sa_i]);

        this_result->SGNSRoot[sa_i] = divisor_a * root_diff;

        this_result->SGNSHybrid[sa_i] = (unity_diff * root_diff);
    }
}


void Fill_Last_with_Zero(Results_Type* this_result)
{
    for (int32_t sa_i = 0; sa_i <= Current.Analysis.upper_process_bin[Current.Analysis.number]; ++sa_i)
    {
        this_result->LastUnity[sa_i] = 0;
        this_result->LastRoot[sa_i] = 0;
        this_result->Skewed[sa_i] = 0;
    }
}


void Zero_FFT_unity_results(Results_Type* this_result)
{
    for (int32_t sa_i = 0; sa_i <= Current.Analysis.upper_process_bin[Current.Analysis.number]; ++sa_i)
    {
        this_result->Unity[sa_i] = 0;
        this_result->Root[sa_i] = 0;
    }
}


void Post_Process_FFT_Results(FFT_Proc_Rec* this_FFT_plan, Results_Type* this_result)
{
    for (int32_t sa_i = 0; sa_i <= Current.Analysis.upper_process_bin[Current.Analysis.number]; ++sa_i)
    {
        double sc_y = this_FFT_plan->DComplex[sa_i].magnitude();

        this_result->LastRoot[sa_i] = sc_y;
        this_result->Root[sa_i] += sc_y;

        this_result->Skewed[sa_i] = sc_y * Skewing_Gain[sa_i << Current.FFT.bit_shift_from_max];

        sc_y *= sc_y;
        this_result->LastUnity[sa_i] = sc_y;
        this_result->Unity[sa_i] += sc_y;
    }
}


void Process_Data_to_History(FFT_Proc_Rec* this_FFT_plan, Results_Type* this_result)
{
    Zero_FFT_unity_results(this_result);

    this_FFT_plan->Task.analyses_performed = 0;

    for (int32_t this_analysis_block_number = 0; this_analysis_block_number < (process.analysis_blocks[Current.FFT.bit_length] + Global.last_codec_block); ++this_analysis_block_number)
    {
        int32_t this_block_start = floor(process.actual_analysis_blocks_start[Current.FFT.bit_length] + this_analysis_block_number * process.FFT_underlap_length[Current.FFT.bit_length]);
        this_FFT_plan->Task.block_start = std::min(std::max(this_block_start, process.limits.minstart),process.limits.maxend-Current.FFT.length);
        this_FFT_plan->Task.analyses_performed++;

        if (this_FFT_plan->Task.Fill_FFT_Proc(this_FFT_plan) == 0)
        {
            Fill_Last_with_Zero(this_result);
        }
        else
        {
            if (FFTW_Initialised())
                FFTW.Execute_R2C_New_Array(FFTW.Plans[this_FFT_plan->FFT->bit_length],&this_FFT_plan->DReal[0],&this_FFT_plan->DReal[0]);
            else
                FFT_DIT_Real(this_FFT_plan);

            for (int32_t sa_i = 0; sa_i <= this_FFT_plan->FFT->length_half; ++sa_i)
            {
                double sc_y = this_FFT_plan->DComplex[sa_i].magsquared();
                this_result->LastUnity[sa_i] = sc_y;
                this_result->Unity[sa_i] += sc_y;
            }
        }
    }

    if (Global.last_codec_block)
    {
        Add_to_Unity(this_result);
        this_FFT_plan->Task.analyses_performed++;
    }

    double divisor = OneOver[this_FFT_plan->Task.analyses_performed];

    for (int32_t sa_i = 0; sa_i <= this_FFT_plan->FFT->length_half; sa_i++)
    {
        this_result->History[sa_i] += divisor * this_result->Unity[sa_i];
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
        Current.Channel = this_channel;
        process.Channel_Data[Current.Channel].maximum_bits_to_remove = settings.static_maximum_bits_to_remove;
        process.Channel_Data[Current.Channel].min_FFT_result.btr = settings.static_maximum_bits_to_remove;
        process.Channel_Data[Current.Channel].min_FFT_result.analysis = 7;

        process.dynamic_maximum_bits_to_remove = int(std::max(0.0, floor(AudioData.Channel_Log2_RMS[Current.Channel] - settings.dynamic_minimum_bits_to_keep)));

        if (process.dynamic_maximum_bits_to_remove < process.Channel_Data[Current.Channel].maximum_bits_to_remove)
        {
            process.Channel_Data[Current.Channel].maximum_bits_to_remove = process.dynamic_maximum_bits_to_remove;
            process.Channel_Data[Current.Channel].min_FFT_result.btr = process.dynamic_maximum_bits_to_remove;
            process.Channel_Data[Current.Channel].min_FFT_result.analysis = 8;
        }

        for (this_analysis_number = 1; this_analysis_number <= PRECALC_ANALYSES; ++this_analysis_number)
        {
            Current.Analysis.number = this_analysis_number;
            Current.FFT = FFT_PreCalc_Data_Rec[Current.Analysis.bits[Current.Analysis.number]];

            FFT_Proc_Rec this_FFT_plan;
            this_FFT_plan.NumberOfBitsNeeded = Current.FFT.bit_length;
            this_FFT_plan.FFT = &FFT_PreCalc_Data_Rec[Current.FFT.bit_length];
            this_FFT_plan.FFT_Array = &FFT_Array;

            results.this_FFT_result.btr = 99;
            results.this_FFT_result.spreading = Max_dB;
            results.this_FFT_result.analysis = Current.Analysis.number;

            Results_Type* this_result = &results.WAVE[Current.Analysis.number][Current.Channel];

            if (settings.analysis[this_analysis_number].active)
            {
                Zero_FFT_unity_results(this_result);

                this_FFT_plan.Task.analyses_performed = 0;

                for (this_analysis_block_number = 0; this_analysis_block_number <= process.analysis_blocks[Current.FFT.bit_length]; ++this_analysis_block_number)
                {
                    int32_t this_block_start = floor(process.actual_analysis_blocks_start[this_FFT_plan.FFT->bit_length] + this_analysis_block_number * process.FFT_underlap_length[this_FFT_plan.FFT->bit_length]);
                    this_FFT_plan.Task.block_start = std::min(std::max(this_block_start, process.limits.minstart),process.limits.maxend-this_FFT_plan.FFT->length);
                    this_FFT_plan.Task.analyses_performed++;

                    if ((this_analysis_block_number == 0) && (!Global.first_codec_block))
                    {
                        Add_to_Unity(this_result);

                        spreading.old_minimum = process.FFT_spreading[Current.Analysis.number][Current.Channel].old_minimum;
                        spreading.new_minimum = process.FFT_spreading[Current.Analysis.number][Current.Channel].new_minimum;
                        spreading.alt_average = process.FFT_spreading[Current.Analysis.number][Current.Channel].alt_average;
                    }
                    else
                        if (FillFFT_Input_From_WAVE(&this_FFT_plan) == 0)
                        {
                            Fill_Last_with_Zero(this_result);

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

                            Post_Process_FFT_Results(&this_FFT_plan, this_result);

                            Spreading_Function(this_result);
                        }

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
                    int32_t this_spreading_index = int(std::min((THRESHOLD_INDEX_SPREAD_RANGE - 1.0), Spreading_result) * THRESHOLD_INDEX_SPREAD);
                    results.this_FFT_result.btr = spreading.threshold_index[this_spreading_index];
                    results.this_FFT_result.start = this_FFT_plan.Task.block_start - Global.Codec_Block.Size;
                    results.this_FFT_result.analysis = Current.Analysis.number;

                    if ((parameters.output.spread != -1))
                    {
                        switch (Spreading_Used)
                        {
                            case 1:
                                ++ process.Alt_Ave_Used[Current.Analysis.number];
                                break;

                            case 2:
                                ++ process.New_Min_Used_History[Current.Analysis.number][process.new_min_bin];
                                ++ process.New_Min_Used[Current.Analysis.number];
                                break;

                            case 3:
                                ++ process.Old_Min_Used_History[Current.Analysis.number][process.old_min_bin];
                                ++ process.Old_Min_Used[Current.Analysis.number];
                                break;

                            default:
                                lossyWAVError("Invalid spreading result",0x99);
                        }

                        if (results.this_FFT_result.btr > settings.static_maximum_bits_to_remove)
                            ++ process.Over_Static[Current.Analysis.number];

                        if (results.this_FFT_result.btr > process.dynamic_maximum_bits_to_remove)
                            ++ process.Over_Dynamic[Current.Analysis.number];
                    }

                    if ((results.this_FFT_result.btr < process.Channel_Data[Current.Channel].min_FFT_result.btr) || (process.Channel_Data[Current.Channel].min_FFT_result.btr == -1))
                        process.Channel_Data[Current.Channel].min_FFT_result = results.this_FFT_result;
                }

                if (Global.last_codec_block)
                {
                    Add_to_Unity(this_result);
                    this_FFT_plan.Task.analyses_performed++;
                }

                process.Analyses_Completed[Current.Analysis.number] += (this_FFT_plan.Task.analyses_performed);

                Add_to_History(&this_FFT_plan, this_result);

                results.saved_FFT_results[Current.Channel][Current.Analysis.number] = results.this_FFT_result;
                process.FFT_spreading[Current.Analysis.number][Current.Channel].old_minimum = spreading.old_minimum;
                process.FFT_spreading[Current.Analysis.number][Current.Channel].new_minimum = spreading.new_minimum;
                process.FFT_spreading[Current.Analysis.number][Current.Channel].alt_average = spreading.alt_average;
            }
        }


        if ((parameters.shaping.active) && (!parameters.shaping.fixed))
        {
            Make_Filter(Current.Channel);
        }


        if (process.Channel_Data[Current.Channel].min_FFT_result.btr < 0)
        {
            process.Channel_Data[Current.Channel].min_FFT_result.btr = 0;
        }


        process.Channel_Data[Current.Channel].calc_bits_to_remove = process.Channel_Data[Current.Channel].min_FFT_result.btr;
        process.Channel_Data[Current.Channel].bits_to_remove = process.Channel_Data[Current.Channel].calc_bits_to_remove;

        if (Current.Channel < Global.Channels)
        {
            Remove_Bits();

            codec_block_dependent_bits_to_remove = std::min(codec_block_dependent_bits_to_remove, process.Channel_Data[Current.Channel].calc_bits_to_remove);
        }
    }

    for (this_channel = 0; this_channel < Global.Channels; ++this_channel)
    {
        Current.Channel = this_channel;
        ++ results.minima[Current.Channel][process.Channel_Data[Current.Channel].min_FFT_result.analysis];

        if (parameters.linkchannels || parameters.midside)
        {
            if (process.Channel_Data[Current.Channel].bits_removed != codec_block_dependent_bits_to_remove)
            {
                if (Current.Channel < Global.Channels)
                {
                    process.Channel_Data[Current.Channel].bits_to_remove = codec_block_dependent_bits_to_remove;
                    Remove_Bits();
                }
            }
        }

        Stats.Incidence.eclip += process.Channel_Data[Current.Channel].Total.eclip;
        Stats.Incidence.sclip += process.Channel_Data[Current.Channel].Total.sclip;
        Stats.Incidence.rclip += process.Channel_Data[Current.Channel].Total.rclip;
        Stats.Incidence.aclip += process.Channel_Data[Current.Channel].Total.aclip;
        Stats.Incidence.round += process.Channel_Data[Current.Channel].Total.round;
        Stats.Incidence.noise += process.Channel_Data[Current.Channel].Total.noise;

        Stats.Count.eclips += process.Channel_Data[Current.Channel].Count.eclips;
        Stats.Count.sclips += process.Channel_Data[Current.Channel].Count.sclips;
        Stats.Count.rclips += process.Channel_Data[Current.Channel].Count.rclips;
        Stats.Count.aclips += process.Channel_Data[Current.Channel].Count.aclips;

        Stats.total_bits_removed += process.Channel_Data[Current.Channel].bits_removed;
        Stats.total_bits_lost += process.Channel_Data[Current.Channel].bits_lost;

        ++ Stats.bits_removed[Current.Channel][process.Channel_Data[Current.Channel].bits_removed];
        ++ Stats.bits_lost[Current.Channel][process.Channel_Data[Current.Channel].bits_lost];

        if (parameters.output.detail)
        {
            bit_removal_history[(Global.blocks_processed - 1) * Global.Channels + Current.Channel] = process.Channel_Data[Current.Channel].bits_removed;
        }

        bits_removed_this_codec_block += OneOver[Global.Channels] * process.Channel_Data[Current.Channel].bits_removed;
    }

    //==========================================================================
    // Post analyse bit removed / correction audio data using uint16_t FFT.
    //==========================================================================
    if (parameters.output.postanalyse)
    {
        for (int32_t this_analysis_number = 1; this_analysis_number <= PRECALC_ANALYSES; ++this_analysis_number)
        {
            Current.Analysis.number = this_analysis_number;

            Analysis_Type* this_analysis = &settings.analysis[this_analysis_number];

            if ((this_analysis->active) && ((parameters.output.longdist) || (this_analysis_number == SHORT_ANALYSIS)))
            {
                Current.FFT = FFT_PreCalc_Data_Rec[this_analysis->FFT.bit_length];

                FFT_Proc_Rec this_FFT_plan;
                this_FFT_plan.NumberOfBitsNeeded = Current.FFT.bit_length;
                this_FFT_plan.FFT_Array = &FFT_Array;
                this_FFT_plan.FFT = &FFT_PreCalc_Data_Rec[Current.FFT.bit_length];

                for (this_channel = 0; this_channel < Global.Channels; ++this_channel)
                {
                    Current.Channel = this_channel;

                    //======================================================================
                    // Create spectrum for BTRD-Codec-Block
                    //======================================================================
                    this_FFT_plan.Task.Fill_FFT_Proc = FillFFT_Input_From_BTRD;
                    Process_Data_to_History(&this_FFT_plan, &results.BTRD[Current.Analysis.number][Current.Channel]);

                    //======================================================================
                    // Create spectrum for CORR-Codec-Block
                    //======================================================================
                    this_FFT_plan.Task.Fill_FFT_Proc = FillFFT_Input_From_CORR;
                    Process_Data_to_History(&this_FFT_plan, &results.CORR[Current.Analysis.number][Current.Channel]);
                }
            }
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


void nProcess_Initialise_Results_Arrays(Results_Type* this_result, int32_t this_analysis)
{
    this_result->History = new double[sizeof(double) * Current.Analysis.length[this_analysis]];
    this_result->Unity = new double[sizeof(double) * Current.Analysis.length[this_analysis]];
    this_result->LastUnity = new double[sizeof(double) * Current.Analysis.length[this_analysis]];
    this_result->SGNSUnity = new double[sizeof(double) * Current.Analysis.length[this_analysis]];
    this_result->Root = new double[sizeof(double) * Current.Analysis.length[this_analysis]];
    this_result->LastRoot = new double[sizeof(double) * Current.Analysis.length[this_analysis]];
    this_result->SGNSRoot = new double[sizeof(double) * Current.Analysis.length[this_analysis]];
    this_result->SGNSHybrid = new double[sizeof(double) * Current.Analysis.length[this_analysis]];
    this_result->Skewed = new double[sizeof(double) * Current.Analysis.length[this_analysis]];

    for (int32_t sa_k = 0; sa_k < Current.Analysis.length[this_analysis]; sa_k++)
    {
        this_result->History[sa_k] = 0.;
        this_result->Unity[sa_k] = 0.;
        this_result->LastUnity[sa_k] = 0.;
        this_result->SGNSUnity[sa_k] = 0.;
        this_result->Root[sa_k] = 0.;
        this_result->LastRoot[sa_k] = 0.;
        this_result->SGNSRoot[sa_k] = 0.;
        this_result->SGNSHybrid[sa_k] = 0.;
        this_result->Skewed[sa_k] = 0.;
    }
    this_result->Analyses_Performed = 0;
}

void nProcess_Init()
{
    int32_t total_overlap_length;

    for (int32_t this_analysis = 1; this_analysis <= PRECALC_ANALYSES; ++this_analysis)
    {
        for (int32_t sa_i = 0; sa_i < int(MAX_FFT_LENGTH_HALF); ++sa_i)
        {
            process.Old_Min_Used_History[this_analysis][sa_i] = 0;
            process.New_Min_Used_History[this_analysis][sa_i] = 0;
        }

        process.Old_Min_Used[this_analysis] = 0;
        process.New_Min_Used[this_analysis] = 0;
        process.Alt_Ave_Used[this_analysis] = 0;
        process.Over_Dynamic[this_analysis] = 0;
        process.Over_Static[this_analysis] = 0;

        if ((settings.analysis[this_analysis].active) || (true))
        {
            Current.FFT = FFT_PreCalc_Data_Rec[settings.analysis[this_analysis].FFT.bit_length];

            if (FFTW_Initialised())
            {
                FFTW.Plans[Current.FFT.bit_length] = FFTW.Plan_DFT_r2c_1d(Current.FFT.length, &FFT_Array.DReal[0], &FFT_Array.DReal[0], 0x69);
            }

            int32_t this_upper_bin = spreading.Bins.Upper[settings.analysis[this_analysis].FFT.bit_length] + 1;

            if ((this_analysis == SHORT_ANALYSIS) || (this_analysis == LONG_ANALYSIS) || (parameters.shaping.hybrid) || (parameters.output.freqdist))
            {
                this_upper_bin = settings.analysis[this_analysis].FFT.length_half;
            }

            Current.Analysis.upper_process_bin[this_analysis] = this_upper_bin;

            for (int32_t sa_i = 0; sa_i < (MAX_CHANNELS); sa_i++)
            {
                nProcess_Initialise_Results_Arrays(&results.WAVE[this_analysis][sa_i], this_analysis);
                nProcess_Initialise_Results_Arrays(&results.BTRD[this_analysis][sa_i], this_analysis);
                nProcess_Initialise_Results_Arrays(&results.CORR[this_analysis][sa_i], this_analysis);
            }
        }
    }

    Global.feedback.eclip = 0;
    Global.feedback.sclip = 0;
    Global.feedback.aclip = 0;
    Global.feedback.noise = 0;
    Global.feedback.rclip = 0;
    Global.feedback.round = 0;

    for (int32_t sa_i = 1; sa_i <= MAX_FFT_BIT_LENGTH; ++sa_i)
    {
        Current.FFT = FFT_PreCalc_Data_Rec[sa_i];

        //========================================================================
        // Calculate actual analysis_time values for FFT_lengths
        //========================================================================
        process.end_overlap_length[Current.FFT.bit_length] = std::min(Current.FFT.length / 2, Global.Codec_Block.Size);
        process.actual_analysis_blocks_start[Current.FFT.bit_length] = -process.end_overlap_length[Current.FFT.bit_length];

        total_overlap_length = Global.Codec_Block.Size + (process.end_overlap_length[Current.FFT.bit_length] << 1) - Current.FFT.length;
        process.FFT_underlap_length[Current.FFT.bit_length] = std::max(1, Current.FFT.length / std::max(2, parameters.fft.underlap));

        process.analysis_blocks[Current.FFT.bit_length] = std::max(0, nRoundEvenInt32(total_overlap_length / process.FFT_underlap_length[Current.FFT.bit_length]));

        process.FFT_underlap_length[Current.FFT.bit_length] = total_overlap_length / std::max(1, process.analysis_blocks[Current.FFT.bit_length]);
    }
}

void nProcess_Cleanup_Results_Arrays(Results_Type* this_result)
{
    if (this_result->History != nullptr)
    {
        delete[] this_result->History;
    }

    if (this_result->Unity != nullptr)
    {
        delete[] this_result->Unity;
    }

    if (this_result->LastUnity != nullptr)
    {
        delete[] this_result->LastUnity;
    }

    if (this_result->SGNSUnity != nullptr)
    {
        delete[] this_result->SGNSUnity;
    }

    if (this_result->Root != nullptr)
    {
        delete[] this_result->Root;
    }

    if (this_result->LastRoot != nullptr)
    {
        delete[] this_result->LastRoot;
    }

    if (this_result->SGNSRoot != nullptr)
    {
        delete[] this_result->SGNSRoot;
    }

    if (this_result->SGNSHybrid != nullptr)
    {
        delete[] this_result->SGNSHybrid;
    }

    if (this_result->Skewed != nullptr)
    {
        delete[] this_result->Skewed;
    }

}

void nProcess_Cleanup()
{
    for (int32_t sa_i = 1; sa_i < (PRECALC_ANALYSES + 1); ++sa_i)
        for (int32_t sa_j = 0; sa_j < MAX_CHANNELS; sa_j++)
        {
            nProcess_Cleanup_Results_Arrays(&results.WAVE[sa_i][sa_j]);
            nProcess_Cleanup_Results_Arrays(&results.BTRD[sa_i][sa_j]);
            nProcess_Cleanup_Results_Arrays(&results.CORR[sa_i][sa_j]);
        }
}
