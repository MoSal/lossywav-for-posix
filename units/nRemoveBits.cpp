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

#include "nRemoveBits.h"

int32_t sample_maxima[BITS_TO_CALCULATE];

struct Removal_type
{
    double this_two_power;
    double this_two_power_recip;
    int64_t this_min_sample;
    int64_t this_max_sample;

} Removal __attribute__ ((aligned(16)));

int64_t Limited_Sample_Int64(int64_t lsv)
{
    return lsv + ((Removal.this_max_sample - lsv) & -(Removal.this_max_sample < lsv)) + ((Removal.this_min_sample - lsv) & -(Removal.this_min_sample > lsv));
}


//=============================================================================================================================
void Store_Proc()
{//============================================================================================================================
    process.Channel_Data[Global.Channel].Count.eclips = 0;
    process.Channel_Data[Global.Channel].Count.sclips = 0;
    process.Channel_Data[Global.Channel].Count.rclips = 0;
    process.Channel_Data[Global.Channel].Count.aclips = 0;

    Removal.this_max_sample = sample_maxima[process.Channel_Data[Global.Channel].bits_removed];

    for (int32_t count = 0; count < AudioData.Size.This; ++count)
    {
        int32_t this_DATA = AudioData.WAVEPTR[THIS_CODEC_BLOCK][Global.Channel][count].Integers[0];

        bool extant_clip = ((this_DATA > Removal.this_max_sample) | (this_DATA == Removal.this_min_sample));

        process.Channel_Data[Global.Channel].Count.eclips += extant_clip;

        double scaled =  this_DATA * settings.scaling_factor;

        bool scaled_clip = ((scaled > Removal.this_max_sample) | (scaled < Removal.this_min_sample)) & (!extant_clip);

        process.Channel_Data[Global.Channel].Count.sclips += scaled_clip;

        int64_t this_BTRD = nRoundEvenInt64(scaled);

        int64_t this_LIMIT = Limited_Sample_Int64(this_BTRD);

        process.Channel_Data[Global.Channel].Count.rclips += ((this_LIMIT != this_BTRD) & (!scaled_clip));

        AudioData.BTRDPTR[THIS_CODEC_BLOCK][Global.Channel][count].Integers[0] = this_LIMIT;
        AudioData.CORRPTR[THIS_CODEC_BLOCK][Global.Channel][count].Integers[0] = this_DATA - nRoundEvenInt64(this_LIMIT * settings.scaling_factor_inv);
    }

    double this_round = 0;

    process.Channel_Data[Global.Channel].Incidence.eclip = (process.Channel_Data[Global.Channel].Count.eclips > 0);
    process.Channel_Data[Global.Channel].Incidence.sclip = (process.Channel_Data[Global.Channel].Count.sclips > 0);

    process.Channel_Data[Global.Channel].Incidence.rclip = (process.Channel_Data[Global.Channel].Count.rclips > parameters.feedback.rclips);
    process.Channel_Data[Global.Channel].Incidence.retry = process.Channel_Data[Global.Channel].Incidence.rclip;

    process.Channel_Data[Global.Channel].Incidence.round = (parameters.feedback.active & (this_round > parameters.feedback.round) & (!process.Channel_Data[Global.Channel].Incidence.retry));
    process.Channel_Data[Global.Channel].Incidence.retry = process.Channel_Data[Global.Channel].Incidence.round;

    process.Channel_Data[Global.Channel].Incidence.noise = false;
    process.Channel_Data[Global.Channel].Incidence.aclip = false;
}//============================================================================================================================


//=============================================================================================================================
void Remove_Bits_Proc_Adaptive_Noise_Shaping_On()
{//============================================================================================================================
    double this_noise = 0;
    double this_round = 0;

    Warped_Lattice_Filter_Init(Global.Channel);

    process.Channel_Data[Global.Channel].Count.eclips = 0;
    process.Channel_Data[Global.Channel].Count.sclips = 0;
    process.Channel_Data[Global.Channel].Count.rclips = 0;
    process.Channel_Data[Global.Channel].Count.aclips = 0;

    Removal.this_max_sample = sample_maxima[process.Channel_Data[Global.Channel].bits_removed];
    Removal.this_two_power = PowersOf.Two[process.Channel_Data[Global.Channel].bits_removed];
    Removal.this_two_power_recip = PowersOf.Two[-process.Channel_Data[Global.Channel].bits_removed];

    double this_DATA_sqr = 0;
    double this_LIMIT_sqr = 0;
    double this_WLFE_sqr = 0;

    double limit_WLFE = npower2(AudioData.Channel_Log2_RMS[Global.Channel] + parameters.feedback.alevel);

    for (int32_t count = 0; count < AudioData.Size.This; ++count)
    {
        int32_t this_DATA = AudioData.WAVEPTR[THIS_CODEC_BLOCK][Global.Channel][count].Integers[0];

        bool extant_clip = (this_DATA > Removal.this_max_sample) | (this_DATA < Removal.this_min_sample);
        process.Channel_Data[Global.Channel].Count.eclips += extant_clip;

        double scaled =  this_DATA * settings.scaling_factor;

        bool scaled_clip = (scaled > Removal.this_max_sample) | (scaled < Removal.this_min_sample);
        process.Channel_Data[Global.Channel].Count.sclips += scaled_clip;

        double this_WLFE = Warped_Lattice_Filter_Evaluate(Global.Channel);
        this_WLFE_sqr += (this_WLFE * this_WLFE);
        process.Channel_Data[Global.Channel].Count.aclips += (fabs(this_WLFE) > limit_WLFE);

        int64_t this_SHAPED = nRoundEvenInt64((scaled + this_WLFE) * Removal.this_two_power_recip) * Removal.this_two_power;

        int64_t this_LIMIT = Limited_Sample_Int64(this_SHAPED);
        process.Channel_Data[Global.Channel].Count.rclips += ((this_LIMIT != this_SHAPED) & (!scaled_clip));

        Warped_Lattice_Filter_Update(Global.Channel, this_LIMIT - scaled);

        AudioData.BTRDPTR[THIS_CODEC_BLOCK][Global.Channel][count].Integers[0] = this_LIMIT;
        AudioData.CORRPTR[THIS_CODEC_BLOCK][Global.Channel][count].Integers[0] = this_DATA - nRoundEvenInt64(this_LIMIT * settings.scaling_factor_inv);

        this_DATA_sqr += (scaled * scaled);
        this_LIMIT_sqr += (this_LIMIT * this_LIMIT);
    }

    this_noise = (0.5d * nlog2(this_WLFE_sqr)) - Global.bits_per_sample;

    this_round = (nlog2(this_LIMIT_sqr) - nlog2(this_DATA_sqr));

    process.Channel_Data[Global.Channel].Incidence.eclip = (process.Channel_Data[Global.Channel].Count.eclips > 0);
    process.Channel_Data[Global.Channel].Incidence.sclip = (process.Channel_Data[Global.Channel].Count.sclips > 0);

    process.Channel_Data[Global.Channel].Incidence.rclip = (process.Channel_Data[Global.Channel].Count.rclips > parameters.feedback.rclips);
    process.Channel_Data[Global.Channel].Incidence.retry = process.Channel_Data[Global.Channel].Incidence.rclip;

    process.Channel_Data[Global.Channel].Incidence.round = (parameters.feedback.active & (this_round > parameters.feedback.round) & (!process.Channel_Data[Global.Channel].Incidence.retry));
    process.Channel_Data[Global.Channel].Incidence.retry |= process.Channel_Data[Global.Channel].Incidence.round;

    process.Channel_Data[Global.Channel].Incidence.noise = (parameters.feedback.active & (this_noise > parameters.feedback.noise) & (!process.Channel_Data[Global.Channel].Incidence.retry));
    process.Channel_Data[Global.Channel].Incidence.retry |= process.Channel_Data[Global.Channel].Incidence.noise;

    process.Channel_Data[Global.Channel].Incidence.aclip = (parameters.feedback.active & (process.Channel_Data[Global.Channel].Count.aclips > parameters.feedback.aclips) & (!(process.Channel_Data[Global.Channel].Incidence.retry)));
    process.Channel_Data[Global.Channel].Incidence.retry |= process.Channel_Data[Global.Channel].Incidence.aclip;
}//============================================================================================================================


//=============================================================================================================================
void Remove_Bits_Proc_Shaping_Off()
{//============================================================================================================================
    process.Channel_Data[Global.Channel].Count.eclips = 0;
    process.Channel_Data[Global.Channel].Count.sclips = 0;
    process.Channel_Data[Global.Channel].Count.rclips = 0;
    process.Channel_Data[Global.Channel].Count.aclips = 0;

    Removal.this_max_sample = sample_maxima[process.Channel_Data[Global.Channel].bits_removed];
    Removal.this_two_power = PowersOf.Two[process.Channel_Data[Global.Channel].bits_removed];
    Removal.this_two_power_recip = PowersOf.Two[-process.Channel_Data[Global.Channel].bits_removed];

    double this_DATA_sqr = 0;
    double this_LIMIT_sqr = 0;

    for (int32_t count = 0; count < AudioData.Size.This; ++count)
    {
        int32_t this_DATA = AudioData.WAVEPTR[THIS_CODEC_BLOCK][Global.Channel][count].Integers[0];

        bool extant_clip = ((this_DATA > Removal.this_max_sample) | (this_DATA == Removal.this_min_sample));

        process.Channel_Data[Global.Channel].Count.eclips += extant_clip;

        double scaled =  this_DATA * settings.scaling_factor;

        bool scaled_clip = ((scaled > Removal.this_max_sample) | (scaled < Removal.this_min_sample)) & (!extant_clip);

        process.Channel_Data[Global.Channel].Count.sclips += scaled_clip;

        int64_t this_BTRD = nRoundEvenInt64(scaled * Removal.this_two_power_recip) * Removal.this_two_power;

        int64_t this_LIMIT = Limited_Sample_Int64(this_BTRD);

        process.Channel_Data[Global.Channel].Count.rclips += ((this_LIMIT != this_BTRD) & (!scaled_clip));

        AudioData.BTRDPTR[THIS_CODEC_BLOCK][Global.Channel][count].Integers[0] = this_LIMIT;
        AudioData.CORRPTR[THIS_CODEC_BLOCK][Global.Channel][count].Integers[0] = this_DATA - nRoundEvenInt64(this_LIMIT * settings.scaling_factor_inv);

        this_DATA_sqr += (scaled * scaled);
        this_LIMIT_sqr += (this_LIMIT * this_LIMIT);
    }

    double this_round = (nlog2(this_LIMIT_sqr) - nlog2(this_DATA_sqr));

    process.Channel_Data[Global.Channel].Incidence.eclip = (process.Channel_Data[Global.Channel].Count.eclips > 0);
    process.Channel_Data[Global.Channel].Incidence.sclip = (process.Channel_Data[Global.Channel].Count.sclips > 0);

    process.Channel_Data[Global.Channel].Incidence.rclip = (process.Channel_Data[Global.Channel].Count.rclips > parameters.feedback.rclips);
    process.Channel_Data[Global.Channel].Incidence.retry = process.Channel_Data[Global.Channel].Incidence.rclip;

    process.Channel_Data[Global.Channel].Incidence.round = (parameters.feedback.active & (this_round > parameters.feedback.round) & (!process.Channel_Data[Global.Channel].Incidence.retry));
    process.Channel_Data[Global.Channel].Incidence.retry = process.Channel_Data[Global.Channel].Incidence.round;

    process.Channel_Data[Global.Channel].Incidence.noise = false;
    process.Channel_Data[Global.Channel].Incidence.aclip = false;
}//============================================================================================================================


//=============================================================================================================================
void Remove_Bits()
{//============================================================================================================================
    process.Channel_Data[Global.Channel].bits_removed = process.Channel_Data[Global.Channel].bits_to_remove;

    process.Channel_Data[Global.Channel].Total.eclip = 0;
    process.Channel_Data[Global.Channel].Total.sclip = 0;
    process.Channel_Data[Global.Channel].Total.rclip = 0;
    process.Channel_Data[Global.Channel].Total.aclip = 0;
    process.Channel_Data[Global.Channel].Total.noise = 0;
    process.Channel_Data[Global.Channel].Total.round = 0;

    process.Channel_Data[Global.Channel].Incidence.retry = true;

    while ((process.Channel_Data[Global.Channel].bits_removed >= 0) && (process.Channel_Data[Global.Channel].Incidence.retry))
    {
        if (process.Channel_Data[Global.Channel].bits_removed == 0)
            Store_Proc();
        else
            if (parameters.shaping.active)
                Remove_Bits_Proc_Adaptive_Noise_Shaping_On();
            else
                Remove_Bits_Proc_Shaping_Off();

        process.Channel_Data[Global.Channel].Total.eclip += process.Channel_Data[Global.Channel].Incidence.eclip;
        process.Channel_Data[Global.Channel].Total.sclip += process.Channel_Data[Global.Channel].Incidence.sclip;
        process.Channel_Data[Global.Channel].Total.rclip += process.Channel_Data[Global.Channel].Incidence.rclip;
        process.Channel_Data[Global.Channel].Total.aclip += process.Channel_Data[Global.Channel].Incidence.aclip;
        process.Channel_Data[Global.Channel].Total.noise += process.Channel_Data[Global.Channel].Incidence.noise;
        process.Channel_Data[Global.Channel].Total.round += process.Channel_Data[Global.Channel].Incidence.round;

        process.Channel_Data[Global.Channel].bits_removed -= (process.Channel_Data[Global.Channel].Incidence.retry);
    }

    if (process.Channel_Data[Global.Channel].bits_removed < 0)
        process.Channel_Data[Global.Channel].bits_removed = 0;

    process.Channel_Data[Global.Channel].bits_lost = process.Channel_Data[Global.Channel].bits_to_remove - process.Channel_Data[Global.Channel].bits_removed;

}//============================================================================================================================

//=============================================================================================================================
void nRemoveBits_Init()
{//============================================================================================================================
    Removal.this_max_sample = PowersOf.TwoM1[Global.bits_per_sample - 1];
    Removal.this_min_sample = -PowersOf.TwoInt32[Global.bits_per_sample - 1];


    for (int32_t fbi_i = 0; fbi_i < BITS_TO_CALCULATE; fbi_i++)
    {
        sample_maxima[fbi_i] = Removal.this_max_sample - (PowersOf.TwoM1[fbi_i]);
    }

    settings.static_maximum_bits_to_remove = std::max(0, Global.bits_per_sample - settings.static_minimum_bits_to_keep);
}//============================================================================================================================

