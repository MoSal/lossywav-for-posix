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

#include "nRemoveBits.h"

struct Removal_Type
{
    double this_two_power;
    double this_two_power_recip;
    int64_t this_min_sample;
    int64_t this_max_sample;

};

Removal_Type RemovalBits[BITS_TO_CALCULATE + 1] __attribute__ ((aligned(16)));

Channel_Data_Type* this_channel_data;

Removal_Type* this_removal;


//=============================================================================================================================
int64_t Limited_Sample_Int64(int64_t lsv)
//=============================================================================================================================
{
    return lsv + ((this_removal->this_max_sample - lsv) & -(this_removal->this_max_sample < lsv)) + ((this_removal->this_min_sample - lsv) & -(this_removal->this_min_sample > lsv));
}
//=============================================================================================================================


//=============================================================================================================================
void Store_Proc()
{//============================================================================================================================
    this_channel_data->Count.eclips = 0;
    this_channel_data->Count.sclips = 0;
    this_channel_data->Count.rclips = 0;
    this_channel_data->Count.aclips = 0;

    for (int32_t count = 0; count < AudioData.Size.This; ++count)
    {
        int64_t this_DATA = AudioData.WAVEPTR[THIS_CODEC_BLOCK][Current.Channel][count].Integers[0];

        bool extant_clip = ((this_DATA > this_removal->this_max_sample) | (this_DATA == this_removal->this_min_sample));

        this_channel_data->Count.eclips += extant_clip;

        double scaled =  this_DATA * settings.scaling_factor;

        bool scaled_clip = ((scaled > this_removal->this_max_sample) | (scaled < this_removal->this_min_sample)) & (!extant_clip);

        this_channel_data->Count.sclips += scaled_clip;

        int64_t this_BTRD = nRoundEvenInt64(scaled);

        int64_t this_LIMIT = Limited_Sample_Int64(this_BTRD);

        this_channel_data->Count.rclips += ((this_LIMIT != this_BTRD) & (!scaled_clip));

        AudioData.BTRDPTR[THIS_CODEC_BLOCK][Current.Channel][count].Integers[0] = this_LIMIT;
        AudioData.CORRPTR[THIS_CODEC_BLOCK][Current.Channel][count].Integers[0] = this_DATA - nRoundEvenInt64(this_LIMIT * settings.scaling_factor_inv);
    }

    double this_round = 0;

    this_channel_data->Incidence.eclip = (this_channel_data->Count.eclips > 0);
    this_channel_data->Incidence.sclip = (this_channel_data->Count.sclips > 0);

    this_channel_data->Incidence.rclip = (this_channel_data->Count.rclips > parameters.feedback.rclips);
    this_channel_data->Incidence.retry = this_channel_data->Incidence.rclip;

    this_channel_data->Incidence.round = (parameters.feedback.active & (this_round > parameters.feedback.round) & (!this_channel_data->Incidence.retry));
    this_channel_data->Incidence.retry = this_channel_data->Incidence.round;

    this_channel_data->Incidence.noise = false;
    this_channel_data->Incidence.aclip = false;
}//============================================================================================================================


//=============================================================================================================================
void Remove_Bits_Proc_Adaptive_Noise_Shaping_On()
{//============================================================================================================================
    double this_noise = 0;
    double this_round = 0;

    Warped_Lattice_Filter_Init(Current.Channel);

    this_channel_data->Count.eclips = 0;
    this_channel_data->Count.sclips = 0;
    this_channel_data->Count.rclips = 0;
    this_channel_data->Count.aclips = 0;

    double this_DATA_sqr = 0;
    double this_LIMIT_sqr = 0;
    double this_WLFE_sqr = 0;

    double limit_WLFE = npower2(AudioData.Channel_Log2_RMS[Current.Channel] + parameters.feedback.alevel);

    for (int32_t count = 0; count < AudioData.Size.This; ++count)
    {
        int64_t this_DATA = AudioData.WAVEPTR[THIS_CODEC_BLOCK][Current.Channel][count].Integers[0];

        bool extant_clip = (this_DATA > this_removal->this_max_sample) | (this_DATA < this_removal->this_min_sample);
        this_channel_data->Count.eclips += extant_clip;

        double scaled =  this_DATA * settings.scaling_factor;

        bool scaled_clip = (scaled > this_removal->this_max_sample) | (scaled < this_removal->this_min_sample);
        this_channel_data->Count.sclips += scaled_clip;

        double this_WLFE = Warped_Lattice_Filter_Evaluate(Current.Channel);
        this_WLFE_sqr += (this_WLFE * this_WLFE);
        this_channel_data->Count.aclips += (fabs(this_WLFE) > limit_WLFE);

        int64_t this_SHAPED = nRoundEvenInt64((scaled + this_WLFE) * this_removal->this_two_power_recip) * this_removal->this_two_power;

        int64_t this_LIMIT = Limited_Sample_Int64(this_SHAPED);

        this_channel_data->Count.rclips += ((this_LIMIT != this_SHAPED) & (!scaled_clip));

        Warped_Lattice_Filter_Update(Current.Channel, this_LIMIT - scaled);

        AudioData.BTRDPTR[THIS_CODEC_BLOCK][Current.Channel][count].Integers[0] = this_LIMIT;
        AudioData.CORRPTR[THIS_CODEC_BLOCK][Current.Channel][count].Integers[0] = this_DATA - nRoundEvenInt64(this_LIMIT * settings.scaling_factor_inv);

        this_DATA_sqr += (scaled * scaled);
        this_LIMIT_sqr += (this_LIMIT * this_LIMIT);
    }

    this_noise = (0.5 * nlog2(this_WLFE_sqr)) - Global.bits_per_sample;

    this_round = (nlog2(this_LIMIT_sqr) - nlog2(this_DATA_sqr));

    this_channel_data->Incidence.eclip = (this_channel_data->Count.eclips > 0);
    this_channel_data->Incidence.sclip = (this_channel_data->Count.sclips > 0);

    this_channel_data->Incidence.rclip = (this_channel_data->Count.rclips > parameters.feedback.rclips);
    this_channel_data->Incidence.retry = this_channel_data->Incidence.rclip;

    this_channel_data->Incidence.round = (parameters.feedback.active & (this_round > parameters.feedback.round) & (!this_channel_data->Incidence.retry));
    this_channel_data->Incidence.retry |= this_channel_data->Incidence.round;

    this_channel_data->Incidence.noise = (parameters.feedback.active & (this_noise > parameters.feedback.noise) & (!this_channel_data->Incidence.retry));
    this_channel_data->Incidence.retry |= this_channel_data->Incidence.noise;

    this_channel_data->Incidence.aclip = (parameters.feedback.active & (this_channel_data->Count.aclips > parameters.feedback.aclips) & (!(this_channel_data->Incidence.retry)));
    this_channel_data->Incidence.retry |= this_channel_data->Incidence.aclip;
}//============================================================================================================================


//=============================================================================================================================
void Remove_Bits_Proc_Shaping_Off()
{//============================================================================================================================
    this_channel_data->Count.eclips = 0;
    this_channel_data->Count.sclips = 0;
    this_channel_data->Count.rclips = 0;
    this_channel_data->Count.aclips = 0;

    double this_DATA_sqr = 0;
    double this_LIMIT_sqr = 0;

    for (int32_t count = 0; count < AudioData.Size.This; ++count)
    {
        int64_t this_DATA = AudioData.WAVEPTR[THIS_CODEC_BLOCK][Current.Channel][count].Integers[0];

        bool extant_clip = ((this_DATA > this_removal->this_max_sample) | (this_DATA == this_removal->this_min_sample));

        this_channel_data->Count.eclips += extant_clip;

        double scaled =  this_DATA * settings.scaling_factor;

        bool scaled_clip = ((scaled > this_removal->this_max_sample) | (scaled < this_removal->this_min_sample)) & (!extant_clip);

        this_channel_data->Count.sclips += scaled_clip;

        int64_t this_BTRD = nRoundEvenInt64(scaled * this_removal->this_two_power_recip) * this_removal->this_two_power;

        int64_t this_LIMIT = Limited_Sample_Int64(this_BTRD);

        this_channel_data->Count.rclips += ((this_LIMIT != this_BTRD) & (!scaled_clip));

        AudioData.BTRDPTR[THIS_CODEC_BLOCK][Current.Channel][count].Integers[0] = this_LIMIT;
        AudioData.CORRPTR[THIS_CODEC_BLOCK][Current.Channel][count].Integers[0] = this_DATA - nRoundEvenInt64(this_LIMIT * settings.scaling_factor_inv);

        this_DATA_sqr += (scaled * scaled);
        this_LIMIT_sqr += (this_LIMIT * this_LIMIT);
    }

    double this_round = (nlog2(this_LIMIT_sqr) - nlog2(this_DATA_sqr));

    this_channel_data->Incidence.eclip = (this_channel_data->Count.eclips > 0);
    this_channel_data->Incidence.sclip = (this_channel_data->Count.sclips > 0);

    this_channel_data->Incidence.rclip = (this_channel_data->Count.rclips > parameters.feedback.rclips);
    this_channel_data->Incidence.retry = this_channel_data->Incidence.rclip;

    this_channel_data->Incidence.round = (parameters.feedback.active & (this_round > parameters.feedback.round) & (!this_channel_data->Incidence.retry));
    this_channel_data->Incidence.retry = this_channel_data->Incidence.round;

    this_channel_data->Incidence.noise = false;
    this_channel_data->Incidence.aclip = false;
}//============================================================================================================================


//=============================================================================================================================
void Remove_Bits()
{//============================================================================================================================
    this_channel_data = &process.Channel_Data[Current.Channel];

    this_channel_data->bits_removed = this_channel_data->bits_to_remove;

    this_channel_data->Total.eclip = 0;
    this_channel_data->Total.sclip = 0;
    this_channel_data->Total.rclip = 0;
    this_channel_data->Total.aclip = 0;
    this_channel_data->Total.noise = 0;
    this_channel_data->Total.round = 0;

    this_channel_data->Incidence.retry = true;

    while ((this_channel_data->bits_removed >= 0) && (this_channel_data->Incidence.retry))
    {
        this_removal = &RemovalBits[this_channel_data->bits_removed];

        if (this_channel_data->bits_removed == 0)
            Store_Proc();
        else
            if (parameters.shaping.active)
                Remove_Bits_Proc_Adaptive_Noise_Shaping_On();
            else
                Remove_Bits_Proc_Shaping_Off();

        this_channel_data->Total.eclip += this_channel_data->Incidence.eclip;
        this_channel_data->Total.sclip += this_channel_data->Incidence.sclip;
        this_channel_data->Total.rclip += this_channel_data->Incidence.rclip;
        this_channel_data->Total.aclip += this_channel_data->Incidence.aclip;
        this_channel_data->Total.noise += this_channel_data->Incidence.noise;
        this_channel_data->Total.round += this_channel_data->Incidence.round;

        this_channel_data->bits_removed -= (this_channel_data->Incidence.retry);
    }

    if (this_channel_data->bits_removed < 0)
        this_channel_data->bits_removed = 0;

    this_channel_data->bits_lost = this_channel_data->bits_to_remove - this_channel_data->bits_removed;

}//============================================================================================================================

//=============================================================================================================================
void nRemoveBits_Init()
{//============================================================================================================================
    for (int32_t fbi_i = 0; fbi_i < (Global.bits_per_sample - 1); fbi_i++)
    {
        RemovalBits[fbi_i].this_max_sample = PowersOf.TwoM1[Global.bits_per_sample - 1] - (PowersOf.TwoM1[fbi_i]);
        RemovalBits[fbi_i].this_min_sample = -PowersOf.TwoInt32[Global.bits_per_sample - 1];
        RemovalBits[fbi_i].this_two_power = PowersOf.TwoX[TWO_OFFSET + fbi_i];
        RemovalBits[fbi_i].this_two_power_recip = PowersOf.TwoX[TWO_OFFSET + -fbi_i];
    }

    settings.static_maximum_bits_to_remove = std::max(0, Global.bits_per_sample - settings.static_minimum_bits_to_keep);
}//============================================================================================================================
