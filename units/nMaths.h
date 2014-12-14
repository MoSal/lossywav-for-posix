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

#ifndef nMaths_h_
#define nMaths_h_

#include <cmath>
#include <stdlib.h>
#include "nCore.h"

//=====================================================================

static const double log2_10_over_10 = log2_10 / 10.0d;
static const double log2_10_over_20 = log2_10 / 20.0d;

inline int32_t nRoundEvenInt32(double val)
{
    int32_t temp_int = std::abs(val);
    double temp_frac = (std::abs(val) - temp_int);
    temp_frac += temp_frac;

    return int32_t     (temp_frac == 1.0d ?
                   (val > 0 ? temp_int + (temp_int & 1) : -(temp_int + (temp_int & 1))) :
                   (val > 0 ? temp_int + int(temp_frac) : -(temp_int + int(temp_frac))));
}

inline int32_t nRoundOddInt32(double val)
{
    int32_t temp_int = std::abs(val);
    double temp_frac = (std::abs(val) - temp_int);
    temp_frac += temp_frac;

    return int32_t     (temp_frac == 1.0d ?
                   (val > 0 ? temp_int + ((temp_int & 1) ^ 1) : -(temp_int + ((temp_int & 1) ^ 1))) :
                   (val > 0 ? temp_int + int(temp_frac) : -(temp_int + int(temp_frac))));
}


//=====================================================================


inline int64_t nRoundEvenInt64(double val)
{
    int64_t temp_int = std::abs(val);
    double temp_frac = (std::abs(val) - temp_int);
    temp_frac += temp_frac;
    return int64_t (temp_frac == 1.0d ?
                   (val > 0 ? temp_int + (temp_int & 1) : -(temp_int + (temp_int & 1))) :
                   (val > 0 ? temp_int + int(temp_frac) : -(temp_int + int(temp_frac))));
}


inline int64_t nRoundOddInt64(double val)
{
    int64_t temp_int = std::abs(val);
    double temp_frac = (std::abs(val) - temp_int);
    temp_frac += temp_frac;

    return int64_t (temp_frac == 1.0d ?
                   (val > 0 ? temp_int + ((temp_int & 1) ^ 1) : -(temp_int + ((temp_int & 1) ^ 1))) :
                   (val > 0 ? temp_int + int(temp_frac) : -(temp_int + int(temp_frac))));
}


//=====================================================================


inline double nlog2(int32_t fl_2v)
{
    return (fl_2v > 0 ? log2(fl_2v) : 0);
}

inline double nlog2(int64_t fl_2v)
{
    return (fl_2v > 0 ? log2(fl_2v) : 0);
}

inline double nlog2(double fl_2v)
{
    return (fl_2v > 0 ? log2(fl_2v) : 0);
}


//=====================================================================


inline double nlog10(int32_t fl_10v)
{
    return (fl_10v > 0 ? log2(fl_10v) * log10_2 : 0);
}

inline double nlog10(int64_t fl_10v)
{
    return (fl_10v > 0 ? log2(fl_10v) * log10_2 : 0);
}

inline double nlog10(double fl_10v)
{
    return (fl_10v > 0 ? log2(fl_10v) * log10_2 : 0);
}


//=====================================================================


inline double nroot(int32_t fr_v)
{
    return (fr_v > 0 ? sqrt(fr_v) : 0);
}

inline double nroot(int64_t fr_v)
{
    return (fr_v > 0 ? sqrt(fr_v) : 0);
}

inline double nroot(double fr_v)
{
    return (fr_v > 0 ? sqrt(fr_v) : 0);
}


//=====================================================================


inline double nsqrd(double val)
{
    return val * val;
}


//=====================================================================


inline double nround10(double nr_v, int32_t nr_d)
{
    return (fabs(nr_d) > 12 ? nr_v : nRoundEvenInt64(nr_v * PowersOf.Ten[nr_d]) * PowersOf.Ten[-nr_d]);
}


//=====================================================================


inline double npower2(double fp_2v)
{
    return pow(2, fp_2v);
}


//=====================================================================


inline double npowerx(double fp_2v, double fp_2p)
{
    return pow(fp_2v, fp_2p);
}


//=====================================================================


struct CubicInterp_Rec
{
    double y0, y1, y2, y3, x, last;
};

inline double CubicInterp(CubicInterp_Rec* CubicData)
{
    double x_p2 = CubicData->x * CubicData->x;

//    double r0 = -0.5 * CubicData->y0 + 1.5 * CubicData->y1 - 1.5 * CubicData->y2 + 0.5 * CubicData->y3;
//    double r1 = CubicData->y0 - 2.5 * CubicData->y1 + 2 * CubicData->y2 - 0.5 * CubicData->y3;
//    double r2 = -0.5 * CubicData->y0 + 0.5 * CubicData->y2;
//    double r3 = CubicData->y1;

//    CubicData->last = (r0 * (CubicData->x * x_p2) + r1 * x_p2 + r2 * CubicData->x + r3);

//    return CubicData->last;

    //double a0 = CubicData->y3 - CubicData->y2 - CubicData->y0 + CubicData->y1;
    //double a1 = CubicData->y0 - CubicData->y1 - a0;
    //double a2 = CubicData->y2 - CubicData->y0;
    //return CubicData->last = a0 * CubicData->x * x_p2 + a1 * x_p2 + a2 * CubicData->x + CubicData->y1;

    return CubicData->last = (CubicData->y3 - CubicData->y2 - CubicData->y0 + CubicData->y1) * CubicData->x * x_p2 + (CubicData->y0 - CubicData->y3 + CubicData->y2 + CubicData->y0) * x_p2 + (CubicData->y2 - CubicData->y0) * CubicData->x + CubicData->y1;
//    return CubicData->last = CubicData->y1 + 0.5d * (((((3 * (CubicData->y1 - CubicData->y2) + CubicData->y3 - CubicData->y0) * CubicData->x) + (CubicData->y0 + CubicData->y0 - 4 * (CubicData->y1 - CubicData->y2) - CubicData->y1 - CubicData->y3)) * CubicData->x + (CubicData->y2 - CubicData->y0)) * CubicData->x);
}


inline double CubicInterpMaxLimit(CubicInterp_Rec* CubicData)
{
    double x_p2 = CubicData->x * CubicData->x;

//    double r0 = -0.5 * CubicData->y0 + 1.5 * CubicData->y1 - 1.5 * CubicData->y2 + 0.5 * CubicData->y3;
//    double r1 = CubicData->y0 - 2.5 * CubicData->y1 + 2 * CubicData->y2 - 0.5 * CubicData->y3;
//    double r2 = -0.5 * CubicData->y0 + 0.5 * CubicData->y2;
//    double r3 = CubicData->y1;
//    CubicData->last = (r0 * (CubicData->x * x_p2) + r1 * x_p2 + r2 * CubicData->x + r3);
//    return CubicData->last;

//    double a0 = CubicData->y3 - CubicData->y2 - CubicData->y0 + CubicData->y1;
//    double a1 = CubicData->y0 - CubicData->y1 - a0;
//    double a2 = CubicData->y2 - CubicData->y0;
//    CubicData->last = a0 * CubicData->x * x_p2 + a1 * x_p2 + a2 * CubicData->x + CubicData->y1;

    CubicData->last = (CubicData->y3 - CubicData->y2 - CubicData->y0 + CubicData->y1) * CubicData->x * x_p2 + (CubicData->y0 - CubicData->y3 + CubicData->y2 + CubicData->y0) * x_p2 + (CubicData->y2 - CubicData->y0) * CubicData->x + CubicData->y1;

//    double rx = (CubicData->y1 - CubicData->y2);
//    double r0 = ((3.0d * rx + CubicData->y3 - CubicData->y0) * CubicData->x);
//    double r1 = (CubicData->y0 + CubicData->y0 - 4.0d * rx - CubicData->y1 - CubicData->y3);
//    CubicData->last = CubicData->y1 + 0.5d * ((r0 + r1) * x_p2 + (CubicData->y2 - CubicData->y0) * CubicData->x);

//    CubicData->last = CubicData->y1 + 0.5d * (((((3 * (CubicData->y1 - CubicData->y2) + CubicData->y3 - CubicData->y0) * CubicData->x) + (CubicData->y0 + CubicData->y0 - 4 * (CubicData->y1 - CubicData->y2) - CubicData->y1 - CubicData->y3)) * CubicData->x + (CubicData->y2 - CubicData->y0)) * CubicData->x);

    return CubicData->last -= (CubicData->last - CubicData->y2) * (CubicData->last > CubicData->y2);
}


//=====================================================================

inline double dB_Amplitude_Ratio(double dB_Value)
{
    return pow(2, dB_Value * log2_10_over_10);
}

inline double dB_Power_Ratio(double dB_Value)
{
    return pow(2, dB_Value * log2_10_over_20);
}

inline double Bark(double Frequency_Value)
{
    if (Frequency_Value < 421.6d)
        return Frequency_Value * 0.01d;
    else
        return (26.81d / (1.0d + (1960.0d / Frequency_Value)) - 0.53d);
//    return 13.0d * atan(0.00076d * frequency) + 3.5d * atan(nsqrd(frequency / 7500.0d));
}

inline double Frequency_from_Bark(double Bark_Value)
{
    if (Bark_Value < 4.216d)
        return Bark_Value * 100.0d;
    else
        return (1960.0d / (26.81d / (Bark_Value + 0.53d) - 1.0d));
}

inline double Critical_Bandwidth(double Bark_Value)
{
    return 52548.0d / ((Bark_Value * Bark_Value) - (52.56d * Bark_Value) + 690.39d);
}

#endif // nMaths_h_
