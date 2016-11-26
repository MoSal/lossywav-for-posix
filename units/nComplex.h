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

===========================================================================**/

#ifndef nComplex_h_
#define nComplex_h_

struct tDComplex
{
    double Re, Im;

    inline tDComplex operator + (int32_t A)
    {
        tDComplex result;
        result.Re = Re + A;
        result.Im = Im;
        return result;
    }

    inline tDComplex operator + (double A)
    {
        tDComplex result;
        result.Re = Re + A;
        result.Im = Im;
        return result;
    }

    inline tDComplex operator + (tDComplex A)
    {
        tDComplex result;
        result.Re = Re + A.Re;
        result.Im = Im + A.Im;
        return result;
    }

    inline tDComplex& operator += (tDComplex A)
    {
        Re += A.Re;
        Im += A.Im;
        return *this;
    }

    inline tDComplex operator - (int32_t A)
    {
        tDComplex result;
        result.Re = Re - A;
        result.Im = Im;
        return result;
    }

    inline tDComplex operator - (double A)
    {
        tDComplex result;
        result.Re = Re - A;
        result.Im = Im;
        return result;
    }

    inline tDComplex operator - (tDComplex A)
    {
        tDComplex result;
        result.Re = Re - A.Re;
        result.Im = Im - A.Im;
        return result;
    }

    inline tDComplex& operator -= (tDComplex A)
    {
        Re -= A.Re;
        Im -= A.Im;
        return *this;
    }

    inline tDComplex operator * (tDComplex A)
    {
        tDComplex result;
        result.Re = Re * A.Re - Im * A.Im;
        result.Im = Re * A.Im + Im * A.Re;
        return result;
    }

    inline tDComplex& operator *= (tDComplex A)
    {
        double Wv = Re;
        Re = Re * A.Re - Im * A.Im;
        Im = Wv * A.Im + Im * A.Re;
        return *this;
    }

    inline tDComplex operator * (int32_t b)
    {
        tDComplex result;
        result.Re = Re * b;
        result.Im = Im * b;
        return result;
    }

    inline tDComplex operator * (double b)
    {
        tDComplex result;
        result.Re = Re * b;
        result.Im = Im * b;
        return result;
    }

    inline tDComplex& operator *= (int32_t A)
    {
        Re *= A;
        Im *= A;
        return *this;
    }

    inline tDComplex& operator *= (double A)
    {
        Re *= A;
        Im *= A;
        return *this;
    }

    inline tDComplex operator / (tDComplex b)
    {

        double Wv = b.Re * b.Re + b.Im * b.Im;

        if (Wv != 0)
        {
            tDComplex result;
            Wv = 1 / Wv;
            result.Re = (Re * b.Re + Im * b.Im) * Wv;
            result.Im = (Im * b.Re - Re * b.Im) * Wv;
            return result;
        }
        else
        {
            // return complex zero - not "correct". Should throw an exception.
            tDComplex result;
            result.Re = 0;
            result.Im = 0;
            return result;
        }
    }

    inline tDComplex& operator /= (tDComplex b)
    {
        double Wv = b.Re * b.Re + b.Im * b.Im;

        if (Wv != 0)
        {
            double Wx = Re;
            Wv = 1 / Wv;
            Re = (Re * b.Re + Im * b.Im) * Wv;
            Im = (Im * b.Re - Wx * b.Im) * Wv;
            return *this;
        }
        else
        {
            // return complex zero - not "correct". Should throw an exception.
            Re = 0;
            Im = 0;
            return *this;
        }
    }

    inline double magsquared()
    {
        return (Re * Re + Im * Im);
    }

    inline double magnitude()
    {
        return sqrt(Re * Re + Im * Im);
    }

    inline tDComplex conj()
    {
        tDComplex result;
        result.Re = Re;
        result.Im = -Im;
        return result;
    }

    inline tDComplex divided_by_i()
    {
        tDComplex result;
        result.Re = Im;
        result.Im = -Re;
        return result;
    }

    inline tDComplex multiplied_by_i()
    {
        tDComplex result;
        result.Re = -Im;
        result.Im = Re;
        return result;
    }

};

inline tDComplex DComplex(double a, double b)
{
    tDComplex result;
    result.Re = a;
    result.Im = b;
    return result;
}

inline tDComplex complex_exp(double a)
{
    tDComplex result;
#ifdef HAVE_SINCOS
    sincos(a, &result.Im, &result.Re);
#else
    result.Im = sin(a);
    result.Re = cos(a);
#endif
    return result;
}

#endif // nComplex_h_
