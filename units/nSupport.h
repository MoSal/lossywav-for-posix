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
//    by Tyge Løvset (tycho), Aug. 2012
//==============================================================================

#ifndef nSupport_h_
#define nSupport_h_

#include <ctime>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <fstream>

// DELPHI EMULATION

template <typename T, int32_t Low, int32_t High> struct range_array
{
    int32_t low() const { return Low; }
    int32_t high() const { return High; }
    size_t size() const { return High - Low + 1; }

    T& operator[](int32_t index) { return value[index - Low]; }
    const T& operator[](int32_t index) const { return value[index - Low]; }

    T value[High - Low + 1];
};


typedef unsigned char* PByte;


inline std::string IntToStr(int32_t number)
{
    std::ostringstream ss;
    ss << number;
    return ss.str();
}


inline std::string floattostrf(double number, int32_t digits)
{
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(digits) << std::setfill('0') << std::left << number;
    return ss.str();
}


inline bool FileIsReadOnly(const std::string& fileName_in)
{
    DWORD fAttr = GetFileAttributesA(fileName_in.c_str());
    if (fAttr == INVALID_FILE_ATTRIBUTES)
        return false;  //something is wrong with your path!
    if (fAttr & FILE_ATTRIBUTE_READONLY)
        return true;   // file is read only!
    return false;    // file is not read only!
}


inline bool DirectoryExists(const std::string& dirName_in)
{
    DWORD fAttrib = GetFileAttributesA(dirName_in.c_str());
    if (fAttrib == INVALID_FILE_ATTRIBUTES)
        return false;  //something is wrong with your path!
    if (fAttrib & FILE_ATTRIBUTE_DIRECTORY)
        return true;   // this is a directory!
    return false;    // this is not a directory! } {
}


inline bool FileExists(std::string file_name)
{
    std::fstream filetest(file_name.c_str(), std::ios_base::in | std::ios_base::binary);
    if (filetest.is_open())
    {
        filetest.close();
        return true;
    }
    return false;
}

#endif // nSupport_h
