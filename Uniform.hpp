//     Copyright 2014-2015 John Wiedenhoeft, Eric Brugel
//
//     This file is part of HaMMLET.
// 
//     HaMMLET is free software: you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.
// 
//     HaMMLET is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
// 
//     You should have received a copy of the GNU General Public License
//     along with HaMMLET.  If not, see <http://www.gnu.org/licenses/>.


#ifndef UNIFORM_HPP
#define UNIFORM_HPP

#include "Distribution.hpp"

class Uniform: public Distribution{
    double a, b;
    public:

    Uniform(double a, double b):a(a), b(b){}

    double pdf(double *obs, size_t n);
    double logPdf(double *obs, size_t n);
    double pdf(WaveletNode *obs, size_t n);
    double logPdf(WaveletNode *obs, size_t n);
};

#endif
