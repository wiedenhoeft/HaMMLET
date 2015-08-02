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


#include "Uniform.hpp"

double Uniform::pdf(double *obs, size_t n){
    double x = *obs;
    if(x >= a && x <= b){
        return 1/(b-a);
    }
    else{
        return 0;
    }
}

double Uniform::logPdf(double *obs, size_t n){
    return log(this->pdf(obs, n));
}

double Uniform::pdf(WaveletNode *obs, size_t n){
    return pow(1/(b-a), n);
}
double Uniform::logPdf(WaveletNode *obs, size_t n){
    return n*log(1/(b-a));
}

