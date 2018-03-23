/***************************************************************************
**                                                                        **
**  Spring Mass Damper vibration solution            					  **
**  Copyright (C) 2015 James Longino 						              **
**                                                                        **
**  This program is free software: you can redistribute it and/or modify  **
**  it under the terms of the GNU General Public License as published by  **
**  the Free Software Foundation, either version 3 of the License, or     **
**  (at your option) any later version.                                   **
**                                                                        **
**  This program is distributed in the hope that it will be useful,       **
**  but WITHOUT ANY WARRANTY; without even the implied warranty of        **
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         **
**  GNU General Public License for more details.                          **
**                                                                        **
**  You should have received a copy of the GNU General Public License     **
**  along with this program.  If not, see http://www.gnu.org/licenses/.   **
**                                                                        **
****************************************************************************
**           Author: James Longino                                   	  **
**  Contact: shot5114@gmail.com                        					  ** **						                                                  **
****************************************************************************/

#define _USE_MATH_DEFINES
#include <cmath>
#include "numeric_algs.hpp"

using std::vector;
using std::sin;

std::vector<double> derivative(const std::vector<double> &signal,const double &dt)
{
    vector<double>::size_type n = signal.size();
    vector<double> signal_prime;
    for (unsigned i = 0; i != n; ++i)
    {
        if (i !=n-1)
            signal_prime.push_back((signal[i+1] - signal[i])/dt);
        else
            // the last value of the derivative of the signal  is undefined so instead duplicate the final value
            signal_prime.push_back(signal_prime[i-1]);
    }
    return signal_prime;


}


vector<double> sweptSine(vector<double> t_vec, double f2, double amp)
{
    auto endIt = t_vec.end();
    --endIt;
    double t_end = *endIt;
    vector<double> y;
    for (auto e : t_vec){
        double temp = (M_PI * f2 * e*e) / t_end;
        temp = sin(temp);
        y.push_back(amp * temp);
    }

    return y;
}


std::vector<double> timeVector(double t_end, double dt)
{
    double fs = 1/dt;
    t_end =t_end * fs;
    std::vector<double> t_vec;
    double temp = 0;

    // to partially avoid round off error count by ones and multiply by time step
    for (int i = 0; i <= t_end; ++i) {
        temp = i * dt;
        t_vec.push_back(temp);
    }
    return t_vec;
}
