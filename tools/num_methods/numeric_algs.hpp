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

#ifndef NUMERIC_ALGS_HPP
#define NUMERIC_ALGS_HPP
#include <vector>
#include <string>
#include <iomanip>
#include <sstream>


std::vector<double> derivative(const std::vector<double> &signal, const double &dt);
std::vector<double> sweptSine(std::vector<double> t_vec, double f2, double amp);
std::vector<double> timeVector(double t_end, double dt);

template <typename T>
std::string toStringWPrecision(const T val, const int n = 6) {
    std::ostringstream out;
    out << std::setprecision(n) << val;
    return out.str();
}


#endif // NUMERIC_ALGS_HPP

