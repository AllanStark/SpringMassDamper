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

#ifndef INTEGRATORS_HPP
#define INTEGRATORS_HPP
#include <vector>

template <typename tF, typename C>
std::vector<double> rk4integrator(std::vector<tF> funcs, const std::vector<double> &state_vec,
                                  const double &dt, const C &sim_class){
    using vec_size_t = std::vector<double>::size_type;
    vec_size_t n = funcs.size();

    // calc K1 values
    std::vector<double> k1;
    for (auto e : funcs) {
        k1.push_back(e(state_vec, sim_class));
    }

    // calc K2 values
    std::vector<double> k2;
    std::vector<double> temp;
    for (vec_size_t i = 0; i != n; ++i) {
        double state_val_temp = state_vec[i] + dt/2 * k1[i];
        temp.push_back(state_val_temp);
    }

    for (auto e : funcs) {
        k2.push_back(e(temp, sim_class));
    }

    // calc K3 values
    std::vector<double> k3;
    temp = {};
    for (vec_size_t i = 0; i != n; ++i) {
        double state_val_temp = state_vec[i] + dt/2 * k2[i];
        temp.push_back(state_val_temp);
    }

    for (auto e : funcs) {
        k3.push_back(e(temp, sim_class));
    }

    // calc K4 values
    std::vector<double> k4;
    temp = {};
    for (vec_size_t i = 0; i != n; ++i) {
        double state_val_temp = state_vec[i] + dt * k3[i];
        temp.push_back(state_val_temp);
    }

    for (auto e : funcs) {
        k4.push_back(e(temp, sim_class));
    }


    // Calc new state vector
    std::vector<double> state_vec_new;
    for (vec_size_t i = 0; i != n; ++i) {
        double state_val_temp = state_vec[i] + dt/6 *(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
        state_vec_new.push_back(state_val_temp);
    }

    return state_vec_new;
}




#endif // INTEGRATORS_HPP

