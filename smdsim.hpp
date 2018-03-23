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

#ifndef SMDSIM_HPP
#define SMDSIM_HPP

#include <vector>
#include <string>

#include "tools/num_methods/integrators.hpp"
#include "tools/plot_tools/figure.hpp"

struct ExactSolData
{
    double k = 0;
    double m = 0;
    double c = 0;

    double omega_n = 0;
    double omega_d = 0;
    double zeta = 0;

    double amp = 0;
    double phi = 0;

    double crit_c1 = 0;
    double crit_c2 = 0;

    double over_c1 = 0;
    double over_c2 = 0;

    double exp_power = 0;

};

class SMDSim
{
    using vec_size_t = std::vector<unsigned>::size_type;
    using rhsODE = double (*)(const std::vector<double> &state_vec, const SMDSim &sim_class);
    friend double x1Dot(const std::vector<double> &state_vec, const SMDSim &sim_class);
    friend double x2Dot(const std::vector<double> &state_vec, const SMDSim &sim_class);
public:
    SMDSim();
    ~SMDSim();

    void smdICSol(const std::string &damp_type);

private:
    // variables
    ExactSolData m_exact;
    std::vector<rhsODE> m_funcs;
    std::vector<double> m_state_vec;
    vec_size_t m_x1_ind = 0;
    vec_size_t m_x2_ind = 1;
    double m_dt;

    double m_mass = 0;
    double m_stiffness = 0;
    double m_damping = 0;

    std::string m_damp_type = "";

    // functions
    void initICSimulation(const std::string &damp_type);
    void exactICSol(std::vector<double> &x1_exact_sol, std::vector<double> &x2_diff_sol,
                  const std::vector<double> &t_vec, const std::string &damp_type);
    void setDamping(const std::string &damp_type);
    double underDampedSol(double t);
    double criticallyDampedSol(double t);
    double overDampedSol(double t);

    //plotting
    void plotData(const FigureData &plot_data, const std::string &file_name = "");

    void plotICsol(const std::vector<double>& x1_sol, const std::vector<double> &x2_sol,
                    const std::vector<double> &x1_exact_sol, const std::vector<double> &x2_diff_sol,
                    const std::vector<double> &t_vec);




};

// plant functions
double x1Dot(const std::vector<double> &state_vector, const SMDSim &sim_class);
double x2Dot(const std::vector<double> &state_vector, const SMDSim &sim_class);


#endif // SMDSIM_HPP
