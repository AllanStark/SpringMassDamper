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

#include <QVector>
#include <cmath>
#include <algorithm>
#include <functional>

#include "smdsim.hpp"
#include "tools/num_methods/numeric_algs.hpp"
#include "tools/num_methods/integrators.hpp"

using std::vector;
using std::sqrt;
using std::atan2;
using std::pow;
using std::exp;
using std::transform;
using std::minus;

SMDSim::SMDSim()
{
}

SMDSim::~SMDSim()
{
}

void SMDSim::smdICSol(const std::string &damp_type)
{
    //array of function pointers
    m_funcs.push_back(x1Dot);
    m_funcs.push_back(x2Dot);

    //initialize the simulation
    initICSimulation(damp_type);


    //Time vector
    m_dt = 0.001;
    double t_end = 3;
    vector<double> t_vec = timeVector(t_end, m_dt);

    //solution loop
    vector<double> x1_sol;
    vector<double> x2_sol;
    vector<double> temp;
    for (auto e : t_vec)
    {
        x1_sol.push_back(m_state_vec[m_x1_ind]);
        x2_sol.push_back(m_state_vec[m_x2_ind]);
        temp = rk4integrator(m_funcs, m_state_vec, m_dt, *this);
        m_state_vec = temp;
    }

    vector<double> x1_exact_sol;
    vector<double> x2_diff_sol;
    exactICSol(x1_exact_sol, x2_diff_sol, t_vec, damp_type);
    plotICsol(x1_sol, x2_sol, x1_exact_sol, x2_diff_sol, t_vec);

}

void SMDSim::initICSimulation(const std::string &damp_type)
{
    m_state_vec = {};

    //system parameters
    m_mass = 1;
    m_stiffness = 1000;
    m_damping = 0;

    //Initial conditions
    double x1_0 = 0.5;  //[m]
    double x2_0 = 0.0;  //[m/s]
    m_state_vec.push_back(x1_0);
    m_state_vec.push_back(x2_0);

    setDamping(damp_type);
}

/**
 * sets the damping value in m_sys_params to achieve the requested damping
 * @brief SMDSim::determinDamping
 * @param damp_type: Valid inputs: "under", "crit", "over" anything else will default to criticaly damped
 * @param m_sys_params
 */
void SMDSim::setDamping(const std::string &damp_type)
{
     m_exact.k = m_stiffness;
     m_exact.m = m_mass;
     m_exact.omega_n = std::sqrt(m_exact.k/m_exact.m);

    // find c so that zeta = 1;
    double c = 2 * m_exact.omega_n * m_exact.m;

    if (damp_type == "under"){
        m_damping = c / 10;
        m_damp_type = "Under_Damped";
    } else if (damp_type == "crit") {
        m_damping = c;
        m_damp_type = "Critically_Damped";
    } else if  (damp_type == "over") {
        m_damping = c * 2;
        m_damp_type = "Over_Damped";
    } else
        setDamping("crit");

    m_exact.c = m_damping;
}

/**
 * @brief SMDSim::plotData
 * @param plot_data
 */
void SMDSim::plotData(const FigureData &plot_data, const std::string &file_name)
{
        Figure* fig = new Figure();
        fig->plotData(plot_data);
        fig->show();
        if (!file_name.empty()){
            QString temp = QString::fromStdString(file_name);
            fig->saveImage(temp);
        }
}



/**
 * @brief SMDSim::exactICSol
 * @param x1_exact_sol
 * @param x2_diff_sol
 * @param t_vec
 */
void SMDSim::exactICSol(vector<double> &x1_exact_sol, vector<double> &x2_diff_sol,
                        const vector<double> &t_vec, const std::string &damp_type)
{
    initICSimulation(damp_type);
    double x0 = m_state_vec[m_x1_ind];
    double v0 = m_state_vec[m_x2_ind];
    m_exact.zeta = m_exact.c / (2*m_exact.omega_n*m_exact.m);
    m_exact.omega_d = m_exact.omega_n * sqrt(1 - m_exact.zeta * m_exact.zeta);

    double temp1 = x0 * m_exact.omega_d;
    double temp2 = v0 +x0*m_exact.zeta*m_exact.omega_n;

    double (SMDSim::*exactSol)(double);
    if (m_exact.zeta<1)
    {
        double temp = pow(temp1, 2) + pow(temp2, 2);
        m_exact.amp = sqrt(temp)/m_exact.omega_d;
        m_exact.phi = atan2(temp1, temp2);
        exactSol = &SMDSim::underDampedSol;
    }

    else if (m_exact.zeta == 1)
    {
        m_exact.crit_c1 = x0;
        m_exact.crit_c2 = temp2;
        exactSol = &SMDSim::criticallyDampedSol;
    }
    else
    {
        double temp = sqrt(pow(m_exact.zeta,2) - 1) ;
        double temp3 = temp + m_exact.zeta;
        temp3 = x0 * m_exact.omega_n * temp3 + v0;
        double temp4 = 2 * m_exact.omega_n * temp;
        m_exact.over_c1 = temp3 / temp4;

        temp3 = temp - m_exact.zeta;
        temp3 = x0 * m_exact.omega_n  * temp3 - v0;
        m_exact.over_c2 = temp3/temp4;

        exactSol = &SMDSim::overDampedSol;
    }

    for (auto e: t_vec)
    {
        m_exact.exp_power = -m_exact.zeta * m_exact.omega_n * e;
        double temp = (this->*exactSol)(e);
        x1_exact_sol.push_back(temp);
    }

    x2_diff_sol = derivative(x1_exact_sol, m_dt);
}


/**
 * @brief SMDSim::underDampedSol
 * @param t
 * @return
 */
double SMDSim::underDampedSol(double t)
{
    double temp = m_exact.omega_d * t + m_exact.phi;
    double x = m_exact.amp*exp(m_exact.exp_power)*sin(temp);
    return x;
}

/**
 * @brief SMDSim::criticallyDampedSol
 * @param t
 * @return
 */
double SMDSim::criticallyDampedSol(double t)
{
    double temp1 = m_exact.crit_c1*exp(m_exact.exp_power);
    double temp2 = m_exact.crit_c2*t*exp(m_exact.exp_power);
    double x = temp1 + temp2;
    return x;
}

/**
 * @brief SMDSim::overDampedSol
 * @param t
 * @return
 */
double SMDSim::overDampedSol(double t)
{
    double temp1 = exp(m_exact.exp_power);
    double temp2 = m_exact.omega_n * sqrt(pow(m_exact.zeta,2) - 1) * t;
    double temp3 = m_exact.over_c1 * exp(temp2);
    double temp4 = m_exact.over_c2 * exp(-temp2);
    double x = temp1 * (temp3  + temp4);
    return x;
}
// * @param x1_exact_sol


/**
 * @brief SMDSim::plotICsol
 * @param x1_sol
 * @param x2_sol
 * @param x1_exact_sol
 * @param x2_diff_sol
 * @param t_vec
 */
void SMDSim::plotICsol(const vector<double>& x1_sol, const vector<double> &x2_sol,
                    const vector<double> &x1_exact_sol, const vector<double> &x2_diff_sol,
                    const vector<double> &t_vec)
{
    FigureData plot_data1;
    plot_data1.addPlotData(x1_sol);
    plot_data1.addPlotData(x1_exact_sol);
    std::vector<double> error(x1_sol.size());
    transform(x1_sol.begin(), x1_sol.end(), x1_exact_sol.begin(), error.begin(), minus<double>());
    plot_data1.addPlotData(error);
    vector<vec_size_t> plot0_ind{0, 1};
    vector<vec_size_t> plot1_ind{2};
    plot_data1.addTVec(t_vec);
    plot_data1.addPlotInds(plot0_ind);
    plot_data1.addPlotInds(plot1_ind);
    plot_data1.addXLabel("");
    plot_data1.addXLabel("Time (s)");
    plot_data1.addYLabel("Pos (m)");
    plot_data1.addYLabel("Pos Error (m)");
    plot_data1.addLegend("Sim");
    plot_data1.addLegend("Exact");
    std::string zetaStr = toStringWPrecision(m_exact.zeta);
    std::string ttl1 = "Position Sol, \u03B6 = " + zetaStr;
    plot_data1.addTitle(ttl1);

    FigureData plot_data2;
    plot_data2.addPlotData(x2_sol);
    plot_data2.addPlotData(x2_diff_sol);
    plot_data2.addTVec(t_vec);
    error = std::vector<double>(x2_sol.size());
    transform(x1_sol.begin(), x1_sol.end(), x1_exact_sol.begin(), error.begin(), minus<double>());
    plot_data2.addPlotData(error);
    plot_data2.addPlotInds(plot0_ind);
    plot_data2.addPlotInds(plot1_ind);
    plot_data2.addXLabel("");
    plot_data2.addXLabel("Time (s)");
    plot_data2.addYLabel("Vel (m/s)");
    plot_data2.addYLabel("Vel Error (m/s)");
    plot_data2.addLegend("Sim");
    plot_data2.addLegend("Exact");
    std::string ttl2 = "Velocity Sol, \u03B6 = " + zetaStr;
    plot_data2.addTitle(ttl2);

    std::string path = "/home/james/Documents/SpringMassDamper/doc/images";
    std::string file_name1 = path + "/" + m_damp_type + "_Pos.png";
    std::string file_name2 = path + "/" + m_damp_type + "_Vel.png";

    plotData(plot_data1, file_name1);
    plotData(plot_data2, file_name2);
}

// Non class member double tfunctions
/**
 * @brief x1Dot
 * @param inputdouble t
 * @param state_vector
 * @param sys_params
 * @return
 */
double x1Dot(const vector<double> &state_vector, const SMDSim &sim_class)
{
    double x2 = state_vector[sim_class.m_x2_ind];
    double x1_dot = x2;
    return x1_dot;
}

/**
 * @brief x2Dotfriend x1Dot(const vector<double> &state_vec, const SMDSim &sim_class);
 * @param input
 * @param state_vector
 * @param sys_params
 * @return
 */
double x2Dot(const vector<double> &state_vector, const SMDSim &sim_class)
{
    double x1 = state_vector[sim_class.m_x1_ind];
    double x2 = state_vector[sim_class.m_x2_ind];

    double m = sim_class.m_mass;
    double c = sim_class.m_damping;
    double k = sim_class.m_stiffness;
    double F = 0;
    double x2_dot = - c/m * x2 -k/m * x1 + F/m;
    return x2_dot;
}
