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

#include "plotdata.hpp"
#include "figure.hpp"

void FigureData::addTVec(std::vector<double> t_vec)
{
    QVector<double> temp = QVector<double>::fromStdVector(t_vec);
    m_t_vec.push_back(temp);
}

void FigureData::addPlotData(std::vector<double> plot_data)
{
    QVector<double> temp = QVector<double>::fromStdVector(plot_data);
    m_plot_data.push_back(temp);
}

void FigureData::addPlotInds(std::vector<vec_size_t> plot_inds)
{
    QVector<vec_size_t> temp = QVector<vec_size_t>::fromStdVector(plot_inds);
    m_plot_inds.push_back(temp);
}

void FigureData::addLegend(std::string legend_data)
{
    QString temp = QString::fromStdString(legend_data);
    m_legend_data.push_back(temp);
}

void FigureData::addXLabel(std::string x_labels)
{
    QString temp = QString::fromStdString(x_labels);
    m_x_labels.push_back(temp);
}

void FigureData::addYLabel(std::string y_labels)
{
    QString temp = QString::fromStdString(y_labels);
    m_y_labels.push_back(temp);
}

void FigureData::addTitle(std::string title)
{
    m_title = QString::fromStdString(title);
}

