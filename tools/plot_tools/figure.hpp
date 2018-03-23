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

#ifndef FIGURE_H
#define FIGURE_H

#include <QObject>
#include <QWidget>

#include "qcustomplot.h"

class FigureData
{
/* Data structure to control plotting using the figure class
 *
 * t_vec: vector for xaxis data.  If it only has one vector element all vectors in plot_data will be plotted vs that
 * one vector. If there are multiple vector elements there must be the same number of elements as in plot_data and they
 * will be plotted vs each other.
 *
 * plot_data: vector of vector elements.  Length of the inner vector elements must match the length of the
 * t_vec vector elements
 *
 * plot_inds: each inner vector element holds the indecies of plot_data (and t_vec if multiple time vectors) to be
 * plotted in the subplots. ie. the indicies in the first element will be used to plot data in the first subplot.
 * Number of outer vector elements determines the number of subplots.
 *
 * legend_data: one entry for each element in plot_data, used to build the legend that will be shown
 *
 * x_lables, y_lables: x label and ylabel used to annotate the subplots
 *
 * title: sting used to title the plot
 *
 * use of add functions will add the data to the approprate vector and tansition from std library to QT
 * library types
 */
    using vec_size_t = std::vector<unsigned>::size_type;
    friend class Figure;

public:
    void addTVec(std::vector<double> t_vec);
    void addPlotData(std::vector<double> plot_data);
    void addPlotInds(std::vector<vec_size_t> plot_inds);
    void addLegend(std::string legend_data);
    void addXLabel(std::string x_labels);
    void addYLabel(std::string y_labels);
    void addTitle(std::string title);
private:
    using vec_size = QList<unsigned>::size_type;
    QVector<QVector<double>> m_t_vec;
    QVector<QVector<double>> m_plot_data;
    QVector<QVector<vec_size_t>> m_plot_inds;
    QVector<QString> m_legend_data;
    QVector<QString> m_x_labels;
    QVector<QString> m_y_labels;
    QString m_title;
};

class Figure : public QWidget
{
    using vec_size_t = std::vector<unsigned>::size_type;
private:
    Q_OBJECT

    // member variables
    QCustomPlot *m_plot_window = new QCustomPlot;
    int m_width = 600;
    int m_height = 400;
    QVector<QPen> m_color_map;

    // member functions
    void createSubPlots(int num_rows);
    QCPLegend* createLegend(QCPAxisRect *axis_rect);
    void handleData(const FigureData &plot_data, const vec_size_t &num_plots);
    void handleLabels(const FigureData &plot_data);
    void handleTitle(const FigureData &plot_data);
    void padAxes(QCPAxisRect *axis_rect);

public:
    explicit Figure(QWidget *parent = 0);
    ~Figure();
//    PlotData m_plot_data;

    void setSize(int w, int h);
    void plotData(const FigureData &plot_data);
    void showGrid(bool show_grid);
    bool saveImage(QString file_name);

};

#endif // FIGURE_H
