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

//#include <algorithm>

#include "figure.hpp"

using std::max_element;
using std::abs;

/**
 * @brief Figure::Figure
 * @param parent
 * Constructor
 */
Figure::Figure(QWidget *parent) : QWidget(parent)
{
    QVBoxLayout *layout = new QVBoxLayout(this);
    layout->addWidget(m_plot_window);
    resize(m_width, m_height);
    setAttribute(Qt::WA_DeleteOnClose, true);

    // Define Default Colors for plot lines
    m_color_map.push_back(QPen(QColor(0, 114, 189)));
    m_color_map.push_back(QPen(QColor(217, 83, 25)));
    m_color_map.push_back(QPen(QColor(237, 177, 32)));
    m_color_map.push_back(QPen(QColor(126, 47, 142)));
    m_color_map.push_back(QPen(QColor(119, 172, 48)));
    m_color_map.push_back(QPen(QColor(77, 190, 238)));
    m_color_map.push_back(QPen(QColor(162, 20, 47)));
}

Figure::~Figure()
{

}

/**
 * @brief Figure::setSize
 * @param w window width in pixels
 * @param h window height in pixels
 * Sets the size of the window. w and h are used to set the member parameters.
 */
void Figure::setSize(int w, int h)
{
    m_width = w;
    m_height = h;
    this->resize(m_width, m_height);
}

/**
 * @brief Figure::plotData
 * @param plot_data:
 */
void Figure::plotData(const FigureData &plot_data)
{
    vec_size_t num_plots = plot_data.m_plot_inds.size();
    createSubPlots(num_plots);
    handleData(plot_data, num_plots);
    handleLabels(plot_data);
    handleTitle(plot_data);
}

/**
 * @brief Figure::showGrid
 * @param show_grid
 */
void Figure::showGrid(bool show_grid)
{
    vec_size_t num_plots = m_plot_window->axisRectCount();
    for (vec_size_t i = 0; i != num_plots; ++i)
    {
        QCPAxisRect *axis_rect = m_plot_window->axisRect(i);
        axis_rect->axis(QCPAxis::atBottom)->grid()->setVisible(show_grid);
        axis_rect->axis(QCPAxis::atLeft)->grid()->setVisible(show_grid);
        foreach(QCPAxis *axis, axis_rect->axes())
        {
            axis->setLayer("axes");
            axis->grid()->setLayer("grid");
        }
    }
}

/**
 * @brief Figure::saveImage
 * @param file_name
 * @return bool
 * Filename should contain the appropriate filename extension, returns status
 */
bool Figure::saveImage(QString file_name)
{
    QPixmap img(grab());
    bool status = img.save(file_name);
    return status;
}

/**
 * @brief Figure::createSubPlots
 * @param num_rows
 */
void Figure::createSubPlots(int num_rows)
{
    m_plot_window->plotLayout()->clear();

    for (int i = 0; i != num_rows; ++i)
    {
        QCPAxisRect *axis_rect =  new QCPAxisRect(m_plot_window, false);
        axis_rect->setupFullAxesBox(true);
        axis_rect->axis(QCPAxis::atTop)->setTicks(false);
        axis_rect->axis(QCPAxis::atRight)->setTicks(false);
        m_plot_window->plotLayout()->addElement(i,0,axis_rect);
    }
    showGrid(true);

}

/**
 * @brief Figure::createLegend
 * @param axis_rect
 * @param plot_inds
 * @param plot_data
 * @return
 */
QCPLegend *Figure::createLegend(QCPAxisRect *axis_rect)
{
    QCPLegend *legend = new QCPLegend;
    axis_rect->insetLayout()->addElement(legend, Qt::AlignRight | Qt::AlignTop);
    legend->setLayer("legend");
    return legend;
}

/**
 * @brief Figure::handleData
 * @param plot_data
 * @param num_plots
 */
void Figure::handleData(const FigureData &plot_data, const vec_size_t &num_plots)
{
    // iterate over subplots
    for(vec_size_t i = 0; i != num_plots; ++i)
    {
        QCPAxisRect *axis_rect{m_plot_window->axisRect(i)};
         QCPLegend *legend = nullptr;

        // create lines
        QVector<vec_size_t> plot_inds = plot_data.m_plot_inds[i];
        vec_size_t num_lines = plot_inds.size();
        for (vec_size_t ii = 0; ii != num_lines; ++ii)
        {
            QVector<double> plot_t_vec;
            vec_size_t plot_ind = plot_inds[ii];

            if (plot_data.m_t_vec.size() == 1)
                plot_t_vec = plot_data.m_t_vec[0];
            else
                plot_t_vec = plot_data.m_t_vec[plot_ind];

            QCPGraph *graph{m_plot_window->addGraph(axis_rect->axis(QCPAxis::atBottom),
                                                    axis_rect->axis(QCPAxis::atLeft))};

            int map_length = m_color_map.size();
            int map_index = ii%map_length;
            graph->setPen(m_color_map[map_index]);
            graph->setData(plot_t_vec, plot_data.m_plot_data[plot_ind]);

            if (plot_ind < plot_data.m_legend_data.size() && !plot_data.m_legend_data[plot_ind].isEmpty())
            {
                if (!legend){
                    legend = createLegend(axis_rect);
                }

                 graph->setName(plot_data.m_legend_data[plot_ind]);
                 legend->addItem(new QCPPlottableLegendItem(legend, graph));
                }
                if (ii == 0){
                    graph->rescaleAxes();
                }

                else
                    graph->rescaleAxes(true);
            }
        padAxes(axis_rect);
        }

}


/**
 * @brief Figure::handleLabels
 * @param plot_data
 * @param num_plots
 * @param axis_rect
 */
void Figure::handleLabels(const FigureData &plot_data)
{
     vec_size_t num_plots = m_plot_window->axisRectCount();

     for(vec_size_t i = 0; i != num_plots; ++i)
     {
        if (i < plot_data.m_x_labels.size() && !plot_data.m_x_labels[i].isEmpty())
            m_plot_window->axisRect(i)->axis(QCPAxis::atBottom)->setLabel(plot_data.m_x_labels[i]);

        if (i < plot_data.m_y_labels.size() && !plot_data.m_y_labels[i].isEmpty())
            m_plot_window->axisRect(i)->axis(QCPAxis::atLeft)->setLabel(plot_data.m_y_labels[i]);

     }
}

/**
 * @brief Figure::handleTitle
 * @param plot_data
 */
void Figure::handleTitle(const FigureData &plot_data)
{
    if(!plot_data.m_title.isEmpty()) {
        m_plot_window->plotLayout()->insertRow(0);
        m_plot_window->plotLayout()->addElement(0, 0, new QCPPlotTitle(m_plot_window, plot_data.m_title));
        setWindowTitle(plot_data.m_title);
    }
}

/**
 * @brief extendAxes
 * @param axis_rect
 * provides extra padding on top and below plot data
 */
void Figure::padAxes(QCPAxisRect *axis_rect)
{
    double upper = axis_rect->axis(QCPAxis::atLeft)->range().upper;
    double lower = axis_rect->axis(QCPAxis::atLeft)->range().lower;
    double scale = 0.1;
    double add = 0;

    if (std::abs(upper) >= std::abs(lower))
        add = scale * std::abs(upper);
    else
        add = scale * std::abs(lower);

    upper = upper + add;
    lower = lower - add;

    axis_rect->axis(QCPAxis::atLeft)->setRange(lower, upper);
}

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
