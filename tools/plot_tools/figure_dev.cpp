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

#include <QApplication>
#include <QtWidgets>
#include <cmath>

#include "figure.hpp"

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);
    QVector<double> a;
    QVector<double> b;
    QVector<double> c;
    QVector<double> t_vec;

    int n = 1000;
    for (int i = 0; i <= n; ++i)
    {
        double t = i*0.01;
        t_vec.push_back(t);
        a.push_back(t*t);
        b.push_back(t*t+5);
        c.push_back(sin(t));
    }

    QVector<unsigned> plot_ind1{0, 1};
    QVector<unsigned> plot_ind2{2};


    PlotData plot_data;
    plot_data.t_vec.push_back(t_vec);

    plot_data.plot_data.push_back(a);
    plot_data.plot_data.push_back(b);
    plot_data.plot_data.push_back(c);

    plot_data.plot_inds.push_back(plot_ind1);
    plot_data.plot_inds.push_back(plot_ind2);

    plot_data.legend_data.push_back("");
    plot_data.legend_data.push_back("b");
//    plot_data.legend_data.push_back("c");

    plot_data.x_labels.push_back("Time (s)");
    plot_data.x_labels.push_back("Time (s)");

    plot_data.y_labels.push_back("a");
    plot_data.y_labels.push_back("b");

    plot_data.title = "Plot Test";

    QString filename = "/home/james/Documents/AnimationProject/doc/Images/test.png";


    Figure fig;
    fig.plotData(plot_data);
    fig.show();
    fig.saveImage(filename);
    return app.exec();
}


