QT += core gui
QT += widgets printsupport

QMAKE_CXXFLAGS = -std=c++11 -Wall  -Wextra -pedantic

SOURCES += \
    main.cpp \
    smdsim.cpp \
    tools/controls/control_algs.cpp \
    tools/num_methods/numeric_algs.cpp \
    tools/plot_tools/figure.cpp \
    tools/plot_tools/qcustomplot.cpp



HEADERS += \
    smdsim.hpp \
    tools/controls/control_algs.hpp \
    tools/num_methods/numeric_algs.hpp \
    tools/num_methods/integrators.hpp \
    tools/plot_tools/figure.hpp \
    tools/plot_tools/qcustomplot.h

