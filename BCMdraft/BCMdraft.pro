TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
QMAKE_CXXFLAGS_RELEASE += -std=c++11 -lpthread
QMAKE_CXXFLAGS_RELEASE += -fopenmp
QMAKE_LFLAGS_RELEASE += -fopenmp

QMAKE_CXXFLAGS_RELEASE += -L/opt/mpc-1.0.1/lib


QMAKE_CXXFLAGS_DEBUG += -fopenmp
QMAKE_LFLAGS_DEBUG += -fopenmp

QMAKE_CXXFLAGS_DEBUG += -std=c++11 -lpthread


QMAKE_CXXFLAGS_RELEASE -=-O2

QMAKE_CXXFLAGS_RELEASE +=-O3



SOURCES +=Source/VariablesNormalDistribution.cpp\
    Source/Variables.cpp\
    Source/Experiment.cpp \
    Source/ABC_BCM.cpp \
    Source/PosteriorLevenbergMarquardt.cpp\
    Source/MCMC.cpp \
    Source/BCM_CD137.cpp \
    Source/MatrixInverse.cpp \
    main.cpp \
    Source/matrixCholesky.cpp
HEADERS += \
    Include/ABC_BCM.h \
    Include/MatrixInverse.h \
    Include/PosteriorLevenbergMarquardt.h \
    Include/BCM_CD137.h \
    Include/Experiment.h \
    Include/MCMC.h \
    Include/Variables.h \
    Include/VariablesNormalDistribution.h

win32{
LIBS += -L$$PWD/bin -lcygblas \
        -L$$PWD/bin -lcyglapack
} else {
LIBS += -L$$PWD/bin -lblas \
        -L$$PWD/bin -llapack
}
