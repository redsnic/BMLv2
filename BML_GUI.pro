#-------------------------------------------------
#
# Project created by QtCreator 2017-09-25T18:07:06
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = BML_GUI
TEMPLATE = app
QMAKE_CXXFLAGS += -std=c++11
QT += webengine
QT += webenginewidgets
QT += widgets
QT += concurrent

# remove possible other optimization flags
QMAKE_CXXFLAGS  -= -O
QMAKE_CXXFLAGS  -= -O1
QMAKE_CXXFLAGS  -= -O2

# add the desired -O3 if not present
QMAKE_CXXFLAGS  *= -O3

SOURCES += main.cpp\
        mainmenu.cpp \
    Algorithm.cpp \
    GenerateReplicate.cpp \
    InitTree.cpp \
    IOHelper.cpp \
    ParseDataMatrix.cpp \
    SearchDAGs.cpp \
    SearchTrees.cpp \
    StructLearn.cpp \
    Utils.cpp \
    MAFtoMatrixTranslator.cpp \
    helpviewer.cpp \
    maftranslator.cpp \
    StringTokenizer.cpp

HEADERS  += mainmenu.h \
    DataStructures.h \
    GenerateReplicate.h \
    InitTree.h \
    IOHelper.h \
    MAFtoMatrixTranslator.h \
    ParseDataMatrix.h \
    SearchDAGs.h \
    SearchTrees.h \
    StructLearn.h \
    Utils.h \
    helpviewer.h \
    Algorithm.h \
    maftranslator.h \
    StringTokenizer.h

FORMS    += mainmenu.ui \
    helpviewer.ui \
    maftranslator.ui
