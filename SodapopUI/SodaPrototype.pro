QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++11
CONFIG += c++14

# The following define makes your compiler emit warnings if you use
# any Qt feature that has been marked deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    cellUI.cpp \
    geneListUI.cpp \
    geneUI.cpp \
    main.cpp \
    mainwindow.cpp \
    summaryUI.cpp

HEADERS += \
    cellUI.h \
    geneListUI.h \
    geneUI.h \
    mainwindow.h ui_mainwindow.h \
    summaryUI.h \
    ui_cellUI.h \
    ui_geneListUI.h \
    ui_geneUI.h \
    ui_mainwindow.h \
    ui_summaryUI.h

FORMS += \
    cellUI.ui \
    geneListUI.ui \
    geneUI.ui \
    mainwindow.ui \
    summaryUI.ui

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
