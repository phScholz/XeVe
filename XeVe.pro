#-------------------------------------------------
#
# Project created by QtCreator 2013-10-01T15:05:24
#
#-------------------------------------------------

QT       += core gui

TARGET = bin/XeVe
TEMPLATE = app


SOURCES += src/main.cpp\
        src/xeve_window.cpp \
    src/XiaEventReader.cpp \
    src/polynom.cpp \
    src/addback.cpp \
    src/advanced_strings.cpp \
    src/evaluate.cpp

HEADERS  += src/xeve_window.h \
    src/XiaEventReader.h \
   src/polynom.h \
    src/event.h \
    src/addback.h \
    src/dgf_format.h \
    src/advanced_strings.h \
	src/evaluate.h

FORMS    += src/xeve_window.ui

OBJECTS_DIR = ./obj
MOC_DIR     = ./moc
UI_DIR = ./ui

SYSTEM +=helios

contains( SYSTEM, helios ) {
    message(System hostname has been defined as helios ...)
    INCLUDEPATH += -I/ikpv2/include -mcmodel=large
    LIBS += -L/ikpv2/lib -lmfile -lm -mcmodel=large -lpthread
} else {
    INCLUDEPATH +=-I/usr/include -mcmodel=large
    LIBS+=   -L/usr/lib -lmfile -lm -mcmodel=large -lpthread
}
