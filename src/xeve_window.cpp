#include "xeve_window.h"
#include "ui_xeve_window.h"

XeVe_Window::XeVe_Window(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::XeVe_Window)
{
    ui->setupUi(this);
}

XeVe_Window::~XeVe_Window()
{
    delete ui;
}
