#ifndef XEVE_WINDOW_H
#define XEVE_WINDOW_H

#include <QMainWindow>

namespace Ui {
class XeVe_Window;
}

class XeVe_Window : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit XeVe_Window(QWidget *parent = 0);
    ~XeVe_Window();
    
private:
    Ui::XeVe_Window *ui;
};

#endif // XEVE_WINDOW_H
