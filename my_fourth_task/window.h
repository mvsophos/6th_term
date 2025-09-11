#ifndef WIND_HPP
#define WIND_HPP

#include <QtWidgets/QtWidgets>
#include "thread.h"

#define DEF_ABS_MIN 1e-50

class Window : public QWidget {
    Q_OBJECT
private:
    QWidget *widget;
    pthread_cond_t cond = PTHREAD_COND_INITIALIZER;
    pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
    int func_id;
    const char *f_name;
    const char *program_name;
    int mode = 0;
    double f_max;
    double F_max, PF_max, RES_max;
    int parameter = 0;
    double norma = 0;

    Args *args = nullptr;
    pthread_t *tid = nullptr;
    Dataset data;
    const int task = 8;
    bool threads_created = false;

public:
    Window(QWidget *parent);
    ~Window();

    QSize minimumSizeHint() const;
    QSize sizeHint() const;
    QPointF l2g(double x_loc, double y_loc, double y_min, double y_max);

    int parse_command_line(int argc, char *argv[]);
    void set_f();
    void draw_one_element(QPainter *painter, double x1, double y1, double x2, double y2, double x3, double y3, QColor color);
    void draw_f(QPainter *painter);
    void draw_Pf(QPainter *painter);
    void draw_residual(QPainter *painter);
    void draw_txt(QPainter *painter);
    bool msr_ready();
    void msr_wait();

public slots:
    void change_func();
    void change_mode();
    void zoom_in();
    void zoom_out();
    void twice_n();
    void halve_n();
    void increase_p();
    void decrease_p();
    void twice_m();
    void halve_m();
    void close();

protected:
    void paintEvent(QPaintEvent *event);
};



#endif
