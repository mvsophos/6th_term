#ifndef WINDOW_HPP
#define WINDOW_HPP

#include <QtWidgets/QtWidgets>
#include <cmath>
//#include "cheba.hpp"

//#include "solver.hpp"

class Window : public QWidget {
    Q_OBJECT

private:
    int func_id = 0;
    int mode = 0;
    int pamyat = 40000000;
    const char *f_name = nullptr;
    int zoom = 0;
    double a0;
    double b0;
    double h;
    int /* size_t */ n;
    int /* size_t */ p = 0;

    int graph_1 = 1;
    int graph_2 = 0;
    int graph_res = 0;

    // x_1 - это массив точек в которых приближается функция (для Чебышева), x_2 - массив точек в которых приближается функция (для второго метода)
    double *x_1, *y_1;
    double *x_2, *y_2;
    double *alpha, *z;     // mas_4n = c для 2го приближения
    double *ksi, *c, *v, *a1, *c1, *d1;     // для второго приближения

    double (*f)(double);
    double (*dd)(double);

public:
    Window (QWidget *parent);
    QSize minimumSizeHint() const;
    QSize sizeHint() const;
    int parse_command_line(int argc, char *argv[]);
    QPointF l2g(double x_loc, double y_loc, double y_min, double y_max);

    void set_func();
    void set_mode();
    void recalculate();
    void init_memory();
    int enough_memory();
    int delete_memory();
    double max_func(int func_id, double (*f)(double));

public slots:
    void change_func();
    void change_approx();
    void zoom_in();
    void zoom_out();
    void increase_n();
    void decrease_n();
    void increase_p();
    void decrease_p();

protected:
    void paintEvent(QPaintEvent *event);
};



#endif