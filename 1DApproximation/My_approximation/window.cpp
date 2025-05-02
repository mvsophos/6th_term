#include <QPainter>
#include <cstdio>
#include <cmath>

#include "window.hpp"
#include "solver.hpp"
/* #include "cheba.hpp"
#include "spline.hpp" */



#define eps 1e-15
#define fifty 50
#define DEFAULT_A -10
#define DEFAULT_B 10
#define DEFAULT_N 10
#define L2G(X, Y) l2g((X), (Y), min_y, max_y)



Window::Window (QWidget *parent)
    : QWidget(parent) {
    a0 = DEFAULT_A;
    b0 = DEFAULT_B;
    n  = DEFAULT_N;

    func_id = 0;
}

QSize Window::minimumSizeHint() const {
    return QSize(100, 100);
}

QSize Window::sizeHint() const {
    return QSize(1000, 1000);
}

// раньше вместо a0 вводил a
int Window::parse_command_line(int argc, char *argv[]) {
    if (argc != 5   || (sscanf(argv[1], "%lf", &a0) != 1) 
                    || (sscanf(argv[2], "%lf", &b0) != 1)
                    || (sscanf(argv[3], "%d",  &n) != 1)
                    || (sscanf(argv[4], "%d",  &func_id) != 1)
                    || (b0 - a0 < 1.e-6) || (func_id < 0) || (func_id > 6) || (n < 1)) {
        printf("Usage: ./graph  a  b  n  func_id \n");
        return -1;
    }
    //a = a0; b = b0;
    return 0;
}

double f_0(double /* x */) { return 1; }
double f_1(double x) { return x; }
double f_2(double x) { return x * x; }
double f_3(double x) { return x * x * x; }
double f_4(double x) { return x * x * x * x; }
double f_5(double x) { return exp(x); }
double f_6(double x) { return 1 / (25 * x * x  +  1); }

void Window::set_func() {
    switch (func_id) {
    case 0:
        f_name = "k = 0, f(x) = 1";
        f = f_0;
        dd = dd_0;
        break;
    case 1:
        f_name = "k = 1, f(x) = x";
        f = f_1;
        dd = dd_1;
        break;
    case 2:
        f_name = "k = 2, f(x) = x ** 2";
        f = f_2;
        dd = dd_2;
        break;
    case 3:
        f_name = "k = 3, f(x) = x ** 3";
        f = f_3;
        dd = dd_3;
        break;
    case 4:
        f_name = "k = 4, f(x) = x ** 4";
        f = f_4;
        dd = dd_4;
        break;
    case 5:
        f_name = "k = 5, f(x) = exp(x)";
        f = f_5;
        dd = dd_5;
        break;
    case 6:
        f_name = "k = 6, f(x) = 1 / (25 * x ** 2 + 1)";
        f = f_6;
        dd = dd_6;
        break;
    }
}

void Window::set_mode() {
    switch (mode) {
    case 0:
        graph_1   = 1;
        graph_2   = 0;
        graph_res = 0;
        break;
    case 1:
        graph_1   = 0;
        graph_2   = 1;
        graph_res = 0;
        break;
    case 2:
        graph_1   = 1;
        graph_2   = 1;
        graph_res = 0;
        break;
    case 3:
        graph_1   = 0;
        graph_2   = 0;
        graph_res = 1;
        break;
    }
}



double Window::max_func(int func_id, double (*f)(double)) {
    switch (func_id) {
    case 0:
        return 1.;
        break;
    case 1:
        return f(b0);
        break;
    case 2:
        return fmax(f(a0), f(b0));
        break;
    case 3:
        return f(b0);
        break;
    case 4:
        return fmax(f(a0), f(b0));
        break;
    case 5:
        return f(b0);
        break;
    case 6:
        if (a0 * b0 > 0) {
            if (b0 > 0)  return f(a0);
            else        return f(b0);
        }
        else f(0);
    }
    return f(0);        // никогда не выполняется, просто для компиляции
}

// init_memory необходима для формирования массивов из n чисел, если n изменяется

void Window::recalculate() {
    //double max_f = 0;
    if (n <= fifty) {
        for (int i = 0; i < n; i++) {
            if (p != 0 && i == n / 2) {
                //max_f = max_func(func_id, f);
                y_1[i] = f(x_1[i]) /*  + p * 0.01 * max_f */;
            }
            else {
                y_1[i] = f(x_1[i]);
            }
        }
        solve_1(x_1, y_1, n, alpha, z);
    }

    for (int i = 0; i < n; i++) {
        if (p != 0 && i == n / 2) {
            //max_f = max_func(func_id, f);
            y_2[i] = f(x_2[i]) /* + p * 0.01 * max_f */;
        }
        else {
            y_2[i] = f(x_2[i]);
        }
    }

    solve_2(x_2, y_2, mas_4n, n, dd(x_2[0]), dd(x_2[n - 1]));
}

// надо это чуть изменить, x_1 должно быть равномерной
void Window::init_memory() {
    if (n <= fifty) {
        for (int i = 0; i < n; i++) {
            x_1[i] = (a0 + b0) / 2 + ((b0 - a0) / 2) * cos(M_PI * (2 * i + 1) / (2 * n));
            z[i] = 2 * cos(M_PI * (2 * i + 1) / (2 * n));
        }
    }

    double dx = (b0 - a0) / (n - 1);
    for (int i = 0; i < n; i++) {
        x_2[i] = a0 + dx * i;
    }

    recalculate();
}





int Window::enough_memory() {
  // Выделение памяти
  if ( (x_1    = new double[pamyat]) == nullptr
    || (y_1    = new double[pamyat]) == nullptr
    || (alpha  = new double[pamyat]) == nullptr
    || (z      = new double[pamyat]) == nullptr
    || (x_2    = new double[pamyat]) == nullptr
    || (y_2    = new double[pamyat]) == nullptr
    || (mas_4n = new double[4 * pamyat]) == nullptr)
    {
        if (x_1    != nullptr) delete [] x_1;
        if (y_1    != nullptr) delete [] y_1;
        if (alpha  != nullptr) delete [] alpha;
        if (z      != nullptr) delete [] z;
        if (x_2    != nullptr) delete [] x_2;
        if (y_2    != nullptr) delete [] y_2;
        if (mas_4n != nullptr) delete [] mas_4n;
        QMessageBox::warning (0, "Error", "Not enough memory!");
        return -2;
    }
    return 0;
}

int Window::delete_memory() {
    if (x_1 != nullptr) delete[] x_1;
    if (y_1 != nullptr) delete[] y_1;
    if (alpha != nullptr) delete[] alpha;
    if (z != nullptr) delete[] z;
    if (x_2 != nullptr) delete[] x_2;
    if (y_2 != nullptr) delete[] y_2;
    if (mas_4n != nullptr) delete[] mas_4n;
    return 0;
}




void Window::change_func() {
    func_id = (func_id + 1) % 7;
    set_func();
    recalculate();
    update();
}

void Window::change_approx() {
    mode = (mode + 1) % 4;
    switch (mode) {
    case 0:
        graph_1 = 1;
        graph_2 = 0;
        graph_res = 0;
        break;
    case 1:
        graph_1 = 0;
        graph_2 = 1;
        graph_res = 0;
        break;
    case 2:
        graph_1 = 1;
        graph_2 = 1;
        graph_res = 0;
        break;
    case 3:
        graph_1 = 0;
        graph_2 = 0;
        graph_res = 1;
        break;
    }
    update();
}

void Window::zoom_in() {
    /* double center = (a + b) / 2;
    double h = (b - a) / 4;
    a = center - h;
    b = center + h; */
    zoom += 1;
    // чистим память
    update();
}

void Window::zoom_out() {
    /* double center = (a + b) / 2;
    double h = (b - a);
    a = center - h;
    b = center + h; */
    if (zoom > 0) {
        zoom -= 1;
        update();
    }
    // чистим память
}

void Window::increase_n() {         // вероятно нельзя до бесконечности повышать
    if (2 * n >= pamyat) n = pamyat;
    else n *= 2;
    init_memory();
    update();
}

void Window::decrease_n() {
    if (n >= 4) {
        n /= 2;
        init_memory();
        update();
    }
}

void Window::increase_p() {
    p += 1;
    recalculate();
    update();
}

void Window::decrease_p() {
    p -= 1;
    recalculate();
    update();
}

QPointF Window::l2g(double x_loc, double y_loc, double y_min, double y_max) {
    //double a1 = 0.5 * ()
    double a = 0.5 * (a0 + b0 - (b0 - a0) / (1ll << zoom));
    double b = 0.5 * (a0 + b0 + (b0 - a0) / (1ll << zoom));
    //if (y_max - y_min < 1e-5) y_max = 1e-5;
    double x_gl = (x_loc - a) / (b - a) * width();
    double y_gl = (y_max - y_loc) / (y_max - y_min) * height();
    return QPointF(x_gl, y_gl);
}

void Window::paintEvent(QPaintEvent *) {
    QPainter painter(this);
    //int M = 1080;
    set_func();

    //a = a0; b = b0;     // этого не было, а все функции Pf_ были от a и b

    double x1, x2, y1, y2;
    double max_y = f(a0), min_y = f(a0);
    /* double delta_y, delta_x = (b - a) / width(); */

    double delta_y, delta_x = (b0 - a0) / (1ll << zoom) / width (); // MONITOR_SIZE_W
    double a = 0.5 * (a0 + b0 - (b0 - a0) / (1ll << zoom));
    double b = 0.5 * (a0 + b0 + (b0 - a0) / (1ll << zoom));

    QPen pen_black(Qt::black, 2, Qt::SolidLine);
    QPen pen_red(Qt::red, 3, Qt::SolidLine);
    QPen pen_blue(Qt::blue, 3, Qt::SolidLine);

    painter.setPen(pen_black);

    char buf[100];          // для всяких надписей в самом окне

    int graph_is_calculated = (n <= fifty);            // надо чтобы не было слишком много точек, не надо строить

    for (x1 = a; x1 - b < 1e-6; x1 += delta_x) {
        y1 = f(x1);
        if (y1 < min_y) min_y = y1;
        if (y1 > max_y) max_y = y1;
    }

    printf("Maximum = %10.3e;    ", max_y);
    sprintf(buf, "Maximum = %10.3e", max_y);
    painter.drawText(0, 40, buf);


    // находим максимумы функций, чтобы понять где их расположить в окне
    if (!(graph_1 || graph_2)) max_y = min_y = 0;
    for (x1 = a; x1 - b < 1.e-6; x1 += delta_x) {
        if (graph_1 && graph_is_calculated) {
            y1 = Pf_1(x1, a0, b0, n, x_1, alpha);          //
            if (y1 < min_y) min_y = y1;
            if (y1 > max_y) max_y = y1;
        }
        if (graph_2) {
            y1 = Pf_2(x1, a0, b0, n, x_2, mas_4n);          //
            if (y1 < min_y) min_y = y1;
            if (y1 > max_y) max_y = y1;
        }
        if (graph_res) {
            if (graph_is_calculated) {
                y1 = fabs( f(x1) - Pf_1(x1, a0, b0, n, x_1, alpha) );          //
                if (y1 < min_y) min_y = y1;
                if (y1 > max_y) max_y = y1;
            }

            y1 = fabs( f(x1) - Pf_2(x1, a0, b0, n, x_2, mas_4n) );                   //
            if (y1 < min_y) min_y = y1;
            if (y1 > max_y) max_y = y1;
        }
    }

    printf("Max = %10.3e\n",        max_y);              // эта информация будет в консоли
    sprintf(buf, "Max = %10.3e",    max_y);
    painter.drawText(0,  60, buf);
    sprintf(buf, "n = %d", n);
    painter.drawText(0,  80, buf);
    sprintf(buf, "p = %d", p);
    painter.drawText(0, 100, buf);
    sprintf(buf, "[%lf, %lf]", a, b);
    painter.drawText(0, 120, buf);

    printf("\na = %lf, b = %lf\n", a, b);


    if (fabs(max_y - min_y) < 1e-16)    delta_y = 1.e-16;
    else                                delta_y = (max_y - min_y) / 10;

    min_y -= delta_y;       // это чтобы график выводился с запасом сверху и снизу
    max_y += delta_y;

    // рисуем оси
    painter.setPen(pen_black);
    painter.drawLine(L2G(a, 0), L2G(b, 0));
    if (a < 0 && b > 0) painter.drawLine(L2G(0, min_y), L2G(0, max_y));
    
    //painter.drawText( width() - 100, 20, "Legend:" );
    painter.drawText( 10, height() - 100, "Legend:" );




    if (graph_1 || graph_2) {
        painter.setPen(pen_black);
        painter.drawText( 10, height() - 80, "Function" );

        x1 = a;
        y1 = f(a);
        for (x2 = a + delta_x; x2 - b < 1.e-6; x2 += delta_x) {
            y2 = f(x2);
            painter.drawLine(L2G(x1, y1), L2G(x2, y2));
            x1 = x2; y1 = y2;
        }
        x2 = b;
        y2 = f(x2);
        painter.drawLine(L2G(x1, y1), L2G(x2, y2));
    }

    if (graph_1 && n <= fifty) {           // мб нельзя до бесконечности увеличивать
        painter.setPen(pen_blue);
        painter.drawText( 10, height() - 60, "Approximation 1" );

        x1 = a;
        y1 = Pf_1(x1, a0, b0, n, x_1, alpha);
        for (x2 = a + delta_x; x2 - b < 1.e-6; x2 += delta_x) {
            y2 = Pf_1(x2, a0, b0, n, x_1, alpha);
            painter.drawLine(L2G(x1, y1), L2G(x2, y2));
            x1 = x2; y1 = y2;
        }
        x2 = b;
        y2 = Pf_1(x2, a0, b0, n, x_1, alpha);
        painter.drawLine(L2G(x1, y1), L2G(x2, y2));
    }

    if (graph_2) {
        painter.setPen(pen_red);    
        painter.drawText( 10, height() - 40, "Approximation 2" );

        x1 = a;
        y1 = Pf_2(x1, a0, b0, n, x_2, mas_4n);
        for (x2 = a + delta_x; x2 - b < 1.e-6; x2 += delta_x) {
            y2 = Pf_2(x2, a, b, n, x_2, mas_4n);
            painter.drawLine(L2G(x1, y1), L2G(x2, y2));
            x1 = x2; y1 = y2;
        }
        x2 = b;
        y2 = Pf_2(x2, a, b, n, x_2, mas_4n);
        painter.drawLine(L2G(x1, y1), L2G(x2, y2));
    }

    if (graph_res) {
        if (n < fifty) {
            painter.setPen(pen_blue);
            painter.drawText( 10, height() - 80, "Difference 1" );

            x1 = a;
            y1 = fabs(f(x1) - Pf_1(x1, a0, b0, n, x_1, alpha));
            for (x2 = a + delta_x; x2 - b < 1.e-6; x2 += delta_x) {
                y2 = fabs(f(x2) - Pf_1(x2, a0, b0, n, x_1, alpha));
                painter.drawLine(L2G(x1, y1), L2G(x2, y2));
                x1 = x2; y1 = y2;
            }
            x2 = b;
            y2 = fabs(f(x2) - Pf_1(x2, a0, b0, n, x_1, alpha));
            painter.drawLine(L2G(x1, y1), L2G(x2, y2));
        }

        
        painter.setPen(pen_red);
        painter.drawText( 10, height() - 60, "Difference 2" );

        x1 = a;
        y1 = fabs(f(x1) - Pf_2(x1, a0, b0, n, x_2, mas_4n));
        for (x2 = a + delta_x; x2 - b < 1.e-6; x2 += delta_x) {
            y2 = fabs(f(x2) - Pf_2(x2, a0, b0, n, x_2, mas_4n));
            painter.drawLine(L2G(x1, y1), L2G(x2, y2));
            x1 = x2; y1 = y2;
        }
        x2 = b;
        y2 = fabs(f(x2) - Pf_2(x2, a0, b0, n, x_2, mas_4n));
        painter.drawLine(L2G(x1, y1), L2G(x2, y2));
    }

    // выводим название функции
    painter.setPen("blue");
    painter.setPen("black");
    painter.drawText(0, 20, f_name);
}

