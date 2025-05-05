#include <QPainter>
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <math.h>
#include <string>
#include "chebyshov_approximation.hpp"
#include "spline_approximation.hpp"
#include "window.h"

#define EPSILON 1e-15
#define DEFAULT_A -10
#define DEFAULT_B 10
#define DEFAULT_N 10
#define DEFAULT_K 0
#define L2G(X,Y) (l2g((X), (Y), min_y, max_y))

static
double f_0(double /*x*/)
{
    return 1;
}

static
double f_1(double x)
{
    return x;
}

static
double f_2(double x)
{
    return x * x;
}

static
double f_3(double x)
{
    return x * x * x;
}

static
double f_4(double x)
{
    return x * x * x * x;
}

static
double f_5(double x)
{
    return exp(x);
}

static
double f_6(double x)
{
    return 1 / (25 * x * x + 1);
}

Window::Window(QWidget *parent) : QWidget(parent)
{
    a = DEFAULT_A;
    b = DEFAULT_B;
    n = DEFAULT_N;
    func_id = DEFAULT_K;
    first = true;
    currentApproximation = SPLINE;
    change_func();
}

void Window::calculateChebyshevApproximation(double max_y, double min_y)
{
    double xm;
    double norm = (fabs(min_y) < fabs(max_y) ? fabs(max_y) : fabs(min_y));
    if (approxData.chebyshevCalculated) return;
    if (n <= 50)
    {
        if (!approxData.alpha) approxData.allocate(n);

        for (int m = 0; m < n; m++)
        {
            xm = ((a + b) + (b - a) * cos(PI * (2 * m + 1) / (2 * n))) / 2;
            approxData.F[m] = ((m != n / 2) ? f(xm) : f(xm) + p * 0.1 * norm);
        }

        ChebyshovAproximation(n, approxData.F, approxData.alpha, approxData.g, approxData.g2, approxData.x);
        approxData.chebyshevCalculated = true;
    }
}

void Window::calculateSplineApproximation(double max_y, double min_y)
{
    if (approxData.splineCalculated) return;
    if (!approxData.c) approxData.allocate(n);

    double step = (b - a) / (n - 1);
    for (int m = 0; m < n; m++)
    {
        approxData.z[m] = a + m * step;
        approxData.F[m] = (m != n / 2) ? f(approxData.z[m]) : f(approxData.z[m]) + p * 0.1 * (fabs(min_y) < fabs(max_y) ? fabs(max_y) : fabs(min_y));
    }

    CalculateDiferences(approxData.dF, approxData.z, approxData.F, n);
    double left_derivative = derivative(a, func_id);
    double right_derivative = derivative(b, func_id);
    CalculateParametrs(approxData.A, approxData.g, approxData.dF, approxData.z, left_derivative, right_derivative, n);
    CalculateCoeficients(approxData.c, approxData.z, approxData.g, approxData.dF, approxData.F, n);
    approxData.splineCalculated = true;
}

void Window::clearApproximationData()
{
    approxData.clear();
}

QSize Window::minimumSizeHint() const
{
    return QSize(100, 100);
}

QSize Window::sizeHint() const
{
    return QSize(1000, 1000);
}

int Window::parse_command_line(int argc, char *argv[])
{
    if (argc == 1)
        return 0;

    if (argc == 2)
        return -1;

    if (argc != 5 || sscanf(argv[1], "%lf", &a) != 1 || sscanf(argv[2], "%lf", &b) != 1 ||
        sscanf(argv[3], "%d", &n) != 1 || sscanf(argv[4], "%d", &func_id) != 1)
    {
        printf("Usage ./basic_graph a b n func_id");
        return -2;
    }
    return 0;
}

void Window::change_func()
{
    if (!first)
    {
        func_id = (func_id + 1) % 7;
    }
    else
    {
        first = false;
    }
    clearApproximationData();
    update_function();
    update();
}

void Window::update_function()
{
    switch (func_id)
    {
        case 0:
            f_name = "f (x) = 1";
            f = f_0;
            break;
        case 1:
            f_name = "f (x) = x";
            f = f_1;
            break;
        case 2:
            f_name = "f (x) = x*x";
            f = f_2;
            break;
        case 3:
            f_name = "f (x) = x * x * x";
            f = f_3;
            break;
        case 4:
            f_name = "f (x) = x*x*x*x";
            f = f_4;
            break;
        case 5:
            f_name = "f (x) = e^x";
            f = f_5;
            break;
        case 6:
            f_name = "f (x) = 1/(25*x*x+1)";
            f = f_6;
            break;
    }
}

void Window::toggle_approximation()
{
    if (currentApproximation == SPLINE)
        currentApproximation = CHEBYSHEV;
    else if (currentApproximation == CHEBYSHEV)
        currentApproximation = BOTH;
    else if (currentApproximation == BOTH)
        currentApproximation = ERRORS;
    else
        currentApproximation = SPLINE;

    update();
}

void Window::zoom_in()
{
    double center = (a + b) / 2.0;
    double half_length = (b - a) / 4.0;
    a = center - half_length;
    b = center + half_length;
    clearApproximationData();
    update();
}

void Window::zoom_out()
{
    double center = (a + b) / 2.0;
    double half_length = (b - a);
    a = center - half_length;
    b = center + half_length;
    clearApproximationData();
    update();
}

void Window::increase_points()
{
    n *= 2;
    clearApproximationData();
    update();
}

void Window::decrease_points()
{
    if (n > 3)
    {
        n /= 2;
        clearApproximationData();
        update();
    }
}

void Window::point_down()
{
    p--;
    clearApproximationData();
    update();
}

void Window::point_up()
{
    p++;
    clearApproximationData();
    update();
}

void Window::keyPressEvent(QKeyEvent *event)
{
    if (event->key() == Qt::Key_0)
        change_func();
    else if (event->key() == Qt::Key_1)
        toggle_approximation();
    else if (event->key() == Qt::Key_2)
        zoom_in();
    else if (event->key() == Qt::Key_3)
        zoom_out();
    else if (event->key() == Qt::Key_4)
        increase_points();
    else if (event->key() == Qt::Key_5)
        decrease_points();
    else if (event->key() == Qt::Key_6)
        point_up();
    else if (event->key() == Qt::Key_7)
        point_down();
}

QPointF Window::l2g(double x_loc, double y_loc, double y_min, double y_max)
{
    double x_gl = (x_loc - a) / (b - a) * width();
    if (fabs(y_max - y_min) < 1e-15)
        y_max+=0.0000000000000005;
    double y_gl = (y_max - y_loc) / (y_max - y_min) * height();
    return QPointF(x_gl, y_gl);
}

void drawSmiley(QPainter &painter, int x, int y, int size)
{
    // Рисуем лицо
    painter.drawEllipse(x, y, size, size);
    /*
    // Рисуем глаза
    painter.setBrush(Qt::black);
    painter.drawEllipse(x + size / 4, y + size / 4, size / 8, size / 8); // Левый глаз
    painter.drawEllipse(x + 3 * size / 4 - size / 8, y + size / 4, size / 8, size / 8); // Правый глаз
       
    // Рисуем улыбку
    QPainterPath smile;
    smile.moveTo(x + size / 4, y + size / 2);
    smile.quadTo(x + size / 2, y + 3 * size / 4, x + 3 * size / 4, y + size / 2);
    painter.drawPath(smile);
    */
}

void Window::paintEvent(QPaintEvent *)
{
    QPainter painter(this);
    int M = 1080;
    update_function();
    
    double x1, x2, y1, y2;
    double max_y, min_y;
    double delta_y, delta_x = (b - a) / M;
    double max;
    const char* prefix = "max{|Fmax||,|Fmin|} = ";
    char* strochka = new char[strlen(prefix) + 20];
    QPen pen_black(Qt::black, 2, Qt::SolidLine);
    QPen pen_red(Qt::red, 0, Qt::SolidLine);
    QPen pen_blue(Qt::blue, 3, Qt::SolidLine);
    QPen pen_green(Qt::green, 3, Qt::SolidLine);
    QPen pen_grid(Qt::lightGray, 0, Qt::DotLine);
    
    max_y = min_y = f(a);
    for (x1 = a; x1 <= b; x1 += delta_x)
    {
        y1 = f(x1);
        if (y1 < min_y) min_y = y1;
        if (y1 > max_y) max_y = y1;
    }
    
    delta_y = 0.1 * (max_y - min_y);
    min_y -= delta_y;
    max_y += delta_y;

    painter.setPen(pen_grid);
    double grid_step_x = (b - a) / 20.0;
    double grid_step_y = (max_y - min_y) / 20.0;
    
    for (double x = a; x <= b; x += grid_step_x)
    {
        QPointF start = l2g(x, min_y, min_y, max_y);
        QPointF end = l2g(x, max_y, min_y, max_y);
        painter.drawLine(start, end);
    }
    
    if (fabs(grid_step_y) < 0.0000001)
    {
        for (double x = a; x <= b; x += grid_step_x)
        {
            QPointF start = l2g(x, -min_y, -min_y, max_y);
            QPointF end = l2g(x, max_y, -min_y, max_y);
            painter.drawLine(start, end);
        }
        grid_step_y = (max_y + min_y) / 20.0;
        for (double y = -min_y; y <= max_y; y += grid_step_y)
        {
            QPointF start = l2g(a, y, -min_y, max_y);
            QPointF end = l2g(b, y, -min_y, max_y);
            painter.drawLine(start, end);
        }
    }
    
    for (double y = min_y; y <= max_y; y += grid_step_y)
    {
        QPointF start = l2g(a, y, min_y, max_y);
        QPointF end = l2g(b, y, min_y, max_y);
        painter.drawLine(start, end);
    }

    // Chebyshev approximation
    if ((n <= 50) && (currentApproximation == CHEBYSHEV || currentApproximation == BOTH || currentApproximation == ERRORS))
    {
        calculateChebyshevApproximation(max_y, min_y);
        painter.setPen(pen_blue);
        
        max_y = min_y = 0;
        for (x1 = a; x1 - b < 1.e-16; x1 += delta_x)
        {
            y1 = currentApproximation != ERRORS ? 
                ChebyshovValue(x1, a, b, n, approxData.alpha) :
                fabs(ChebyshovValue(x1, a, b, n, approxData.alpha) - f(x1));
            if (y1 < min_y) min_y = y1;
            if (y1 > max_y) max_y = y1;
        }
        
        delta_y = 0.01 * (max_y - min_y);
        max = (fabs(min_y) < fabs(max_y) ? fabs(max_y) : fabs(min_y));
        prefix = "max{|Fmax||,|Fmin|} = ";
        min_y -= delta_y;
        max_y += delta_y;
        sprintf(strochka, "%s%e", prefix, max);
        printf("%s\n", strochka);
        painter.drawText(0, 110, strochka);

        x1 = a;
        y1 = currentApproximation != ERRORS ? 
            ChebyshovValue(x1, a, b, n, approxData.alpha) :
            fabs(ChebyshovValue(x1, a, b, n, approxData.alpha) - f(x1));
        
        for (x2 = x1 + delta_x; x2 - b < 1.e-16; x2 += delta_x)
        {
            y2 = currentApproximation != ERRORS ? 
                ChebyshovValue(x2, a, b, n, approxData.alpha) :
                fabs(ChebyshovValue(x2, a, b, n, approxData.alpha) - f(x2));
            painter.drawLine(L2G(x1, y1), L2G(x2, y2));
            x1 = x2, y1 = y2;
        }
    }

    // Spline approximation
    if (currentApproximation == SPLINE || currentApproximation == BOTH || currentApproximation == ERRORS)
    {
        calculateSplineApproximation(max_y, min_y);
        painter.setPen(pen_green);
        
        max_y = min_y = 0;
        for (x1 = a; x1 - b < 1.e-16; x1 += delta_x)
        {
            y1 = currentApproximation != ERRORS ?
                SplineValue(x1, a, b, n, approxData.c, approxData.z) :
                fabs(SplineValue(x1, a, b, n, approxData.c, approxData.z) - f(x1));
            if (y1 < min_y) min_y = y1;
            if (y1 > max_y) max_y = y1;
        }
        
        delta_y = 0.01 * (max_y - min_y);
        max = (fabs(min_y) < fabs(max_y) ? fabs(max_y) : fabs(min_y));
        prefix = "max{|Fmax||,|Fmin|} = ";
        min_y -= delta_y;
        max_y += delta_y;
        sprintf(strochka, "%s%e", prefix, max);
        printf("%s\n", strochka);
        painter.drawText(0, 130, strochka);

        x1 = a;
        y1 = currentApproximation != ERRORS ?
            SplineValue(x1, a, b, n, approxData.c, approxData.z) :
            fabs(SplineValue(x1, a, b, n, approxData.c, approxData.z) - f(x1));
        
        for (x2 = x1 + delta_x; x2 - b < 1.e-16; x2 += delta_x)
        {
            y2 = currentApproximation != ERRORS ?
                SplineValue(x2, a, b, n, approxData.c, approxData.z) :
                fabs(SplineValue(x2, a, b, n, approxData.c, approxData.z) - f(x2));
            painter.drawLine(L2G(x1, y1), L2G(x2, y2));
            x1 = x2, y1 = y2;
        }
    }

    painter.setPen(pen_black);
    painter.drawLine(L2G(a, 0), L2G(b, 0));
    painter.drawLine(L2G(0, min_y), L2G(0, max_y));
    
    QFont font2("Arial", 12, QFont::Bold);
    painter.setFont(font2);
    painter.setPen(pen_black);
    painter.drawText(0, 20, f_name);
    
    prefix = "n = ";
    sprintf(strochka, "%s%d", prefix, n);
    painter.drawText(0, 95, strochka);

    prefix = "p = ";
    sprintf(strochka, "%s%d", prefix, p);
    painter.drawText(0, 155, strochka);
    
    if (currentApproximation == CHEBYSHEV)
    {
        painter.setPen(Qt::black);
        painter.drawText(0, 45, "Chebyshev Approximation");
        painter.setBrush(Qt::blue);
        drawSmiley(painter, 210, 30, 20);
    }
    else if (currentApproximation == SPLINE)
    {
        painter.setPen(Qt::black);
        painter.drawText(0, 45, "Spline Approximation");
        painter.setBrush(Qt::green);
        drawSmiley(painter, 170, 30, 20);
    }
    else if (currentApproximation == BOTH)
    {
        painter.setPen(Qt::black);
        painter.drawText(0, 45, "Spline Approximation");
        painter.drawText(0, 70, "Chebyshev Approximation");
        painter.setBrush(Qt::green);
        drawSmiley(painter, 170, 30, 20);
        painter.setBrush(Qt::blue);
        drawSmiley(painter, 210, 55, 20);
    }
    else if (currentApproximation == ERRORS)
    {
        painter.setPen(pen_black);
        painter.drawText(0, 40, "ERRORS");
        painter.setPen(Qt::black);
        painter.drawText(0, 58, "Spline Approximation");
        painter.drawText(0, 80, "Chebyshev Approximation");
        painter.setBrush(Qt::green);
        drawSmiley(painter, 170, 42, 20);
        painter.setBrush(Qt::blue);
        drawSmiley(painter, 210, 62, 20);
    }

    painter.setPen(pen_black);
    QFont font = painter.font();
    font.setPointSize(10);
    painter.setFont(font);

    auto draw_x_mark = [&](double x_value) {
        QPointF mark_pos = l2g(x_value, 0, min_y, max_y);
        painter.drawLine(mark_pos.x(), mark_pos.y() - 5, mark_pos.x(), mark_pos.y() + 5);
        painter.drawText(mark_pos.x() * 0.97 + 5, mark_pos.y() - 20, QString::number(x_value));
    };

    if (!(func_id == 0 && n > 50 && currentApproximation == CHEBYSHEV))
    {
        draw_x_mark(a);
        draw_x_mark(b);

        if (a < 0 && b > 0) draw_x_mark(0);
        if (a < 1 && b > 1) draw_x_mark(1);
        if (a < -1 && b > -1) draw_x_mark(-1);
    }
    delete[] strochka;
}