#include <QPainter>
#include <stdio.h>
#include <math.h>

#include "window.h"
#include "f.h"
#include "solve.h"

#define MONITOR_SIZE_W 1920
#define DEFAULT_A -10
#define DEFAULT_B 10
#define DEFAULT_N 10
#define L2G(X,Y) (l2g ((X), (Y), min_y, max_y))

Window::Window (QWidget *parent)
  : QWidget (parent)
{
  //printf("WINDOW_CONSTRUCTOR\n");
  a = DEFAULT_A;
  b = DEFAULT_B;
  n = DEFAULT_N;

  func_id = 0;
}

QSize Window::minimumSizeHint () const
{
  return QSize (100, 100);
}

QSize Window::sizeHint () const
{
  return QSize (1000, 1000);
}

int Window::parse_command_line (int argc, char *argv[])
{
  if (argc != 5 || (sscanf (argv[1], "%lf", &a) != 1)
                || (sscanf (argv[2], "%lf", &b) != 1)
                || b - a < 1.e-6
                || (sscanf (argv[3], "%lu", &n) != 1)
                || (sscanf (argv[4], "%d", &func_id) != 1) || (func_id < 0) || (func_id > 6)
                || n <= 0)
    return -1;
  if (n > n_mas)
    n = n_mas;
  if (n < 2)
    n = 2;
  return 0;
}

int Window::enough_memory ()
{
  // Выделение памяти
  if (   (x_t01 = new double[n_mas]) == nullptr
      || (y_t01 = new double[n_mas]) == nullptr
      || (alpha = new double[n_mas]) == nullptr
      || (z = new double[n_mas]) == nullptr
      || (x_t02 = new double[n_mas]) == nullptr
      || (y_t02 = new double[n_mas]) == nullptr
      || (mas_4n = new double[4*n_mas]) == nullptr) 
    {
      if (x_t01 != nullptr) delete[] x_t01;
      if (y_t01 != nullptr) delete[] y_t01;
      if (alpha != nullptr) delete[] alpha;
      if (z != nullptr) delete[] z;
      if (x_t02 != nullptr) delete[] x_t02;
      if (y_t02 != nullptr) delete[] y_t02;
      if (mas_4n != nullptr) delete[] mas_4n;
      QMessageBox::warning (0, "Error", 
                            "Not enough memory!");
      return -2;
    }
  return 0;
}

/*int Window::reenough_memory (int new_n)
{
  
  double *x_t01_new = nullptr;
  double *y_t01_new = nullptr;
  double *alpha_new = nullptr;
  double *z_new = nullptr;
  double *x_t02_new = nullptr;
  double *y_t02_new = nullptr;
  double *mas_4n_new = nullptr;
  
  // Выделение памяти
  if (   (x_t01_new = new double[new_n]) == nullptr
      || (y_t01_new = new double[new_n]) == nullptr
      || (alpha_new = new double[new_n]) == nullptr
      || (z_new = new double[n_mas]) == nullptr
      || (x_t02_new = new double[new_n]) == nullptr
      || (y_t02_new = new double[new_n]) == nullptr
      || (mas_4n_new = new double[4*new_n]) == nullptr) 
    {
      if (x_t01_new != nullptr) delete[] x_t01_new;
      if (y_t01_new != nullptr) delete[] y_t01_new;
      if (alpha_new != nullptr) delete[] alpha_new;
      if (z_new != nullptr) delete[] z_new;
      if (x_t02_new != nullptr) delete[] x_t02_new;
      if (y_t02_new != nullptr) delete[] y_t02_new;
      if (mas_4n_new != nullptr) delete[] mas_4n_new;
      QMessageBox::warning (0, "Warning", 
                            "Not enough memory!");
  printf("Not enough memory!\n");
      return -2;
    }
  //delete_memory();
  
  x_t01 = x_t01_new;
  y_t01 = y_t01_new;
  alpha = alpha_new;
  z = z_new;
  x_t02 = x_t02_new;
  y_t02 = y_t02_new;
  mas_4n = mas_4n_new;
  n = n_mas = new_n;
  printf("n %d n_mas %d\n", n, n_mas);
  return 0;
}*/

void Window::init_memory ()
{
  //printf("INIT_MEM\n");
  if (n <= n_mas_1)
    {
      for (size_t i = 0; i < n; i++)
        {
          x_t01[i] = (a + b) / 2 + (b - a) / 2 * cos_my (M_PI * (2 * i + 1) * 0.5 / n);
          z[i] = 2 * cos_my (M_PI * (2 * i + 1) * 0.5 / n);
        }
    }
  
  double delta = (b - a) / (n - 1);
  for (size_t i = 0; i < n; i++)
    {
      x_t02[i] = a + delta * i;
    }
  
  /*if (n <= n_mas_1)
    {
      printf("INIT x_t01:");
      for (size_t i = 0; i < n; i++)
        printf(" %10.3e", x_t01[i]);
      printf("\n");
      printf("INIT Z:");
      for (size_t i = 0; i < n; i++)
        printf(" %10.3e", z[i]);
      printf("\n");
    }
  printf("INIT x_t02:");
  for (size_t i = 0; i < n; i++)
    printf(" %10.3e", x_t02[i]);
  printf("\n");*/
  
  recalculate ();
  
  //solve_01 (n, x_t01, y_t01, alpha, z);
  //solve_02 (n, x_t02, y_t02, mas_4n, dd(x_t02[0]), dd(x_t02[n-1]));
}

int Window::delete_memory ()
{
  if (x_t01 != nullptr) delete[] x_t01;
  if (y_t01 != nullptr) delete[] y_t01;
  if (alpha != nullptr) delete[] alpha;
  if (z != nullptr) delete[] z;
  if (x_t02 != nullptr) delete[] x_t02;
  if (y_t02 != nullptr) delete[] y_t02;
  if (mas_4n != nullptr) delete[] mas_4n;
  return 0;
}

void Window::set_func ()
{
  switch (func_id)
    {
      case 0:
        f_name = "k = 0; f (x) = 1";
        f = f_0;
        dd = dd_0;
        break;
      case 1:
        f_name = "k = 1; f (x) = x";
        f = f_1;
        dd = dd_1;
        break;
      case 2:
        f_name = "k = 2; f (x) = x ^ 2";
        f = f_2;
        dd = dd_2;
        break;
      case 3:
        f_name = "k = 3; f (x) = x ^ 3";
        f = f_3;
        dd = dd_3;
        break;
      case 4:
        f_name = "k = 4; f (x) = x ^ 4";
        f = f_4;
        dd = dd_4;
        break;
      case 5:
        f_name = "k = 5; f (x) = exp (x)";
        f = f_5;
        dd = dd_5;
        break;
      case 6:
        f_name = "k = 6; f (x) = 1 / (25 * x * x + 1)";
        f = f_6;
        dd = dd_6;
        break;
    }
}

/// change current function for drawing
void Window::change_func ()
{
  func_id = (func_id + 1) % 7;
  set_func ();
  recalculate ();
  update ();
}

void Window::change_n_up ()
{
  //printf("CHANGE_N_UP: %ld %ld\n", n, n_mas);
  if (n == n_mas)
    return;
  
  if (2 * n > n_mas)
    n = n_mas;
  else
    n *= 2;
  init_memory ();
  update ();
}

void Window::change_n_down ()
{
  if (n >= 4)
    {
      n /= 2;
      init_memory ();
      update ();
    }
}

void Window::change_view ()
{
  view_id = (view_id + 1) % 4;
  set_view ();
  update ();
}
void Window::change_zoom_up ()
{
  if (zoom < (int) (8 * sizeof(int)-1))
    {
      zoom += 1;
      update ();
    }
}
void Window::change_zoom_down ()
{
  if (zoom > 0)
    {
      zoom -= 1;
      update ();
    }
}
void Window::change_disturbance_up ()
{
  disturbance += 1;
  recalculate ();
  update ();
}
void Window::change_disturbance_down ()
{
  disturbance -= 1;
  recalculate ();
  update ();
}

void Window::set_view ()
{
  switch (view_id)
    {
      case 0:
        print_graph_1 = 1;
        print_graph_2 = 0;
        print_graph_diff = 0;
        break;
      case 1:
        print_graph_1 = 0;
        print_graph_2 = 1;
        print_graph_diff = 0;
        break;
      case 2:
        print_graph_1 = 1;
        print_graph_2 = 1;
        print_graph_diff = 0;
        break;
      case 3:
        print_graph_1 = 0;
        print_graph_2 = 0;
        print_graph_diff = 1;
        break;
    }
}


double Window::max_f_ab (int func_id, double (*f) (double))
{
  double p1, p2 = 0;
  switch (func_id)
    {
      case 0:
        return 1.;
        break;
      case 1:
        return f (b);
        break;
      case 2:
        p1 = f (a);
        p2 = f (b);
        return (p1 > p2 ? p1 : p2);
        break;
      case 3:
        return f (b);
        break;
      case 4:
        p1 = f (a);
        p2 = f (b);
        return (p1 > p2 ? p1 : p2);
        break;
      case 5:
        return f (b);
        break;
      case 6:
        if (a * b >= 0)
          {
            if (b > 0)
              return f (a);
            else
              return f (b);
          }
        else
          {
            return f (0);
          }
        break;
    }
  return f (0);
}

void Window::recalculate ()
{
  //printf("RECALC: n = %lu\n", n);
  double max_f = 0;
  if (n <= n_mas_1)
    {
      for (size_t i = 0; i < n; i++)
        {
          
          if (disturbance != 0 && i == n/2)
            {
              max_f = max_f_ab (func_id, f);
              y_t01[i] = f(x_t01[i]) + disturbance * 0.01 * max_f;
            }
          else
            {
              y_t01[i] = f(x_t01[i]);
            }
        }
      /*printf("MUST y_t01 (%d):", func_id);
      for (size_t i = 0; i < n; i++)
        printf(" %10.3e", y_t01[i]);
      printf("\n");*/
      solve_01 (n, x_t01, y_t01, alpha, z);
    }
    
  for (size_t i = 0; i < n; i++)
    {
      if (disturbance != 0 && i == n/2)
        {
          max_f = max_f_ab (func_id, f);
          y_t02[i] = f(x_t02[i]) + disturbance * 0.01 * max_f;
        }
      else
        {
          y_t02[i] = f(x_t02[i]);
        }
    }
  /*printf("MUST y_t02 (%d):", func_id);
  for (size_t i = 0; i < n; i++)
    printf(" %10.3e", y_t02[i]);
  printf("\n");*/
  
  solve_02 (n, x_t02, y_t02, mas_4n, dd(x_t02[0]), dd(x_t02[n-1]));
  
  /*printf("MAS_4N: (fid = %d):", func_id);
  for (size_t i = 0; i < n-1; i++)
    {
      printf(" %d %10.3e %10.3e %10.3e %10.3e\n", i, mas_4n[4*i], mas_4n[4*i+1], mas_4n[4*i+2], mas_4n[4*i+3]);
    }*/
}
  

QPointF Window::l2g (double x_loc, double y_loc, double y_min, double y_max)
{
  double a1 = 0.5 * (a + b - (b - a) / (1ll << zoom));
  double b1 = 0.5 * (a + b + (b - a) / (1ll << zoom));
  double x_gl = (x_loc - a1) / (b1 - a1) * width ();
  double y_gl = (y_max - y_loc) / (y_max - y_min) * height ();
  return QPointF (x_gl, y_gl);
}

/// render graph
void Window::paintEvent (QPaintEvent * /* event */)
{  
  QPainter painter (this);
  double x1, x2, y1, y2;
  double max_y, min_y;
  double delta_y, delta_x = (b - a) / (1ll << zoom) / width (); // MONITOR_SIZE_W
  double a1 = 0.5 * (a + b - (b - a) / (1ll << zoom));
  double b1 = 0.5 * (a + b + (b - a) / (1ll << zoom));
  QPen pen_black(Qt::black, 0, Qt::SolidLine);
  QPen pen_red(Qt::red, 0, Qt::SolidLine);
  QPen pen_blue(Qt::blue, 0, Qt::SolidLine);

  painter.setPen (pen_black);
  
  char buf[100];
  
  int graph_1_is_calculated = (n <= n_mas_1);
  
  // calculate min and max for current function
  max_y = min_y = 0;
  for (x1 = a1; x1 - b1 < 1.e-6; x1 += delta_x)
    {
      y1 = f (x1);
      if (y1 < min_y)
        {
          min_y = y1;
          //printf("fu = %10.3e, min_y1 = %10.3e, x1 = %10.3e\n", y1, min_y, x1);
        }
      if (y1 > max_y)
        {
        max_y = y1;
      //printf("fu = %10.3e, max_y1 = %10.3e, x1 = %10.3e\n", y1, max_y1, x1);
        }
    }
  printf("+------------------------------------------------+\n");
  sprintf(buf, "Max_of_func = %10.3e", max_y);
  painter.drawText (0, 40, buf);
  printf("Max_of_func = %10.3e\n", max_y);
  if (!(print_graph_1 || print_graph_2))
    max_y = min_y = 0;
  for (x1 = a1; x1 - b1 < 1.e-6; x1 += delta_x)
    {
      if (print_graph_1 && graph_1_is_calculated)
        {
          y1 = Pf_01 (x1, a, b, n, x_t01, alpha);
          if (y1 < min_y)
            {
              min_y = y1;
          //printf("a1 = %10.3e, min_y1 = %10.3e, x1 = %10.3e\n", y1, min_y, x1);
            }
          if (y1 > max_y)
            max_y = y1;
        }
      if (print_graph_2)
        {
          y1 = Pf_02 (x1, a, b, n, x_t02, mas_4n);
          if (y1 < min_y)
            {
              min_y = y1;
            }
          //printf("a2 = %10.3e, min_y1 = %10.3e, x1 = %10.3e\n", y1, min_y, x1);
          if (y1 > max_y)
            max_y = y1;
        }
      if (print_graph_diff)
        {
          if (graph_1_is_calculated)
            {
              y1 = fabs (f (x1) - Pf_01 (x1, a, b, n, x_t01, alpha));
              if (y1 < min_y)
                {
                  min_y = y1;
          //printf("d1 = %10.3e, min_y1 = %10.3e, x1 = %10.3e\n", y1, min_y, x1);
                }
              if (y1 > max_y)
                max_y = y1;
            }
          y1 = fabs (f (x1) - Pf_02 (x1, a, b, n, x_t02, mas_4n));
          if (y1 < min_y)
            {
              min_y = y1;
          //printf("d2 = %10.3e, min_y1 = %10.3e, x1 = %10.3e\n", y1, min_y, x1);
            }
          if (y1 > max_y)
            max_y = y1;
        }
    }
  printf("max = %10.3e\n", (max_y > -min_y ? max_y : -min_y));
  sprintf(buf, "max = %10.3e", (max_y > -min_y ? max_y : -min_y));
  //printf("min_y1 = %10.3e, max_y1 = %10.3e %d\n", min_y, max_y, width ());
  //sprintf(buf, "max = %10.3e", (max_y > -min_y ? max_y : -min_y));
  painter.drawText (0, 60, buf);
  sprintf(buf, "Zoom: s = %d", zoom);
  painter.drawText (0, 80, buf);
  sprintf(buf, "n/n_max: %ld/%ld", n, n_mas);
  painter.drawText (0, 100, buf);
  sprintf(buf, "disturbance: p = %d", disturbance);
  painter.drawText (0, 120, buf);
  
  if (fabs (max_y - min_y) < 1e-20)
    delta_y = 1e-20;
  else
    delta_y = 0.01 * (max_y - min_y);
  min_y -= delta_y;
  max_y += delta_y;

  // draw axis
  painter.setPen (pen_black);
  painter.drawLine (L2G(a, 0), L2G(b, 0));
  if (a <= 0 && 0 <= b)
    painter.drawLine (L2G(0, min_y), L2G(0, max_y));

  
  painter.drawText (width() - 100, 20, "Legend:");
  
  if (print_graph_1 || print_graph_2)
    {
      painter.setPen (pen_black);
      painter.drawText (width() - 100, 40, "Function");
      
      x1 = a1;
      y1 = f (x1);
      for (x2 = x1 + delta_x; x2 - b1 < 1.e-6; x2 += delta_x) 
        {
          y2 = f (x2);
          // local coords are converted to draw coords
          painter.drawLine (L2G(x1, y1), L2G(x2, y2));

          x1 = x2, y1 = y2;
        }
      x2 = b1;
      y2 = f (x2);
      painter.drawLine (L2G(x1, y1), L2G(x2, y2));
    }
  
  if (print_graph_1 && n <= n_mas_1)
    {
      painter.setPen (pen_blue);
      painter.drawText (width() - 100, 60, "Aprox 1");
      
      // draw approximated line for graph
      x1 = a1;
      y1 = Pf_01 (x1, a, b, n, x_t01, alpha);
      for (x2 = x1 + delta_x; x2 - b1 < 1.e-6; x2 += delta_x) 
        {
          y2 = Pf_01 (x2, a, b, n, x_t01, alpha);
          // local coords are converted to draw coords
          painter.drawLine (L2G(x1, y1), L2G(x2, y2));
          //painter.drawPoint (L2G(x2, y2));

          x1 = x2, y1 = y2;
        }
      x2 = b1;
      y2 = Pf_01 (x2, a, b, n, x_t01, alpha);
      painter.drawLine (L2G(x1, y1), L2G(x2, y2));
    }
  
  if (print_graph_2)
    {
      painter.setPen (pen_red);
      if (print_graph_1 && n <= n_mas_1)
        {
          painter.drawText (width () - 100, 80, "Aprox 2");
        }
      else
        {
          painter.drawText (width () - 100, 60, "Aprox 2");
        }

      // draw approximated line for graph
      x1 = a1;
      y1 = Pf_02 (x1, a, b, n, x_t02, mas_4n);
      for (x2 = x1 + delta_x; x2 - b1 < 1.e-6; x2 += delta_x) 
        {
          y2 = Pf_02 (x2, a, b, n, x_t02, mas_4n);
          // local coords are converted to draw coords
          painter.drawLine (L2G(x1, y1), L2G(x2, y2));

          x1 = x2, y1 = y2;
        }
      x2 = b1;
      y2 = Pf_02 (x2, a, b, n, x_t02, mas_4n);
      painter.drawLine (L2G(x1, y1), L2G(x2, y2));
    }
  
  if (print_graph_diff)
    {
      if (n <= n_mas_1)
        {
          painter.setPen (pen_blue);
          painter.drawText (width () - 100, 40, "Diff 1");
          
          x1 = a1;
          y1 = fabs (f (x1) - Pf_01 (x1, a, b, n, x_t01, alpha));
          for (x2 = x1 + delta_x; x2 - b1 < 1.e-6; x2 += delta_x) 
            {
              y2 = fabs (f (x2) - Pf_01 (x2, a, b, n, x_t01, alpha));
              // local coords are converted to draw coords
              painter.drawLine (L2G(x1, y1), L2G(x2, y2));
              //painter.drawPoint (L2G(x2, y2));

              x1 = x2, y1 = y2;
            }
          x2 = b1;
          y2 = fabs (f (x2) - Pf_01 (x2, a, b, n, x_t01, alpha));
          painter.drawLine (L2G(x1, y1), L2G(x2, y2));
        }
      
      painter.setPen (pen_red);
      if (n <= n_mas_1)
        {
          painter.drawText (width () - 100, 60, "Diff 2");
        }
      else
        {
          painter.drawText (width () - 100, 40, "Diff 2");
        }

      // draw approximated line for graph
      x1 = a1;
      y1 = fabs (f (x1) - Pf_02 (x1, a, b, n, x_t02, mas_4n));
      for (x2 = x1 + delta_x; x2 - b1 < 1.e-6; x2 += delta_x) 
        {
          y2 = fabs (f (x2) - Pf_02 (x2, a, b, n, x_t02, mas_4n));
  //sprintf(buf, "(%10.3e, %10.3e)", x2, y2);
  //painter.drawText (L2G(x2, y2), buf);
  //printf("pr_d2: x2 = %10.3e, y2 = %10.3e\n", x2, y2);
          // local coords are converted to draw coords
          painter.drawLine (L2G(x1, y1), L2G(x2, y2));

          x1 = x2, y1 = y2;
        }
      x2 = b1;
      y2 = fabs (f (x2) - Pf_02 (x2, a, b, n, x_t02, mas_4n));
      painter.drawLine (L2G(x1, y1), L2G(x2, y2));
    }
  
  
  /*if (print_graph_1 && n <= n_mas_1)
    {
      painter.setPen (pen_blue);
      for (size_t i = 0; i < n; i++) 
        {
          if (a1 <= x_t01[i] && x_t01[i] <= b1)
            y1 = Pf_01 (x_t01[i], a, b, n, x_t01, alpha);
          // local coords are converted to draw coords
          painter.drawEllipse (L2G(x_t01[i], y1), 3, 3);
        }
    }
  if (print_graph_2)
    {
      painter.setPen (pen_red);
      for (size_t i = 0; i < n; i++) 
        {
          if (a1 <= x_t02[i] && x_t02[i] <= b1)
            y1 = Pf_02 (x_t02[i], a, b, n, x_t02, mas_4n);
          // local coords are converted to draw coords
          painter.drawEllipse (L2G(x_t02[i], y1), 3, 3);
        }
    }*/

  // render function name
  painter.setPen ("blue");
  painter.setPen ("black");
  painter.drawText (0, 20, f_name);

}
