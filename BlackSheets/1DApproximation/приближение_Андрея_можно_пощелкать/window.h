#ifndef WINDOW_H
#define WINDOW_H

#include <QtWidgets/QtWidgets>

class Window : public QWidget
{
  Q_OBJECT

private:
  int func_id;
  int view_id = 0;
  int zoom = 0;
  int disturbance = 0;
  const char *f_name;
  double a;
  double b;
  size_t n;
  size_t n_mas_1 = 50;
  size_t n_mas = 40'000'000;
  double (*f) (double);
  double (*dd) (double);
  
  double *x_t01, *y_t01, *alpha, *z;
  double *x_t02, *y_t02, *mas_4n;
  
  int print_graph_1 = 0;
  int print_graph_2 = 0;
  int print_graph_diff = 1;

public:
  Window (QWidget *parent);

  QSize minimumSizeHint () const;
  QSize sizeHint () const;
  
  double max_f_ab (int func_id, double (*f) (double));
  int parse_command_line (int argc, char *argv[]);
  int enough_memory ();
  //int reenough_memory (int n_new);
  int delete_memory ();
  void set_func ();
  void init_memory ();
  void recalculate ();
  void set_view ();
  QPointF l2g (double x_loc, double y_loc, double y_min, double y_max);
public slots:
  void change_func ();
  void change_view ();
  void change_zoom_up ();
  void change_zoom_down ();
  void change_n_down ();
  void change_n_up ();
  void change_disturbance_down ();
  void change_disturbance_up ();

protected:
  void paintEvent (QPaintEvent *event);
};

#endif
