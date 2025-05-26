
#include <QtWidgets/QApplication>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QAction>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QMessageBox>

#include "window.h"

int main (int argc, char *argv[])
{
  QApplication app (argc, argv);

  QMainWindow *window = new QMainWindow;
  QMenuBar *tool_bar = new QMenuBar (window);
  Window *graph_area = new Window (window);
  QAction *action;

  //printf("!!!\n");
  if (graph_area->parse_command_line (argc, argv))
    {
      QMessageBox::warning (0, "Wrong input arguments!", 
                            "Usage ./<a.out> a b n k");
      return -1;
    }
  if (graph_area->enough_memory ())
    {
      QMessageBox::warning (0, "Error", 
                            "Not enough memory!");
      return -2;
    }
  graph_area->set_func ();
  graph_area->set_view ();
  graph_area->init_memory ();
  

  action = tool_bar->addAction ("Change function (0)", graph_area, SLOT (change_func ()));
  action->setShortcut (QString ("0"));
  
  action = tool_bar->addAction ("Change view (1)", graph_area, SLOT (change_view ()));
  action->setShortcut (QString ("1"));
  
  action = tool_bar->addAction ("Zoom+ (2)", graph_area, SLOT (change_zoom_up ()));
  action->setShortcut (QString ("2"));
  action = tool_bar->addAction ("Zoom- (3)", graph_area, SLOT (change_zoom_down ()));
  action->setShortcut (QString ("3"));
  
  action = tool_bar->addAction ("n *= 2 (4)", graph_area, SLOT (change_n_up ()));
  action->setShortcut (QString ("4"));
  action = tool_bar->addAction ("n /= 2 (5)", graph_area, SLOT (change_n_down ()));
  action->setShortcut (QString ("5"));
  
  action = tool_bar->addAction ("p++ (6)", graph_area, SLOT (change_disturbance_up ()));
  action->setShortcut (QString ("6"));
  action = tool_bar->addAction ("p-- (7)", graph_area, SLOT (change_disturbance_down ()));
  action->setShortcut (QString ("7"));

  action = tool_bar->addAction ("E&xit", window, SLOT (close ()));
  action->setShortcut (QString ("Ctrl+X"));

  tool_bar->setMaximumHeight (30);

  window->setMenuBar (tool_bar);
  window->setCentralWidget (graph_area);
  window->setWindowTitle ("Graph");

  window->show ();
  app.exec ();
  
  graph_area->delete_memory ();
  
  delete window;
  return 0;
}
