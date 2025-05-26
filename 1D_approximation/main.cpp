#include <QtWidgets/QApplication>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QAction>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QMessageBox>

#include "window.hpp"

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);

    QMainWindow *window = new QMainWindow;
    QMenuBar *tool_bar  = new QMenuBar(window);
    Window *graph_area  = new Window(window);
    QAction *action;

    if (graph_area->parse_command_line(argc, argv)) {
        QMessageBox::warning(0, "Wrong argument", "Usage: ./a.out  a  b  n  func_id");
        return -1;
    }

    if (graph_area->set_memory()) return 1000;
    graph_area->set_func();
    graph_area->set_mode();
    graph_area->set_arrays();
    

    // Изменить функции, когда напишу их
    action = tool_bar->addAction("Change function [0]", graph_area, SLOT(change_func()));
    action->setShortcut(QString("0"));

    action = tool_bar->addAction("Change view [1]", graph_area, SLOT(change_approx()));
    action->setShortcut(QString("1"));


    action = tool_bar->addAction("Zoom + [2]", graph_area, SLOT(zoom_in()));
    action->setShortcut(QString("2"));

    action = tool_bar->addAction("Zoom - [3]", graph_area, SLOT(zoom_out()));
    action->setShortcut(QString("3"));


    action = tool_bar->addAction("Twice n [4]", graph_area, SLOT(increase_n()));
    action->setShortcut(QString("4"));

    action = tool_bar->addAction("Halve n [5]", graph_area, SLOT(decrease_n()));
    action->setShortcut(QString("5"));


    action = tool_bar->addAction("Increase p [6]", graph_area, SLOT(increase_p()));
    action->setShortcut(QString("6"));

    action = tool_bar->addAction("Decrease p [7]", graph_area, SLOT(decrease_p()));
    action->setShortcut(QString("7"));


    action = tool_bar->addAction("Exit", window, SLOT(close()));
    action->setShortcut(QString("Ctrl+X"));


    tool_bar->setMaximumHeight(30);

    window->setMenuBar(tool_bar);
    window->setCentralWidget(graph_area);
    window->setWindowTitle("Graph");

    window->show();
    app.exec();

    graph_area->clean_memory();

    delete window;
    return 0;
}