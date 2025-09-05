#include <QtWidgets/QApplication>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QAction>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QMessageBox>
#include "window.hpp"
//#include "thread_func.hpp"
#include "func.hpp"



int main(int argc, char *argv[]){
    QApplication app(argc, argv);

    QMainWindow *window = new QMainWindow;
    QMenuBar *tool_bar = new QMenuBar(window);
    Window *graph_area = new Window(window);
    QAction *action;

    if(graph_area->parse_command_line(argc, argv)){
        printf("Usage : %s  a  b  c  d  nx  ny  mx  my  k  eps  maxit  p\n", argv[0]);
        QMessageBox::warning(0, "Wrong input arguments!", "Wrong input arguments!");
        return -1;
    }

    action = tool_bar->addAction("&0 Next func", graph_area, SLOT(next_func()));
    action->setShortcut(QString("0"));

    action = tool_bar->addAction("&1 Next graph", graph_area, SLOT(next_graph()));
    action->setShortcut(QString("1"));

    action = tool_bar->addAction("&2 Inc s", graph_area, SLOT(inc_s()));
    action->setShortcut(QString("2"));

    action = tool_bar->addAction("&3 Dec sc", graph_area, SLOT(dec_s()));
    action->setShortcut(QString("3"));

    action = tool_bar->addAction("&4 Inc n", graph_area, SLOT(inc_n()));
    action->setShortcut(QString("4"));

    action = tool_bar->addAction("&5 Dec n", graph_area, SLOT(dec_n()));
    action->setShortcut(QString("5"));

    action = tool_bar->addAction("&6 Inc p", graph_area, SLOT(inc_p()));
    action->setShortcut(QString("6"));

    action = tool_bar->addAction("&7 Dec p", graph_area, SLOT(dec_p()));
    action->setShortcut(QString("7"));

    action = tool_bar->addAction("&8 Inc m", graph_area, SLOT(inc_m()));
    action->setShortcut(QString("8"));

    action = tool_bar->addAction("&9 Dec m", graph_area, SLOT(dec_m()));
    action->setShortcut(QString("9"));

    action = tool_bar->addAction("E&xit", graph_area, SLOT(close()));
    action->setShortcut(QString("Ctrl+X"));

    tool_bar->setMaximumHeight(30);

    window->setMenuBar(tool_bar);
    window->setCentralWidget(graph_area);
    window->setWindowTitle("Graph");

    window->show();
    app.exec();
    graph_area->close();
    delete window;
    return 0;
}
