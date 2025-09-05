#include <QtWidgets/QApplication>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QAction>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QMessageBox>
#include "window.hpp"
#include "func.hpp"



int main(int argc, char *argv[]) {
    QApplication app(argc, argv);

    QMainWindow *window = new QMainWindow;
    QMenuBar *tool_bar = new QMenuBar(window);
    Window *graph_area = new Window(window);
    QAction *action;

    if (graph_area->parse_command_line(argc, argv)) {
        printf("Usage : %s  a  b  c  d  nx  ny  mx  my  func_id  eps  maxit  p\n", argv[0]);
        QMessageBox::warning(0, "Проверьте аргументы", "Некорректный вызов");
        return -1;
    }

    action = tool_bar->addAction("[0] Change function", graph_area, SLOT(change_func()));
    action->setShortcut(QString("0"));

    action = tool_bar->addAction("[1] Change mode", graph_area, SLOT(change_mode()));
    action->setShortcut(QString("1"));

    action = tool_bar->addAction("[2] Zoom +", graph_area, SLOT(zoom_in()));
    action->setShortcut(QString("2"));

    action = tool_bar->addAction("[3] Zoom -", graph_area, SLOT(zoom_out()));
    action->setShortcut(QString("3"));

    action = tool_bar->addAction("[4] Twice n", graph_area, SLOT(twice_n()));
    action->setShortcut(QString("4"));

    action = tool_bar->addAction("[5] Halve n", graph_area, SLOT(halve_n()));
    action->setShortcut(QString("5"));

    action = tool_bar->addAction("[6] Increase p", graph_area, SLOT(increase_p()));
    action->setShortcut(QString("6"));

    action = tool_bar->addAction("[7] Decrease p", graph_area, SLOT(decrease_p()));
    action->setShortcut(QString("7"));

    action = tool_bar->addAction("[8] Twice m", graph_area, SLOT(twice_m()));
    action->setShortcut(QString("8"));

    action = tool_bar->addAction("[9] Halve m", graph_area, SLOT(halve_m()));
    action->setShortcut(QString("9"));

    action = tool_bar->addAction("Exit", graph_area, SLOT(close()));
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
