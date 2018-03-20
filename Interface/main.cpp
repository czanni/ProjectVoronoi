#include "window.h"

#include <QApplication>
#include <Graph.h>
#include <graphmaker.h>

int main(int argc, char *argv[])
{
    //Q_INIT_RESOURCE(Interface);
    std::unique_ptr<Graph> nGon ;
    int edgeNumber = 4;
    std::unique_ptr <Graph> voronoiGraph ;
    GraphMaker::initialize();
    QApplication app(argc, argv);
    Window window;
    window.show();
    return app.exec();
}
