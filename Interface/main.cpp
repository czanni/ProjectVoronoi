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
    /*Graph debug;
    debug.addPoint({2,3});
    debug.addPoint({1,4});
    debug.addPoint({2,2});
    debug.addPoint({1,1});
    debug.addPoint({3,1});
    debug.addPoint({3,0});

    debug.addEdge({0,1});
    debug.addEdge({1,4});
    debug.addEdge({0,2});
    debug.addEdge({2,3});
    debug.addEdge({2,4});
    debug.addEdge({4,5});

    debug.changeStatus(0, treatment(inside));
    debug.changeStatus(1, treatment(outside));
    debug.changeStatus(2, treatment(inside));
    debug.changeStatus(3, treatment(outside));
    debug.changeStatus(4, treatment(outside));
    debug.changeStatus(5, treatment(outside));

    debug.removeOutsidePoints();*/
    QApplication app(argc, argv);
    Window window;
    window.show();
    return app.exec();
}
