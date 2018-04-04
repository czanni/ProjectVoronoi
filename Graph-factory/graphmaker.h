#ifndef GRAPHMAKER_H
#define GRAPHMAKER_H

#include <geogram_gfx/basic/GLSL.h>
#include <geogram_gfx/basic/GL.h>
#include <geogram_gfx/GLUP/GLUP.h>
#include <geogram_gfx/glup_viewer/glup_viewer.h>
#include <geogram/delaunay/delaunay.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/logger.h>
#include <geogram/numerics/predicates.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include "Graph.h"


namespace GraphMaker {

typedef GEO::vector<GEO::vec2> Polygon;

//Tout mettre en private sauf extractvoronoi et infinitevertex et voronoiIntersection

    GEO::index_t findVertex(GEO::index_t t, GEO::index_t v);
    GEO::vec2 infiniteVertex(GEO::index_t t, GEO::index_t e);
    void getVoronoiCell(GEO::index_t t0, GEO::index_t lv, Polygon& cell);
    GEO::vec2 circumcenter(GEO::index_t t) ;
    void initialize();
    std::unique_ptr<Graph> removeOutsidePoints(Graph &voronoiGraph);
    std::unique_ptr<Graph> makeMorePoints(Graph &inputGraph, float step=5);
    std::unique_ptr <Graph> extractVoronoi(Graph &inputGraph);
    std::map <std::pair<int,int>, bool> voronoiIntersection (Graph& inputGraph);
    void fixOutsidePoints(Graph &inputGraph);
    std::unique_ptr<Graph> extractDelaunay(Graph &inputGraph_small, float step=5);


}


#endif // GRAPHMAKER_H
