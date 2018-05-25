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
#include <set>



namespace GraphMaker {

typedef GEO::vector<GEO::vec2> Polygon;


    //MUST USE THIS FUNCTION BEFORE ANY USE OF THE LIBRARY
    void initialize();

    std::unique_ptr<Graph> extractMedialAxis(Graph & inputGraph);
    void removeOutsidePoints(std::unique_ptr<Graph> voronoiGraph);
    std::unique_ptr<Graph> makeMorePoints(Graph &inputGraph, const float step=9);
    std::unique_ptr<Graph> extractVoronoi(Graph &inputGraph, std::set<std::pair<int,int>> &intersects);
    std::map <std::pair<int,int>, bool> voronoiIntersection (Graph& inputGraph);
    void fixOutsidePoints(Graph &inputGraph, std::set<std::pair<int, int> > &voronoiIntersection);
    std::unique_ptr<Graph> extractDelaunay(Graph &inputGraph_small, float step=50);

    //Auxiliary functions using GEOGRAM
    void getVoronoiCell(GEO::index_t t0, GEO::index_t lv, Polygon& cell);
    GEO::vec2 circumcenter(GEO::index_t t) ;
    GEO::index_t findVertex(GEO::index_t t, GEO::index_t v);
    GEO::vec2 infiniteVertex(GEO::index_t t, GEO::index_t e);

}


#endif // GRAPHMAKER_H
