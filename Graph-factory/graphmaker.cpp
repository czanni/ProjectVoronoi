#include "graphmaker.h"
#include <geogram/delaunay/delaunay.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/logger.h>
#include <geogram/numerics/predicates.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram_gfx/basic/GLSL.h>
#include <geogram_gfx/basic/GL.h>
#include <cmath>
#include <map>
#include "Graph.h"


using namespace GEO;
/*
 *
 * Faire une map qui mappe une couleur à un point pour le parcours de graph, pour vérifier qu'on voit tous les points, et dans le bon ordre.
 * (premier temps : map globale)
 *
 *
 */



namespace GraphMaker {


typedef GEO::vector<GEO::vec2> Polygon;

GEO::Delaunay_var delaunay; //the delaunay polygon


void initialize() {


    GEO::initialize();
    GEO::Logger::instance()->set_quiet(false);

    GEO::CmdLine::import_arg_group("standard");
    GEO::CmdLine::import_arg_group("algo");
    GEO::CmdLine::import_arg_group("gfx");

    GEO::CmdLine::set_arg("sys:assert","abort");

    GEO::Graphics::initialize();
}

GEO::vec2 circumcenter(index_t t) {

    signed_index_t v1 = delaunay->cell_to_v()[3*t];
    signed_index_t v2 = delaunay->cell_to_v()[3*t+1];
    signed_index_t v3 = delaunay->cell_to_v()[3*t+2];
    vec2 p1(delaunay->vertex_ptr(index_t(v1)));
    vec2 p2(delaunay->vertex_ptr(index_t(v2)));
    vec2 p3(delaunay->vertex_ptr(index_t(v3)));
    return Geom::triangle_circumcenter(p1,p2,p3);
}

float distance(GEO::vec2 &p1,GEO::vec2 &p2) {
    float res=0;
    for (int i=0;i<2;++i) {
        res += (p1[i]-p2[i])*(p1[i]-p2[i]);
    }
    return std::sqrt(res);
}

std::unique_ptr<Graph> makeMorePoints(Graph &inputGraph, float const step) {
    /*
     * [DONE]
     * Faire un set pour tester si un edge est déjà parcourru : stocker les edges en (i,j) avec i<j
     * std::set<std::pair<int, int>> s
     * s.find(std::make_pair(i,j)) != s.end()
     * utiliser s.emplace(i,j) évite les recopies
     */

    std::set<std::pair<int, int>> visitedEdges;
    std::unique_ptr<Graph> out = std::make_unique<Graph>();

    for (auto p: inputGraph.getPoints()) {
        out -> addPoint(p);
    }

    //We compute the distance we went from the point on border, then we compare it to the distance to the point adjacent to this one.
    //Since we go in the same direction, if the distance squared is greater, then we went too far and must change vertex.

    for (int currentVertex = 0; currentVertex < inputGraph.numVertex() ;++currentVertex) {

        for (Neighbor nextVertex_N : inputGraph.directAdjacency(currentVertex)) {
            int nextVertex = nextVertex_N.index;
            if (visitedEdges.find(std::make_pair(currentVertex, nextVertex)) == visitedEdges.end()) {

                float distFromConnectedVertex = distance(out -> getPointCoordinate(currentVertex), out -> getPointCoordinate(nextVertex));
                GEO::vec2 direction = ( out -> getPointCoordinate(nextVertex) - out -> getPointCoordinate(currentVertex));

                //If we moved, we check the distance to the current point, to see if we aren't farther away than needed

                //Step is here a density parameter (the density : step is to adjust globally to fit the project requirements)
                int nStep = (int) ((distFromConnectedVertex)/(step)) + 1 ;


                //first vertex of the edge (which is a little different)
                if (nStep == 1) {
                                    out -> addEdge({currentVertex, nextVertex});
                                }
                else {

                    out -> addPoint(inputGraph.getPointCoordinate(currentVertex) + direction/ (float) nStep);
                    out -> addEdge({currentVertex, out -> numVertex() -1 });

                    for (int k=2;k<nStep;++k) {

                        out -> addPoint(inputGraph.getPointCoordinate(currentVertex) + direction*((float) k / (float) nStep));
                        out -> addEdge({out -> numVertex()-1, out -> numVertex() - 2});

                    }

                    out -> addEdge({out -> numVertex() - 1, nextVertex });
                }


                //mark edge as visited
                visitedEdges.emplace(currentVertex, nextVertex);
                visitedEdges.emplace(nextVertex, currentVertex);

            }
        }
    }


    return std::move(out);

}

bool isOut( const Graph& graph, index_t t, index_t e) {
    index_t lv1 = (e+1)%3;
    index_t lv2 = (e+2)%3;
    index_t v1 = index_t(delaunay->cell_to_v()[3*t+lv1]);
    index_t v2 = index_t(delaunay->cell_to_v()[3*t+lv2]);

    return !(graph.existsEdge(v1,v2)); // t et e sur le diagramme de voronoi, v1 v2 sur le graphe, calculer le projeté en produit scalaire et voir le plus proche aux deux edges.

}

/** Param :
 * [OUT] An empty graph that will contain the voronoi graph
 * [IN] the step (in arbitrary units, in function of the chosen unit in the Graph) of the cutting of the Graph, in order ot make it closer to a graph of segment.
 *
 */
std::unique_ptr<Graph> extractVoronoi(Graph &inputGraph, std::set<std::pair<int,int>> &intersects) {
//    std::unique_ptr <Graph> voronoiGraph ;
    //Then, we call the function makeMorePoints to add the necessary points to be closest to the segment voronoi diagram.

 //   auto ptrGraph = makeMorePoints(inputGraph_small, step);
//    Graph inputGraph = *makeMorePoints(inputGraph_small, step);


    //return makeMorePoints(inputGraph_small, step);

    auto voronoiGraph = std::make_unique<Graph>();
    delaunay = Delaunay::create(2,"BDEL2d");
    delaunay->set_vertices(inputGraph.getPoints().size(), &(inputGraph.getPoints().data()->x));
    //CZ : delaunay is a smart pointer, overator-> allow to have access to the dereferenced inner pointer
    //CZ : set_vertices actually recompute the cells !
    // data is the pointer to the dynamic array stored inside the std::vector
    // we can note that once we have done
    //      &points.data()->x
    // we only obtain a pointer to a double (we do not know anymore about the vec2 class !
    // basically this is only valid if the storage for the vec2 class only contains the memory of two doubles and nothing more
    // (which is the case here, this would not be anymore the case if the class vec2 was containing other attribut or was virtual)
    //CZ : we only provide the number of points and not the number of coordinates as delaunay is actually a Delaunay2D.


    //We add the points in the Voronoi Graph
    for(index_t t=0; t<delaunay->nb_cells(); ++t) {
        voronoiGraph -> addPoint(circumcenter(t));
        //Index in graph are the same as inb cell indexes in voronoi
    }
    //Then, we add the link corresponding to each cell

    for(index_t t=0; t<delaunay->nb_cells(); ++t) {
        //function from display_edges
        for(index_t e=0; e<3; ++e) {
            signed_index_t t2 = delaunay->cell_to_cell()[3*t+e];

            bool out = isOut(inputGraph,t,e);
            if(t2 == -1) {
                voronoiGraph -> addInfinite(t); //add to the infinite matrix the vertex connected to the inf vertex.
                //TODO: évaluer quel edge rencontre un edge du graph donné en input et le stocker.
                if (out) {
                    voronoiGraph -> changeStatus(t,treatment::outside);
                }
                else {
                    voronoiGraph -> changeStatus(t,treatment::inside);
                }

            } else if(t2 >signed_index_t(t)) {
                voronoiGraph ->addEdge({(int) t,(int)t2});
                //If the edge exists, we put one of the closest point on the delaunay triangle, i.e the edge of delaunay who's dual it is.
                voronoiGraph -> fixClosest(t,t2, inputGraph.getPointCoordinate(delaunay -> cell_to_v()[3*t + (e+1)%3]));

                if(!out) {
                    intersects.emplace((int)t,(int)t2);
                    intersects.emplace((int)t2,(int)t);
                }
            }
        }
    }
    return std::move(voronoiGraph);
}

std::unique_ptr<Graph> extractDelaunay(Graph &inputGraph, float step) {
    auto delaunayGraph = std::make_unique<Graph>();
    delaunay = Delaunay::create(2,"BDEL2d");
    delaunay->set_vertices(inputGraph.getPoints().size(), &(inputGraph.getPoints().data()->x));
    int newindex=0;
    for(index_t c=0; c<delaunay->nb_cells(); ++c) {
                const signed_index_t* cell = delaunay->cell_to_v() + 3*c;
                for(index_t e=0; e<3; ++e) {
                    signed_index_t v1 = cell[e];
                    signed_index_t v2 = cell[(e+1)%3];
                    delaunayGraph -> addPoint((GEO::vec2) delaunay->vertex_ptr(index_t(v1)));
                    delaunayGraph -> addPoint((GEO::vec2) delaunay->vertex_ptr(index_t(v2)));
                    delaunayGraph -> addEdge({newindex, newindex+1});
                    newindex+=2;
                }
            }

    return std::move(delaunayGraph);
}





std::map <std::pair<int,int>, bool> voronoiIntersection (Graph& inputGraph) {
    delaunay = Delaunay::create(2,"BDEL2d");
    delaunay->set_vertices(inputGraph.getPoints().size(), &(inputGraph.getPoints().data()->x));

    std::map <std::pair<int,int>, bool> intersects;
    for(index_t t=0; t<delaunay->nb_cells(); ++t) {
        //function from display_edges
        for(index_t e=0; e<3; ++e) {
            signed_index_t t2 = delaunay->cell_to_cell()[3*t+e];
            index_t lv1 = (e+1)%3;
            index_t lv2 = (e+2)%3;
            index_t v1 = index_t(delaunay->cell_to_v()[3*t+lv1]);
            index_t v2 = index_t(delaunay->cell_to_v()[3*t+lv2]);

            if (inputGraph.existsEdge(v1,v2)) {
                intersects.emplace(std::make_pair((int)v1,(int)v2),true);
                intersects.emplace(std::make_pair((int)v2,(int)v1),true);
            }
        }
    }
    return intersects;
}

void depthSearch_outside(Graph & voronoiGraph,int previous, int next,  std::vector<bool> & visited, std::set <std::pair<int,int>> voronoiIntersection){

    if  (voronoiGraph.getStatus(previous)==treatment::unknown ) {
        return;
    }
    (visited)[next]=true;
    //If the edge intersects the border we reverse status
    if (voronoiIntersection.find(std::pair<int,int>(previous,next)) != voronoiIntersection.end()) {
        if (voronoiGraph.getStatus(previous)==treatment::outside) {
            voronoiGraph.changeStatus(next, treatment::inside);
        }
        else {
            voronoiGraph.changeStatus(next, treatment::outside);
        }
    }
    //Else, we make it the stame
    else {
        voronoiGraph.changeStatus(next,voronoiGraph.getStatus(previous));
    }
    for (Neighbor i:voronoiGraph.directAdjacency(next)) {
        if (!((visited)[i.index])) {
            depthSearch_outside(voronoiGraph,next,i.index, visited, voronoiIntersection);
        }
    }


}



void fixOutsidePoints(Graph &voronoiGraph, std::set <std::pair<int,int>> & voronoiIntersection) {
     if (voronoiGraph.numVertex()>0) {
         //As we will compute the search on different points and we don't want to visit each points several times, we store visited as a ptr

         std::vector<bool> visited(voronoiGraph.numVertex(), false);

         //Depth first search on outside points
         //Get the list of points to begin with (might be only one, and can be empty)
         std::vector<int> initialPoints = {};
         for (int i=0;i<voronoiGraph.numVertex();++i) {
             if (voronoiGraph.getStatus(i) !=treatment::unknown) {
                 initialPoints.push_back(i);
             }
         }
         for (int i:initialPoints) { //Parcourir sur TOUS les points et puis tester si outside ET not visited
             for (Neighbor k:voronoiGraph.directAdjacency(i)) {
                 depthSearch_outside(voronoiGraph,i,k.index,visited,voronoiIntersection);
             }
         }
     }
     voronoiGraph.removeOutsidePoints();

}




void getVoronoiCell(index_t t0, index_t lv, Polygon& cell) {

    cell.resize(0);
    index_t v = index_t(delaunay->cell_to_v()[3*t0+lv]);
    bool on_border = false;
    index_t t = t0;

    // First, we turn around the vertex v. To do that, we compute
    // lv, the local index of v in the current triangle. Following
    // the standard numerotation of a triangle, edge number lv is
    // not incident to vertex v. The two other edges (lv+1)%3 and
    // (lv+2)%3 of the triangle are indicent to vertex v. By always
    // traversing (lv+1)%3, we turn around vertex v.
    do {
        index_t e = (lv+1)%3;
        signed_index_t neigh_t = delaunay->cell_to_cell()[3*t+e];
        if(neigh_t == -1) {
            on_border = true;
            break;
        }
        cell.push_back(circumcenter(t));
        t = index_t(neigh_t);
        lv = findVertex(t,v);
    } while(t != t0);


    // If one traversed edge is on the border of the convex hull, then
    // we empty the cell, and start turing around the vertex in the other
    // direction, i.e. by traversing this time edge (lv+2)%3 until we
    // reach the other edge on the border of the convex hull that is
    // incident to v.
    if(on_border) {
        cell.resize(0);
        cell.push_back(infiniteVertex(t,(lv + 1)%3));
        for(;;) {
            cell.push_back(circumcenter(t));
            index_t e = (lv+2)%3;
            signed_index_t neigh_t = delaunay->cell_to_cell()[3*t+e];
            if(neigh_t == -1) {
                cell.push_back(infiniteVertex(t, e));
                break;
            }
            t = index_t(neigh_t);
            lv = findVertex(t,v);
        }
    }

    Polygon clipped;
    //convex_clip_polygon(cell, border, clipped);
    cell.swap(clipped);
}

index_t findVertex(index_t t, index_t v) {
    for(index_t lv=0; lv<3; ++lv) {
        if(index_t(delaunay->cell_to_v()[3*t+lv]) == v) {
            return lv;
        }
    }
    geo_assert_not_reached;
}

std::unique_ptr<Graph> extractMedialAxis(Graph & inputGraph) {
    std::unique_ptr<Graph> medialAxis = std::make_unique<Graph>();
    std::set <std::pair<int,int>> edgeIntersects{};

    medialAxis = extractVoronoi(inputGraph, edgeIntersects);
    fixOutsidePoints(*medialAxis, edgeIntersects);
    return std::move(medialAxis);

}

/**
 * Standard graph depth search, works on cells of the voronoi diagram
 */
void graphSearchVoronoiCompute(Graph &inputGraph /*, FUNCTION AS ARGUMENT */ ) {
    //extractVoronoi(inputGraph);
    std::vector<bool> v_visited(delaunay->nb_vertices());
    Polygon cell;
    for(index_t t=0; t<delaunay->nb_cells(); ++t) {
        for(index_t lv=0; lv<3; ++lv) {
            index_t v = index_t(delaunay->cell_to_v()[3*t+lv]);
            if(!v_visited[v]) {

                v_visited[v] = true;
                //obtain the cell
                getVoronoiCell(t,lv,cell);
                for(index_t i=1; i+1<cell.size(); ++i) {
                    /* EXECUTE FUNCTION (on each point of a given cell) */
                }
            }
        }
    }
    //every vertex was visited
}


/**
 * @brief standardDepthFirstSearch_AUX
 * @param inputGraph
 * @param v [the current computed point]
 * @param visited [Vector of visited points]
 * @param processFct applies fct to every point, use identity if you don't want one
 */

template<typename TLambda>
void standardDepthFirstSearch_AUX (Graph &inputGraph, int v, std::vector<bool> visited, TLambda &&processFct) {
    visited[v] = true;
    for(const auto &p :  inputGraph.directAdjacency(v)) {
        processFct(p);
        standardDepthFirstSearch_AUX(inputGraph,p,visited, processFct);

    }
}

/**
 * @brief standardDepthFirstSearch
 * @param inputGraph
 * @param firstVertex
 *
 * Uses auxiliary function
 *
 */

template<typename TLambda>
void standardDepthFirstSearch(Graph &inputGraph, int firstVertex, TLambda &&processFct)
{
    int size = inputGraph.numVertex();
    std::vector<bool> visited(size);
    for (int i = 0; i < size ; i++)
        visited[i] = false;

    standardDepthFirstSearch_AUX(inputGraph, firstVertex , visited, processFct);
}


vec2 infiniteVertex(index_t t, index_t e) {
    index_t lv1 = (e+1)%3;
    index_t lv2 = (e+2)%3;
    index_t v1 = index_t(delaunay->cell_to_v()[3*t+lv1]);
    index_t v2 = index_t(delaunay->cell_to_v()[3*t+lv2]);
    vec2 p1(delaunay->vertex_ptr(v1));
    vec2 p2(delaunay->vertex_ptr(v2));
    vec2 n = normalize(p2-p1);
    n = vec2(n.y, -n.x);
    return 0.5*(p1+p2)+100000.0*n;

}








} //End namespace GraphMaker
