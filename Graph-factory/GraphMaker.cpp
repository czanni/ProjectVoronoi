#include "GraphMaker.h"

#include "Graph.h"
#include "Delaunay_psm.h"

#include <set>


//---------------------------------------------------------

using namespace GEO;

typedef GEO::vector<GEO::vec2> Polygon;

struct DelaunayHelper {

  index_t findVertex(index_t t, index_t v) {
      for(index_t lv=0; lv<3; ++lv) {
          if(index_t(m_delaunay->cell_to_v()[3*t+lv]) == v) {
              return lv;
          }
      }
      geo_assert_not_reached;
  }

  vec2 infiniteVertex(index_t t, index_t e) {
      index_t lv1 = (e+1)%3;
      index_t lv2 = (e+2)%3;
      index_t v1 = index_t(m_delaunay->cell_to_v()[3*t+lv1]);
      index_t v2 = index_t(m_delaunay->cell_to_v()[3*t+lv2]);
      vec2 p1(m_delaunay->vertex_ptr(v1));
      vec2 p2(m_delaunay->vertex_ptr(v2));
      vec2 n = normalize(p2-p1);
      n = vec2(n.y, -n.x);
      return 0.5*(p1+p2)+100000.0*n;

  }

  GEO::vec2 circumcenter(index_t t) {
      signed_index_t v1 = m_delaunay->cell_to_v()[3*t];
      signed_index_t v2 = m_delaunay->cell_to_v()[3*t+1];
      signed_index_t v3 = m_delaunay->cell_to_v()[3*t+2];
      vec2 p1(m_delaunay->vertex_ptr(index_t(v1)));
      vec2 p2(m_delaunay->vertex_ptr(index_t(v2)));
      vec2 p3(m_delaunay->vertex_ptr(index_t(v3)));
      return Geom::triangle_circumcenter(p1,p2,p3);
  }

  void getVoronoiCell(index_t t0, index_t lv, Polygon& cell) {

      cell.resize(0);
      index_t v = index_t(m_delaunay->cell_to_v()[3*t0+lv]);
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
          signed_index_t neigh_t = m_delaunay->cell_to_cell()[3*t+e];
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
              signed_index_t neigh_t = m_delaunay->cell_to_cell()[3*t+e];
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

  std::pair<index_t,index_t> contourPoints(index_t t, index_t e) {
      index_t lv1 = (e+1)%3;
      index_t lv2 = (e+2)%3;
      index_t v1 = index_t(m_delaunay->cell_to_v()[3*t+lv1]);
      index_t v2 = index_t(m_delaunay->cell_to_v()[3*t+lv2]);
      return std::make_pair(v1,v2);
  }
  bool isOut( const Graph& graph, index_t t, index_t e) {
      index_t lv1 = (e+1)%3;
      index_t lv2 = (e+2)%3;
      index_t v1 = index_t(m_delaunay->cell_to_v()[3*t+lv1]);
      index_t v2 = index_t(m_delaunay->cell_to_v()[3*t+lv2]);

      return !(graph.existsEdge(v1,v2)); // t et e sur le diagramme de voronoi, v1 v2 sur le graphe, calculer le projeté en produit scalaire et voir le plus proche aux deux edges.
  }


  Delaunay_var m_delaunay; //the delaunay polygon
};

//------------------------------------------------------------------

namespace GraphMaker {

//------------------------------------------------------------------

void initialize() {
    GEO::initialize();
    GEO::Logger::instance()->set_quiet(false);

    GEO::CmdLine::import_arg_group("standard");
    GEO::CmdLine::import_arg_group("algo");

    GEO::CmdLine::set_arg("sys:assert","abort");
}

//------------------------------------------------------------------


//[IN] the step (in arbitrary units, in function of the chosen unit in the Graph) of the cutting of the Graph, in order ot make it closer to a graph of segment.

std::unique_ptr<Graph> extractVoronoi(ClipperLib::Paths &inputPath, float density)
{
    // Densify contour if need be and create vector of point at the same time
    
    std::vector<GEO::vec2> points;
    std::set<std::pair<int,int>> edge;
    //todo : could reserve before ...
    
    for(auto &path : inputPath)
    {
        int nbPointPrev = points.size();

        GEO::vec2 prevPoint(static_cast<double>(path.back().X),
                            static_cast<double>(path.back().Y));

        for(auto &p : path)
        {
            GEO::vec2 currPoint(static_cast<double>(p.X),
                                static_cast<double>(p.Y));
            
            float segLength = currPoint.distance(prevPoint);
            GEO::vec2 direction = currPoint - prevPoint;
            
            int nStep = static_cast<int>( segLength/density )+1;
            
            if (nStep == 1 || density==0.0f) {
                points.push_back(currPoint);
            } else {
                for (int k=1; k<=nStep; ++k) {
                    points.push_back(prevPoint + direction*(static_cast<float>(k) / static_cast<float>(nStep)));
                }
            }

            prevPoint = currPoint;
        }
        
        int nbPoint = points.size();
        for(int i=nbPointPrev; i+1<nbPoint; ++i) {
            edge.emplace(i,i+1);
            edge.emplace(i+1,i);
        }
        edge.emplace(nbPoint-1,nbPointPrev);
        edge.emplace(nbPointPrev,nbPoint-1);
    }

    // Generate Delaunay triangulation
    
    DelaunayHelper delaunayHelper;
    delaunayHelper.m_delaunay = Delaunay::create(2,"BDEL2d");
    delaunayHelper.m_delaunay->set_vertices(points.size(), &(points.data()->x));

    //TODO : could directly build a simple graph (provided a fiew additionnal
    // data member are added (or passed as additional variable to function call)
    auto voronoiGraph = std::make_unique<Graph>();
        
    //We add the points in the Voronoi Graph
    for(index_t t=0; t<delaunayHelper.m_delaunay->nb_cells(); ++t) {
        voronoiGraph -> addPoint(delaunayHelper.circumcenter(t));
        //Index in graph are the same as inb cell indexes in voronoi
    }
    
    //Then, we add the link corresponding to each cell
        
    for(index_t t=0; t<delaunayHelper.m_delaunay->nb_cells(); ++t) {
        //function from display_edges
        for(index_t e=0; e<3; ++e) {
            signed_index_t t2 = delaunayHelper.m_delaunay->cell_to_cell()[3*t+e];
                
            bool out = (edge.find(delaunayHelper.contourPoints(t,e)) == edge.end());
            
            if(t2 == -1) {
                    voronoiGraph -> addInfinite(t); //add to the infinite matrix the vertex connected to the inf vertex.
                    if (out) {
                        voronoiGraph -> changeStatus(t,treatment::outside);
                    } else {
                        voronoiGraph -> changeStatus(t,treatment::inside);
                    }
                } else if(t2 >signed_index_t(t)) {
                    voronoiGraph ->addEdge({(int) t,(int)t2});

                    //If the edge exists, we put one of the closest point on the delaunay triangle, i.e the edge of delaunay who's dual it is.
                    int closestPoint_i = delaunayHelper.m_delaunay -> cell_to_v()[3*t + (e+1)%3];
                    voronoiGraph -> fixClosest(t,t2, points[closestPoint_i]);
                    voronoiGraph -> fixClosest(t2,t, points[closestPoint_i]);

                    if(!out) {
                        voronoiGraph->setIntersect((int)t,(int)t2);
                    }
                }
            }
        }
        return std::move(voronoiGraph);
    }

//------------------------------------------------------------------

void depthSearch_outside(Graph & voronoiGraph,
                         int previous, int next,
                         std::vector<bool> & visited)
{
    visited[next]=true;
    //If the edge intersects the border we reverse status
    if (voronoiGraph.intersect(previous,next)) {
        if (voronoiGraph.getStatus(previous)==treatment::outside) {
            voronoiGraph.changeStatus(next, treatment::inside);
        } else {
            voronoiGraph.changeStatus(next, treatment::outside);
        }
    } else { // we make it the same
        voronoiGraph.changeStatus(next,voronoiGraph.getStatus(previous));
    }
    for (const auto &i : voronoiGraph.getNeighbors(next)) {
        if (!(visited[i.index]) && voronoiGraph.getStatus(i.index)==treatment::unknown) {
            depthSearch_outside(voronoiGraph,next,i.index, visited);
        }
    }
}

//------------------------------------------------------------------

void fixOutsidePoints(Graph &voronoiGraph)
{
     if (voronoiGraph.numVertex()>0) {
         //As we will compute the search on different points and we don't want to visit each points several times, we store visited as a ptr

         std::vector<bool> visited(voronoiGraph.numVertex(), false);

         //Depth first search on outside points
         //Get the list of points to begin with (might be only one, and can be empty)
         
         for (int i=0;i<voronoiGraph.numVertex();++i)
         {
             if (!visited[i] && voronoiGraph.getStatus(i) != treatment::unknown)
             {
                 for (const auto &k:voronoiGraph.getNeighbors(i))
                 {
                    depthSearch_outside(voronoiGraph, i,k.index,visited);
                 }
             }
         }
     }
    
     voronoiGraph.removeOutsidePoints();
}

//------------------------------------------------------------------

std::unique_ptr<Graph> extractMedialAxis(ClipperLib::Paths& inputPath,
                                         float density)
{
  std::unique_ptr<Graph> medialAxis = extractVoronoi(inputPath, density);
  fixOutsidePoints(*medialAxis);
  return std::move(medialAxis);
}

//------------------------------------------------------------------

} //End namespace GraphMaker

