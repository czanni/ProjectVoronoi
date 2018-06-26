#ifndef GRAPHMAKER_H
#define GRAPHMAKER_H

#include "Graph.h"

#include "clipper.hpp"

namespace GraphMaker {

	using GraphMaker::Graph;

    /**
     * @brief initialize initialize geogram, THIS FUNCTION MUST BE CALLED BEFORE ANY USE OF THE LIBRARY
     */
    void initialize();

    /**
     * Extract approximate version of the medial axis.
     * if density is zero do not subdivide paths, otherwise, use target density during subdivision
     */
    std::unique_ptr<Graph> extractMedialAxis(const ClipperLib::Paths& inputPath,
                                             float density=0.0f);


    template<typename GraphT>
    std::unique_ptr<GraphT>  copyGraph(Graph &graph)
    {
          auto res = std::make_unique<GraphT>();

          for(const auto &p : graph.getPositions())
          {
            res->addVertex(ClipperLib::IntPoint(static_cast<int>(p[0]), static_cast<int>(p[1])));
          }

          int i=0;
          for(const auto &connexions : graph.getNeighbors())
          {
              for(const auto &n : connexions)
              {
                res->addNeighbor(i, n.index,
                                 ClipperLib::IntPoint(static_cast<int>(n.closest[0]), static_cast<int>(n.closest[1]))
                    );
              }
              ++i;
          }

          return res;
    }

}

#endif // GRAPHMAKER_H
