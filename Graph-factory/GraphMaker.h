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
    std::pair<std::unique_ptr<Graph>, std::unique_ptr<Graph> > extractInsideOutsiteMedialAxes(
		    const ClipperLib::Paths & input, float density);

    template<typename GraphT>
    std::unique_ptr<GraphT>  copyGraph(Graph &graph)
    {
          auto res = std::make_unique<GraphT>();
		typedef ClipperLib::cInt cInt;

		//std::cerr << std::endl << "COPYING: ";
		for(const auto &p : graph.getPositions())
		{
			ClipperLib::IntPoint cp(static_cast<cInt>(::round(p[0])), static_cast<cInt>(::round(p[1])));
			res->addVertex(cp);
			//std::cerr << '('<<p<< ") --> " << cp << "\n";
		}

		int i=0;
		for(const auto &connexions : graph.getNeighbors())
		{
			for(const auto &n : connexions)
			{
				res->addNeighbor(i, n.index,
						ClipperLib::IntPoint(static_cast<cInt>(::round(n.closest[0])), static_cast<cInt>(::round(n.closest[1])))
						);
			}
			++i;
		}
		//std::cerr << std::endl;
          return res;
    }

}

#endif // GRAPHMAKER_H
