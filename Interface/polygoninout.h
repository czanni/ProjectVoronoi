#ifndef POLYGONINOUT
#define POLYGONINOUT

#include <QPointF>
#include <string>
#include <vector>
#include <Graph.h>

namespace Medial{

void savePoints(std::string filename, Graph &extPoints, Graph &intPoints);
std::vector<Graph> loadSlices(std::string filename);
Graph loadPoints(std::string filename, bool externe);

}

#endif // POLYGONINOUT

