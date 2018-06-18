#ifndef POLYGONINOUT
#define POLYGONINOUT

#include <Graph.h>

#include <clipper.hpp>

#include <string>
#include <vector>

std::vector<ClipperLib::Paths> loadSlices(const std::string& filename, int &nbSlice);

#endif // POLYGONINOUT

