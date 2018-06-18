#include "polygoninout.h"

#include <fstream>

using namespace ClipperLib;

std::vector<Paths> loadSlices(const std::string &filename, int &nbSlice){

    std::ifstream istrm(filename, std::ios::out);
    istrm >> nbSlice;

    std::vector<Paths> slices(nbSlice);

    float x;
    float y;
    for (int k=0;k<nbSlice;++k){
        int concav;
        istrm >> concav;
        Paths &paths = slices[k];
        paths.resize(concav);
        for (int k=0;k<concav;++k) {
            int nbrePoints; // le nombre de point par contour
            istrm >> nbrePoints;
            paths[k].resize(nbrePoints);
            for(int i = 0 ; i < nbrePoints ; i++){
                istrm >> x >> y;
                paths[k][i] = IntPoint(static_cast<int>(x),
                                       static_cast<int>(y));
            }
        }
    }
    return std::move(slices);
}


