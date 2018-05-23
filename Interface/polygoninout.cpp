#include "polygoninout.h"
#include <fstream>

namespace Medial{

void savePoints(std::string filename, Graph &extPoints, Graph &intPoints) //Pour enregistrer les graphes dans un fichier texte
{
    std::ofstream ostrm(filename, std::ios::out);
    // Écriture des différentes lignes dans le fichier
    // Pour le graphe du contours extérieur
    ostrm << extPoints.getPoints().size() << "\n";
    foreach(const GEO::vec2 &point, extPoints.getPoints()){
        ostrm << point[0] << " " << point[1] << "\n";
    }
    foreach(const std::vector<int> &connexion, extPoints.getConnexions()){
        ostrm << connexion.size() << "\n";
        foreach(int index, connexion){
            ostrm << index << " ";
        }
        ostrm << "\n";
    }
    ostrm << extPoints.getInfiniteConnection().size() << "\n";
    foreach(const int &infiniteConnection, extPoints.getInfiniteConnection()){
        ostrm << infiniteConnection << "\n";
    }
    // Pour le graphe du contours intérieur
    ostrm << intPoints.getPoints().size() << "\n";
    foreach(const GEO::vec2 &point, intPoints.getPoints()){
        ostrm << point[0] << " " << point[1] << "\n";
    }
    foreach(const std::vector<int> &connexion, intPoints.getConnexions()){
        ostrm << connexion.size() << "\n";
        foreach(int index, connexion){
            ostrm << index << " ";
        }
        ostrm << "\n";
    }
    ostrm << intPoints.getInfiniteConnection().size() << "\n";
    foreach(const int &infiniteConnection, intPoints.getInfiniteConnection()){
        ostrm << infiniteConnection << "\n";
    }
}

Graph loadPoints(std::string filename, bool externe){ //Pour charger les graphes à partir d'un fichier texte
    std::ifstream istrm(filename, std::ios::out);
    Graph graph;
    float x;
    float y;
    if (externe){ // Pour le graphe du contours extérieur
        int nbrePoints;
        istrm >> nbrePoints;
        for(int i = 0 ; i < nbrePoints ; i++){
            istrm >> x >> y;
            GEO::vec2 point(x, y);
            graph.addPoint(point);
        }
        int nbreConnexions;
        for(int j = 0 ; j < nbrePoints ; j ++){
            istrm >> nbreConnexions;
            for(int i = 0 ; i < nbreConnexions ; i++){
                istrm >> x;
                std::array<int,2> connexion({(int) x, j});
                graph.addEdge(connexion);
            }
        }
        int nbreInfiniteConnexions;
        istrm >> nbreInfiniteConnexions;
        for(int i = 0 ; i < nbreInfiniteConnexions ; i++){
            istrm >> x;
            graph.addInfinite((int) x);
        }
    }
    else{
        // D'abord on saute la partie du contour extérieur
        int nbrePoints;
        istrm >> nbrePoints;
        for(int i = 0 ; i < nbrePoints ; i++){
            istrm >> x >> y;
        }
        int nbreConnexions;
        for(int j = 0 ; j < nbrePoints ; j ++){
            istrm >> nbreConnexions;
            for(int i = 0 ; i < nbreConnexions ; i++){
                istrm >> x;
            }
        }
        int nbreInfiniteConnexions;
        istrm >> nbreInfiniteConnexions;
        for(int i = 0 ; i < nbreInfiniteConnexions ; i++){
            istrm >> x;
        }
        // Et maintenant, on traite la partie qui nous intéresse qui est celle du contours intérieur
        istrm >> nbrePoints;
        for(int i = 0 ; i < nbrePoints ; i++){
            istrm >> x >> y;
            GEO::vec2 point(x, y);
            graph.addPoint(point);
        }
        for(int j = 0 ; j < nbrePoints ; j ++){
            istrm >> nbreConnexions;
            for(int i = 0 ; i < nbreConnexions ; i++){
                istrm >> x;
                std::array<int,2> connexion({(int) x, j});
                graph.addEdge(connexion);
            }
        }
        istrm >> nbreInfiniteConnexions;
        for(int i = 0 ; i < nbreInfiniteConnexions ; i++){
            istrm >> x;
            graph.addInfinite((int) x);
        }
    }
    return graph;
}


std::vector<Graph> loadSlices(std::string filename){ //Pour charger les graphes à partir d'un fichier texte
    std::ifstream istrm(filename, std::ios::out);
    std::vector<Graph> Slices = std::vector<Graph>();
    int nbSlice;
    istrm >> nbSlice;
    float x;
    float y;
    for (int k=0;k<nbSlice;++k){
        int concav;
        istrm >> concav;
        Graph graph = Graph();
        for (int k=0;k<concav;++k) {
            int nbrePoints; //le nombre de point par concavité0
            istrm >> nbrePoints;
            for(int i = 0 ; i < nbrePoints ; i++){
                istrm >> x >> y;
                GEO::vec2 point(x/5+350, y/5+250);
                graph.addPoint(point);
            }
            for (int i=graph.numVertex()-nbrePoints; i < graph.numVertex()-1; ++i) {
                graph.addEdge({i,i+1});
            }
            graph.addEdge({graph.numVertex()-nbrePoints, graph.numVertex()-1});

        }
        Slices.push_back(graph);
    }
    return Slices;

}

}
