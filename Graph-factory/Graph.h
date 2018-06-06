#ifndef GRAPH_H
#define GRAPH_H


#include <geogram_gfx/basic/GLSL.h>
#include <geogram_gfx/basic/GL.h>
#include <vector>
#include <array>
#include <memory>
#include <stdbool.h>
#include <stdio.h>
#include <map>


enum class treatment { inside, outside, unknown };
struct Neighbor {
    int index;
    GEO::vec2 closest;
};

class Graph
{
    //TODO: marquer les fonctions et input en const s'ils ne sont pas modifiés
    std::vector <std::vector <Neighbor> > m_connexions;
    std::vector <GEO::vec2> m_points;
    std::vector <int> m_infiniteConnections;
    std::vector <treatment> m_pointTreatment;
    std::map<std::pair<int,int>, int> m_neighbors;


public:
    Graph();
    void fixClosest(int i, int j, GEO::vec2 &close);
    std::map<std::pair<int,int>, int> getNeighbors();
    void addNeighbor(int i, int j, int pointIndex);
    void addEdge(const std::array<int,2> & connexion, bool checkIfAlreadyPresent=false );
    void addPoint (const GEO::vec2 &point);
    static std::unique_ptr <Graph> demoGraph(int n);
    std::vector<GEO::vec2> &getPoints();
    void addInfinite (GEO::index_t i);
    void writeToFile();
    int numVertex();
    std::vector<Neighbor> directAdjacency(int v);
    const std::vector<std::vector<Neighbor> > &getConnexions();
    const std::vector <int>& getInfiniteConnection();
    GEO::vec2 & getPointCoordinate(int i);
    bool existsEdge(int i, int k) const;
    void changeStatus(int i, treatment T);
    treatment getStatus(int i);
    void removeOutsidePoints();




};

#endif // GRAPH_H
