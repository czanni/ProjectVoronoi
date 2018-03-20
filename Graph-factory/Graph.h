#ifndef GRAPH_H
#define GRAPH_H


#include <geogram_gfx/basic/GLSL.h>
#include <geogram_gfx/basic/GL.h>
#include <vector>
#include <array>
#include <memory>
#include <stdbool.h>
#include <stdio.h>

class Graph
{

    std::vector <std::vector <int> > m_connexions;
    std::vector <GEO::vec2> m_points;
    std::vector <int> infinite_connection;
public:
    Graph();
    void addEdge(std::array<int,2> connexion );
    void addPoint (GEO::vec2 point);
    static std::unique_ptr <Graph> demoGraph(int n);
    std::vector<GEO::vec2> &getPoints();
    void addInfinite (GEO::index_t i);
    void writeToFile();
    int getSize();
    std::vector <int> directAdjacency(int v);
    const std::vector <std::vector <int> >& getConnexions();
    const std::vector <int>& getInfiniteConnection();
    GEO::vec2 & getPointCoordinate(int i);



};

#endif // GRAPH_H
