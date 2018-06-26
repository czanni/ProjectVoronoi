#ifndef GRAPH_H
#define GRAPH_H

#include "Delaunay_psm.h"

#include <vector>
#include <array>
#include <memory>

namespace GraphMaker {

enum class treatment { inside, outside, unknown };
struct Neighbor {
    int index;
    GEO::vec2 closest;
    bool intersect; //whether or not the edge intersect the surface
};

class Graph
{
public:
    Graph() = default;
    
    int numVertex() const {
        return m_points.size();
    }

    std::vector<GEO::vec2> &getPositions() {
        return m_points;
    }
    const std::vector<std::vector<Neighbor> >& getNeighbors() const {
        return m_connexions;
    }
    const std::vector<Neighbor>& getNeighbors(int v) const {
        return m_connexions[v];
    }
    const std::vector <int>& getInfiniteConnection() const {
        return m_infiniteConnections;
    }

    void setIntersect(int i, int j);
    bool intersect(int i, int j) const;

    void addNeighbor(int i, int j, int pointIndex);
    void addEdge(const std::array<int,2> & connexion, bool checkIfAlreadyPresent=false );
    void addPoint (const GEO::vec2 &point);
    void addInfinite (GEO::index_t i);
    GEO::vec2& getPointCoordinate(int i);
    bool existsEdge(int i, int k) const;

    void fixClosest(int i, int j, const GEO::vec2 &point);
    
    void changeStatus(int i, treatment T) {
        m_pointTreatment[i]=T;
    }
    treatment getStatus(int i) const {
        return m_pointTreatment[i];
    }

    void removeOutsidePoints();

protected:
    std::vector <std::vector <Neighbor> > m_connexions;
    std::vector <GEO::vec2> m_points;
    std::vector <int> m_infiniteConnections;
    std::vector <treatment> m_pointTreatment;
public:
    std::vector <bool> m_visited;
};

} // end of namespace GraphMaker

#endif // GRAPH_H
