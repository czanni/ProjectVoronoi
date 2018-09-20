#include "Graph.h"

#include <cassert>

#include <map>
#include <algorithm>
#include <functional> // std::unary_function

namespace GraphMaker {

/**
 * Edges are stored in form of an adjacency matrix, as edges work in both ways,
 * They are added both ways in the adjacency matrix
 * By default, adding an edges checks the edge isn't already stored in the
 * matrix, you can change the boolean status to simply add the edge to the adjacency matrix
 * without any research
 */
void Graph::addEdge(const std::array<int,2> & connexion, bool checkIfAlreadyPresent)
//adds an edge to the structure. By default, they work in both directions
{
    struct Neighbor from = {.index = connexion[0], .closest=GEO::vec2(0,0)};
    struct Neighbor to = {.index = connexion[1], .closest=GEO::vec2(0,0)};
    assert(from.index<m_connexions.size());
    assert(to.index<m_connexions.size());

    if (checkIfAlreadyPresent) {

        if (std::find_if(m_connexions[from.index].begin(), m_connexions[from.index].end(),
                         std::not1(std::function<bool (Neighbor)>([&connexion](Neighbor i){ return i.index == connexion[1]; }))
                            ) == m_connexions[from.index].end())
        {
            m_connexions[connexion[0]].push_back(to);
        }


        if (std::find_if(m_connexions[to.index].begin(), m_connexions[to.index].end(),
                         std::not1(std::function<bool (Neighbor)>([&connexion](Neighbor i){ return i.index == connexion[0]; }))
                            ) == m_connexions[to.index].end())
        {
            m_connexions[connexion[1]].push_back(from);

        }
    }

    else {
        m_connexions[connexion[0]].push_back(to);
        m_connexions[connexion[1]].push_back(from);
    }

  //END of checking if present
  /*  else {
        m_connexions[to].push_back(from);
        m_connexions[from].push_back(to);
    } //Less operations*/

}

void Graph::fixClosest(int i, int j, const GEO::vec2 &point) {
    for (Neighbor &N : m_connexions[i]){
        if (N.index == j) {
            N.closest = point;
            return;
        }
    }
    return;
}

void Graph::setIntersect(int i, int j) {
    for (Neighbor &N : m_connexions[i]){
        if (N.index == j) {
            N.intersect = true;
            break;
        }
    }
    for (Neighbor &N : m_connexions[j]){
        if (N.index == i) {
            N.intersect = true;
            break;
        }
    }
}
bool Graph::intersect(int i, int j) const {
    for (const Neighbor &N : m_connexions[i]){
        if (N.index == j && N.intersect) {
            return true;
        }
    }
    return false;
}


void Graph::addPoint (const GEO::vec2 &point)
{
    m_points.push_back(point);
    m_connexions.push_back(std::vector<Neighbor>());
    m_pointTreatment.push_back(treatment::unknown);

}

GEO::vec2 & Graph::getPointCoordinate(int i) {
    return m_points[i];
}

void Graph::addInfinite (GEO::index_t i) {
    m_infiniteConnections.push_back(i);
}

bool Graph::existsEdge(int i, int k) const {
    assert(i<m_connexions.size());
    assert(k<m_connexions.size());

    for (auto N : m_connexions[i]) {
        if (N.index == k) {
            return true;
        }
    }
    return false;


}

//TODO merge at the same time point the map to same integer coordinates!
void Graph::removeOutsidePoints(){
	std::map <int, int> oldToNewIndex;
	int i = 0;
	int j = numVertex()-1;
	int originalPos = 0;
	while (i<=j) {
		if (getStatus(i) == treatment::inside){
			//The point is to be kept, then we add its index to the table
			oldToNewIndex.insert(std::pair<int,int> (originalPos,i));
			++i;
			originalPos = i;
		} else {
			if( i < j ) {
				//Change point coordinate
				std::swap(m_points[j], m_points[i]);
				//Change connexions
				std::swap(m_connexions[j], m_connexions[i]);
				//Change status
				std::swap(m_pointTreatment[j], m_pointTreatment[i]);
			}
			originalPos = j;
			--j;
		}
	}

	m_points.resize(i);
	m_connexions.resize(i);
	m_pointTreatment.resize(i);

	for (int k= 0; k<m_connexions.size();++k ) {
		int l =0;
		while (l<m_connexions[k].size()) {
			const auto result = oldToNewIndex.find(m_connexions[k][l].index);
			if (result == oldToNewIndex.end() ) {
				m_connexions[k].erase(m_connexions[k].begin() + l);
			} else {
				m_connexions[k][l].index = result->second;
				++l;
			}

		}
	}
}

} // end of namespace GraphMaker
