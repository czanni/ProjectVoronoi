#include "Graph.h"
#include <vector>
#include <geogram_gfx/basic/GLSL.h>
#include <geogram_gfx/basic/GL.h>
#include <geogram_gfx/GLUP/GLUP.h>
#include <geogram_gfx/glup_viewer/glup_viewer.h>
#include <geogram/delaunay/delaunay.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/logger.h>
#include <geogram/numerics/predicates.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <iostream>
#include <fstream>
#include <stdbool.h>
#include <stdio.h>
#include <functional>   // std::unary_function

#include <cmath>
#define PI 3.14159265

#include <algorithm> // To avoid implementing the find function

#include <cassert> // For assert

Graph::Graph() //Creates empty graph
{

}


/**
 * Edges are stored in form of an adjacency matrix, as edges work in both ways,
 * They are added both ways in the adjacency matrix
 * By default, adding an edges checks the edge isn't already stored in the
 * matrix, you can change the boolean status to simply add the edge to the adjacency matrix
 * without any research
 */
void Graph::addEdge(const std::array<int,2> & connexion /*, bool checkIfAlreadyPresent = true */)
//adds an edge to the structure. By default, they work in both directions
{
    struct Neighbor from = {.index = connexion[0], .closest=GEO::vec2(0,0)};
    struct Neighbor to = {.index = connexion[1], .closest=GEO::vec2(0,0)};
    assert(from.index<m_connexions.size());
    assert(to.index<m_connexions.size());

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
    else {
        assert(false);
    }
  //END of checking if present
  /*  else {
        m_connexions[to].push_back(from);
        m_connexions[from].push_back(to);
    } //Less operations*/

}

void Graph::fixClosest(int i, int j, GEO::vec2 close) {
    for (Neighbor &N : m_connexions[i]){
        if (N.index == j) {
            N.closest = close;
            return;
        }
    return;
    }
}


std::vector <Neighbor> Graph::directAdjacency(int v) {
    return m_connexions[v];
}



std::map<std::pair<int,int>, int> Graph::getNeighbors() {
    return m_neighbors;
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


std::vector <GEO::vec2>& Graph::getPoints() {
    return m_points;
}

int Graph::numVertex() {
    return m_points.size();
}


void Graph::addInfinite (GEO::index_t i) {
    m_infiniteConnections.push_back(i);
}


std::unique_ptr <Graph> Graph::demoGraph(int n) {
   std::unique_ptr <Graph> demo = std::unique_ptr <Graph> (new Graph());
    demo -> addPoint({0.4,0});
    for (int i =1 ; i<n ; i++ ) {
        demo -> addPoint(   GEO::vec2(0.25, 0.25) +
                            GEO::vec2(GEO::Numeric::random_float64()*0.5, GEO::Numeric::random_float64()*0.5) );   //{cos( 2*PI*i/((float) n)) * 0.4,sin( 2*PI*i/((float) n))*0.4}
        demo->addEdge({i-1,i});
    }
    return std::move(demo);
}

void Graph::writeToFile() {
    std::ofstream myfile;
     myfile.open ("/home/ayme/Qt-projects/ProjetVoronoi/output.txt");

     //points
     for(const auto &p :  m_points) {
          std::ostringstream buffer ;
          buffer << p.x << " " << p.y << "\n";
          myfile << buffer.str();
     }     
     myfile.close();
}

const std::vector<std::vector<Neighbor> >& Graph::getConnexions() {
    return m_connexions;
}

bool Graph::existsEdge(int i, int k) const {
    assert(i<m_connexions.size());
    assert(k<m_connexions.size());

    if (std::find_if(m_connexions[i].begin(), m_connexions[i].end(),
                     std::not1(std::function<bool (Neighbor)>([k](Neighbor i){ return i.index == k; }))
                        ) == m_connexions[i].end())
    {
        return false;
    }
    return true;
}

void Graph::changeStatus(int i, treatment T) {
    m_pointTreatment[i]=T;
}

treatment Graph::getStatus(int i) {
    return m_pointTreatment[i];
}



void Graph::removeOutsidePoints(){
    std::map <int, int> oldToNewIndex;
    int i = 0;
    int j = numVertex()-1;
    int swapCount = 0;
    while (i<=j) {
        if (getStatus(i) == treatment::inside){
            //The point is to be kept, then we add its index to the table, if it hasn't been swaped, then it must stay in place, otherwise we put its old index
            if (swapCount == 0) {
                oldToNewIndex.insert(std::pair<int,int> (i,i));
            }
            else {
                oldToNewIndex.insert(std::pair<int,int> (j+1,i));
                swapCount = 0;
            }
            ++i;
        }
        else {
            //Change point coordinate
            std::swap(m_points[j], m_points[i]);

            //Change connexions
            std::swap(m_connexions[j], m_connexions[i]);

            //Change status
            std::swap(m_pointTreatment[j], m_pointTreatment[i]);


            ++swapCount;
            --j;
        }
    }

    m_points.resize(i);
    m_connexions.resize(i);
    m_pointTreatment.resize(i);


    for (int k= 0; k<m_connexions.size();++k ) {
        int l =0;
        while (l<m_connexions[k].size()) {

            if (oldToNewIndex.find(m_connexions[k][l].index) == oldToNewIndex.end() ) {
                m_connexions[k].erase(m_connexions[k].begin() + l);
            }
            else {
                m_connexions[k][l].index = oldToNewIndex.find(m_connexions[k][l].index) -> second;
                ++l;
            }

        }
    }
}



/*
 * Point IN/OUT : partir des points à l'infini
 * Tester si dans le dual, l'edge correspondant à l'edge parcourru est de la form (i, i+1),
 * ou alors construire un set avec tous les edges (pour éviter les problèmes si on a deux contours)
 *
 *
 * Rajouter un param qui correspondra à l'ordre dans le parcours.
 *
 */

const std::vector <int>& Graph::getInfiniteConnection() {
    return m_infiniteConnections;
}






