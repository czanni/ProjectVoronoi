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
void Graph::addEdge(std::array<int,2> connexion /*, bool checkIfAlreadyPresent = true */)
//adds an edge to the structure. By default, they work in both directions
{
    int from = connexion[0];
    int to =  connexion[1];
    assert(from<m_connexions.size());
    assert(to<m_connexions.size());

    if (std::find(m_connexions[from].begin(), m_connexions[from].end(), to )== m_connexions[from].end())
    {
        m_connexions[from].push_back(to);
    }


    if (std::find(m_connexions[to].begin(), m_connexions[to].end(), from )== m_connexions[to].end())
    {
        m_connexions[to].push_back(from);
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

std::vector <int> Graph::directAdjacency(int v) {
    return m_connexions[v];
}


void Graph::addPoint (GEO::vec2 point)
{
    m_points.push_back(point);
    m_connexions.push_back(std::vector<int>());

}

GEO::vec2 & Graph::getPointCoordinate(int i) {
    return m_points[i];
}


std::vector <GEO::vec2>& Graph::getPoints() {
    return m_points;
}

int Graph::getSize() {
    return m_points.size();
}


void Graph::addInfinite (GEO::index_t i) {
    infinite_connection.push_back(i);
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

const std::vector<std::vector<int> >& Graph::getConnexions() {
    return m_connexions;
}
/**
 *  function for deletion, i is the point to delete
 */

/*
 * Faire une fonction pour supprimer une liste de points, avec une map qui donne en temps constant l'index du nouveau point en fonction du précédent.
 */
/*
void Graph::remove(std::vector <int> toRemove) {

    // First, we create the map that will correspond each index to the next
    // We create a vector, same size as the graph, with at each spot will be the shift needed to correspond to the new graph. an index of -1 means to delete
    std::vector <int> shifts (m_points.size(), 0);

    for (const auto &p : toRemove) {
        // We go through every point we need to remove
        for (int i=0;i<shifts.size();++i) {
            // If we have to remove it : we put -1 in it
            if (i==p) {
                shifts[i] = -1;
            }
            //Else, if greater and not -1, we need to increase the shift
            else if (shifts[i]!=-1 && i>p) {
                shifts[i] = shifts[i] + 1;
            }
        }
    }

    //We remove one point
    m_points.pop_back();

    //Then, we reindex all the points in the connections
    for (int k=0;k<m_connexions.size;++k) {
        for (int j=0;j<m_connexions[k].size;++j){

           //If point index greater than i, decrease it by one
            if (m_connexions[k][j] > i) {
                m_points[k][j] = m_points[k][j]-1;
            }
        }

    }
}
*/

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
    return infinite_connection;
}




