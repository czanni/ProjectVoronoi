
#include "polygoninout.h"

#include <Graph.h>
#include <GraphMaker.h>

#if 0
#include <GLUT/glut.h>  // GLUT, includes glu.h and gl.h
#else
#include <GL/glut.h>
#endif

#include <chrono>
#include <ctime>

int n_slices = 0;
std::vector<ClipperLib::Paths> slices;

int SIZE_X = 640;
int SIZE_Y = 480;

int id = 0;
int win;
double zoom = 1.0;

void display()
{
    if( slices[id].size() == 0 ) {
        std::cout << "Slice " << id << " is empty." << std::endl;
    }
    auto medial = GraphMaker::extractMedialAxis(slices[id],60);

    glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Set background color to black and opaque
    glClear(GL_COLOR_BUFFER_BIT);         // Clear the color buffer

    int index = 0;

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glScaled(zoom, zoom, 1.0);

    glBegin(GL_LINES);
    glColor3f(0.2f, 0.2f, 0.6f);
    index = 0;
    for(auto &connexions : medial->getNeighbors()) {
        for(auto &n : connexions) {
          auto p1 = medial->getPointCoordinate(index);
          auto p2 = n.closest;

            glVertex2f(p1[0]/(5*SIZE_Y),
                       p1[1]/(5*SIZE_Y));
            glVertex2f(p2[0]/(5*SIZE_Y),
                       p2[1]/(5*SIZE_Y));
        }
        index+=1;
    }
    glEnd();

    glBegin(GL_LINES);
    glColor3f(0.0f, 1.0f, 0.0f);
    for(auto &path : slices[id])
    {
        auto prevPoint = path.back();
        for(auto it=path.begin(); it!=path.end(); ++it)
        {
            auto currPoint = *it;
            glVertex2f((float)prevPoint.X/(float)(5*SIZE_Y),
                       (float)prevPoint.Y/(float)(5*SIZE_Y));
            glVertex2f((float)currPoint.X/(float)(5*SIZE_Y),
                       (float)currPoint.Y/(float)(5*SIZE_Y));
            prevPoint=currPoint;
        }
    }
    glEnd();

    // Draw medial axis
    glBegin(GL_LINES);
    glColor3f(1.0f, 0.0f, 0.0f);

    index = 0;
    for(auto &connexions : medial->getNeighbors())
    {
        for(auto &n : connexions)
        {
            glVertex2f(medial->getPositions()[n.index][0]/(5*SIZE_Y),
                       medial->getPositions()[n.index][1]/(5*SIZE_Y));
            glVertex2f(medial->getPositions()[index][0]/(5*SIZE_Y),
                       medial->getPositions()[index][1]/(5*SIZE_Y));
        }
        index+=1;
    }
    glEnd();

    glFlush();  // Render now
}

void processNormalKeys(unsigned char key, int /*xx*/, int /*yy*/) {
    switch( key ) {
        case 27: // Escape key
            glutDestroyWindow(win);
            exit (0);
            break;
        case '+':
            zoom *= 1.414141414;
            break;
        case '-':
            zoom /= 1.414141414;
            break;
        case 'b':
            id=(id+n_slices-1)%n_slices;
            break;
        default:
            id=(id+1)%n_slices;
    }
    glutPostRedisplay() ;
}

int main(int argc, char *argv[])
{
    GraphMaker::initialize();

    std::string filename;
    if( argc <= 1 )
        filename = "../kitten-slices.txt";
    else
        filename = argv[1];
//    std::string filename("../test_slice8_case1.txt");
//    std::string filename("../test_slice8_case2.txt");
//    std::string filename("../test_slice_square.txt");
//    std::string filename("../test_slice_autointersect.txt");

    slices = loadSlices(filename, n_slices);
    std::cout << n_slices << std::endl;

    std::vector<int> vec_nb_points;
    vec_nb_points.reserve(slices.size());

    auto start = std::chrono::high_resolution_clock::now();

    for(auto &slice : slices) {
        auto medial = GraphMaker::extractMedialAxis(slice,60);
        vec_nb_points.push_back(medial->numVertex());
    }

    auto end = std::chrono::high_resolution_clock::now();

    int elapsed = std::chrono::duration_cast<std::chrono::microseconds>
    (end-start).count();

    std::cout << "Processing time (microseconds): " << elapsed << std::endl;

    glutInit(&argc, argv);
    glutInitWindowSize(SIZE_X, SIZE_Y);
    win = glutCreateWindow("Display slice medial axis");
    glutDisplayFunc(display);
    glutKeyboardFunc(processNormalKeys);
    glScalef(0.05f, 0.05f, 1.0f);
    glutMainLoop();
}
