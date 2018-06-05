#include "renderarea.h"
#include "polygoninout.h"
#include <graphmaker.h>

#include <QPainter>
#include <QMouseEvent>
#include <QFileDialog>
#include <iostream>


RenderArea::RenderArea(QWidget *parent)
    : QWidget(parent)
{
    ContoursExterieur = Graph();
    ContoursInterieur = Graph();

    VoronoiExterieur = Graph();
    VoronoiInterieur = Graph();

    IsExterieurContours = true;

    setBackgroundRole(QPalette::Base);
    setAutoFillBackground(true);
}

QSize RenderArea::minimumSizeHint() const //Taille minimale de la fenêtre de l'interface
{
    return QSize(100, 100);
}

QSize RenderArea::sizeHint() const
{
    return QSize(800, 800);
}

void RenderArea::mouseDoubleClickEvent(QMouseEvent *event) //Si double click, alors on ajoute un point dans le graphe considéré
{
    QPoint point = event->pos();
    if (IsExterieurContours){
 //       GEO::vec2 test = GEO::vec2((double)<point.x(), (double)point.y());
        this->ContoursExterieur.addPoint(GEO::vec2((double)point.x(), (double)point.y()));
    }
    else{
        this->ContoursInterieur.addPoint(GEO::vec2((double)point.x(), (double)point.y()));
    }
    update();
}

//void RenderArea::previousPoints()
//{
//    if (IsExterieurContours){
//        this->ContoursExterieur.removeLastPoint();
//    }
//    else{
//        this->ContoursInterieur.removeLastPoint();
//    }
//    update();
//}
void RenderArea::previousPoints(){}

//void RenderArea::setPoints(std::vector <GEO::vec2> &points, bool externe)
//{
//    if (externe){
//        foreach(GEO::vec2 point, points){
//            this->ContoursExterieur.addPoint(point);
//        }
//    }
//    else{
//        foreach(GEO::vec2 point, points){
//            this->ContoursInterieur.addPoint(point);
//        }
//    }
//    update();
//}
void RenderArea::setPoints(std::vector <GEO::vec2> &points, bool externe){}

void RenderArea::resetPoints() //Pour enlever tous les points du graphe considéré
{
    if (IsExterieurContours){
        this->ContoursExterieur = Graph();
        this->VoronoiExterieur = Graph();
        }
    else{
        this->ContoursInterieur = Graph();
        this->VoronoiInterieur = Graph();

    }
    update();
}

void RenderArea::savePoints() //Pour enregistrer les graphes dans un fichier texte
{
    QString filename = QFileDialog::getSaveFileName(this, "Enregistrer un fichier texte", QString(), "Texte (*.txt)");//Je saisis sous quel nom je veux enregistrer le fichier
    Medial::savePoints(filename.toStdString(), ContoursExterieur, ContoursInterieur);
}

/*void RenderArea::loadPoints(){ //Pour charger les graphes à partir d'un fichier texte
    QString filename = QFileDialog::getOpenFileName(this, "Charger un fichier texte", QString(), "Texte (*.txt)");
    ContoursExterieur = Medial::loadPoints(filename.toStdString(), true);
    ContoursInterieur = Medial::loadPoints(filename.toStdString(), false);
}*/


void RenderArea::loadPoints(){ //Pour charger les slices à partir d'un fichier texte
    QString filename = QFileDialog::getOpenFileName(this, "Charger un fichier texte", QString(), "Texte (*)");
    slices = Medial::loadSlices(filename.toStdString());
    edgeIntersects = *new std::set <std::pair<int,int>>();
    slice=96;
    ContoursExterieur = slices[slice];
    update();
}

void RenderArea::nxtSlice(){
    if (slice < slices.size()){
        ++slice;
        ContoursExterieur = slices[slice];
        VoronoiExterieur = * GraphMaker::extractMedialAxis(ContoursExterieur);
        update();
    }
    else {
        slice=0;
    }

}
void RenderArea::prevSlice(){
    if (slice < slices.size()){
        --slice;
        ContoursExterieur = slices[slice];
        update();
    }
    else {
        slice=0;
    }
}

void RenderArea::setExterieurContoursBool(bool exterieurContoursBool) //Si la case est coché, on traite le contour extérieur. Sinon, le contour intérieur.
{
    this->IsExterieurContours = exterieurContoursBool;
    update();
}

std::vector <QPointF> ConvertVec2toQPointF(std::vector <GEO::vec2> points){ //Conversion des points de Vec2 à QPointF pour l'affichage dans l'interface
    std::vector <QPointF> pointsF = {};
    foreach(GEO::vec2 point, points){
        pointsF.push_back(QPointF((qreal)point[0], (qreal)point[1]));
    }
    return pointsF;
}

void RenderArea::linkPoints(){ //Pour relier les points du graphe considéré
    if (IsExterieurContours){
        int size = this->ContoursExterieur.numVertex();
        for(int i = 0 ; i < size ; i++){
            std::array<int,2> connexion({i, (i+1) % size});
            this->ContoursExterieur.addEdge(connexion);
        }
    }
    else{
        int size = this->ContoursInterieur.numVertex();
        for(int i = 0 ; i < size ; i++){
            std::array<int,2> connexion({i, (i+1) % size});
            this->ContoursInterieur.addEdge(connexion);
        }
    }
    update();
}

void RenderArea::drawGraph(Graph &graph, QPainter &painter){
    //Conversion des points pour l'affichage
    std::vector <QPointF> pointsF = ConvertVec2toQPointF(graph.getPoints());

    int index = 0;
    foreach(const std::vector<Neighbor>& connexion, graph.getConnexions()){
        foreach(Neighbor i, connexion){
            painter.drawLine(pointsF[i.index], pointsF[index]);
        }
        index++;
    }

    painter.setPen(QPen(Qt::black, 5));
    painter.drawPoints(pointsF.data(), pointsF.size());
}

void RenderArea::drawGraphV(Graph &graph, QPainter &painter){
    //Conversion des points pour l'affichage
    std::vector <QPointF> pointsF = ConvertVec2toQPointF(graph.getPoints());
    std::vector <QPointF> pointsFC = ConvertVec2toQPointF(ContoursExterieur.getPoints());
    int index = 0;


    if (graph.numVertex() > 0) {
        for ( auto pair = graph.getNeighbors().begin(); pair != graph.getNeighbors().end(); ++pair){
            painter.setPen(QPen(QColor(255,127,80), 2));
            painter.drawLine(pointsF[(std::get<0>(pair -> first))]/2 + pointsF[(std::get<1>(pair -> first))]/2, pointsFC[pair->second]);
            painter.setPen(QPen(QColor(0,0,0,2)));

        }
    }

    foreach(const std::vector<Neighbor>& connexion, graph.getConnexions()){
        foreach(Neighbor i, connexion){
            if (edgeIntersects.find(std::pair<int,int>(i.index,index)) != edgeIntersects.end()) {
                    painter.setPen(QPen(Qt::red, 2));
                    painter.drawLine(pointsF[i.index], pointsF[index]);
                    painter.setPen(QPen(Qt::green, 3));

                }
            else {
                painter.drawLine(pointsF[i.index], pointsF[index]);
            }
        }

        index++;
    }
    for (int i=0;i<pointsF.size();++i) {
        if (graph.getStatus(i)==treatment::inside) {
            painter.setPen(QPen(Qt::red, 5));
            painter.drawPoint(pointsF[i]);
            painter.setPen(QPen(Qt::green, 3));
        }
        else if (graph.getStatus(i)==treatment::outside) {
            painter.setPen(QPen(Qt::cyan, 5));
            painter.drawPoint(pointsF[i]);
            painter.setPen(QPen(Qt::green, 3));
        }
        else if (graph.getStatus(i)==treatment::unknown) {
            painter.setPen(QPen(Qt::black, 5));
            painter.drawPoint(pointsF[i]);
            painter.setPen(QPen(Qt::green, 3));
        }
    }

}
/*
void RenderArea::drawGraphV(Graph &graph, QPainter &painter){
    //Conversion des points pour l'affichage
    std::vector <QPointF> pointsF = ConvertVec2toQPointF(graph.getPoints());
    int index = 0;
    //GraphMaker::fixOutsidePoints(graph);

    foreach(const std::vector<int>& connexion, graph.getConnexions()){
        foreach(int i, connexion){
            if (edgeIntersects.find(std::pair<int,int>(i,index)) != edgeIntersects.end()) {
                    painter.setPen(QPen(Qt::red, 2));
                    painter.drawLine(pointsF[i], pointsF[index]);
                    painter.setPen(QPen(Qt::green, 3));

                }
            else {
                painter.drawLine(pointsF[i], pointsF[index]);
            }
        }
        index++;
    }
    for (int i=0;i<pointsF.size();++i) {
        if (graph.getStatus(i)==treatment::inside) {
            painter.setPen(QPen(Qt::red, 5));
            painter.drawPoint(pointsF[i]);
            painter.setPen(QPen(Qt::green, 3));
        }
        else if (graph.getStatus(i)==treatment::outside) {
            painter.setPen(QPen(Qt::cyan, 5));
            painter.drawPoint(pointsF[i]);
            painter.setPen(QPen(Qt::green, 3));
        }
        else if (graph.getStatus(i)==treatment::unknown) {
            painter.setPen(QPen(Qt::black, 5));
            painter.drawPoint(pointsF[i]);
            painter.setPen(QPen(Qt::green, 3));
        }
    }

}*/

void RenderArea::voronoiDiagram(){
    edgeIntersects = std::set <std::pair<int,int>> ();
    if (IsExterieurContours){
        Graph debug = Graph();
        /*debug.addPoint({100,200});
        debug.addPoint({200,100});
        debug.addPoint({300,200});
        debug.addPoint({200,300});
        debug.addPoint({400,100});
        debug.addPoint({500,200});
        debug.addPoint({650,350});
        debug.addPoint({350,350});

        debug.addEdge({0,1});
        debug.addEdge({1,2});
        debug.addEdge({2,3});
        debug.addEdge({3,0});
        debug.addEdge({2,4});
        debug.addEdge({4,5});
        debug.addEdge({5,6});
        debug.addEdge({6,7});
        debug.addEdge({7,2});
        ContoursExterieur = debug;*/
        ContoursExterieur = *GraphMaker::makeMorePoints(ContoursExterieur);

        VoronoiExterieur = * GraphMaker::extractVoronoi(ContoursExterieur, edgeIntersects);
        GraphMaker::fixOutsidePoints(VoronoiExterieur, edgeIntersects);

    }
    else{
        VoronoiInterieur = *GraphMaker::extractVoronoi(ContoursInterieur, edgeIntersects);
    }
    //GraphMaker::fixOutsidePoints(VoronoiExterieur, edgeIntersects);
    update();
}

void RenderArea::paintEvent(QPaintEvent * /* event */)
{
    QPainter painter(this);

    painter.save();

    //Le contour extérieur en bleu
    painter.setPen(QPen(Qt::blue, 3));
    drawGraph(ContoursExterieur, painter);

    //Le diagramme de Voronoi du contours extérieur en vert
    painter.setPen(QPen(Qt::green, 3));
    drawGraphV(VoronoiExterieur, painter);

    //Le contour intérieur en noir
    painter.setPen(QPen(Qt::black, 3));
    drawGraph(ContoursInterieur, painter);


    //Le diagramme de Voronoi du contours intérieur en rouge
    painter.setPen(QPen(Qt::red, 3));
    drawGraph(VoronoiInterieur, painter);

    painter.setPen(QPen(QColor(50,50,50), 2));

    painter.restore();

    painter.setPen(palette().dark().color());
    painter.setBrush(Qt::NoBrush);
}
//! [13]
