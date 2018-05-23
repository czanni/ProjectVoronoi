#ifndef RENDERAREA_H
#define RENDERAREA_H

#include <QBrush>
#include <QPen>
#include <QWidget>
#include <vector>
#include <Graph.h>
#include <set>

//! [0]
class RenderArea : public QWidget
{
    Q_OBJECT

public:
    RenderArea(QWidget *parent = 0);

    QSize minimumSizeHint() const Q_DECL_OVERRIDE;
    QSize sizeHint() const Q_DECL_OVERRIDE;
    virtual void mouseDoubleClickEvent(QMouseEvent *event);

public slots:
    void setPoints(std::vector <GEO::vec2> &points, bool externe);
    void resetPoints();
    void previousPoints();
    void savePoints();
    void loadPoints();
    void setExterieurContoursBool(bool exterieurContoursBool);
    void linkPoints();
    void drawGraph(Graph &graph, QPainter &painter);
    void drawGraphV(Graph &graph, QPainter &painter);
    void voronoiDiagram();
    void nxtSlice();
    void prevSlice();

protected:
    void paintEvent(QPaintEvent *event) Q_DECL_OVERRIDE;

private:
    unsigned int slice=0;
    Graph ContoursExterieur;
    Graph ContoursInterieur;
    Graph VoronoiExterieur;
    Graph VoronoiInterieur;
    std::vector<Graph> slices;
    std::set <std::pair<int,int>> edgeIntersects;


    bool IsExterieurContours;
};


#endif // RENDERAREA_H
