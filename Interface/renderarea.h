#ifndef RENDERAREA_H
#define RENDERAREA_H

#include <QBrush>
#include <QPen>
#include <QWidget>
#include <vector>
#include <Graph.h>

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

protected:
    void paintEvent(QPaintEvent *event) Q_DECL_OVERRIDE;

private:
    Graph ContoursExterieur;
    Graph ContoursInterieur;
    Graph VoronoiExterieur;
    Graph VoronoiInterieur;

    bool IsExterieurContours;
};


#endif // RENDERAREA_H
