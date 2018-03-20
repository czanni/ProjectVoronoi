#ifndef WINDOW_H
#define WINDOW_H

#include <QWidget>

QT_BEGIN_NAMESPACE
class QPushButton;
class QCheckBox;
QT_END_NAMESPACE
class RenderArea;

class Window : public QWidget
{
    Q_OBJECT

public:
    Window();

private slots:
    void saveImage();

private:
    RenderArea *renderArea;
    QPushButton *saveImageButton;
    QPushButton *resetButton;
    QPushButton *previousPointsButton;
    QPushButton *savePointsButton;
    QPushButton *loadPointsButton;
    QCheckBox *exterieurContoursCheckBox;
    QPushButton *linkPointsButton;
    QPushButton *voronoiButton;
};


#endif // WINDOW_H
