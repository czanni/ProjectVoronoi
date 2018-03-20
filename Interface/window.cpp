#include "renderarea.h"
#include "window.h"

#include <QtWidgets>


const int IdRole = Qt::UserRole;


Window::Window()
{
    renderArea = new RenderArea;
    saveImageButton = new QPushButton("Save Image", this);
    resetButton = new QPushButton("Reset Points", this);
    previousPointsButton = new QPushButton("Previous Points", this);
    savePointsButton = new QPushButton("Save Points", this);
    loadPointsButton = new QPushButton("Load Points", this);
    exterieurContoursCheckBox = new QCheckBox(tr("&Exterieur Contour"));
    linkPointsButton = new QPushButton("Link Points", this);
    voronoiButton = new QPushButton("Create Voronoï Diagram", this);

    connect(saveImageButton, SIGNAL(released()), this, SLOT(saveImage()));
    connect(resetButton, SIGNAL(released()), renderArea, SLOT(resetPoints()));
    connect(previousPointsButton, SIGNAL(released()), renderArea, SLOT(previousPoints()));
    connect(savePointsButton, SIGNAL(released()), renderArea, SLOT(savePoints()));
    connect(loadPointsButton, SIGNAL(released()), renderArea, SLOT(loadPoints()));
    connect(exterieurContoursCheckBox, SIGNAL(toggled(bool)), renderArea, SLOT(setExterieurContoursBool(bool)));
    connect(linkPointsButton, SIGNAL(released()), renderArea, SLOT(linkPoints()));
    connect(voronoiButton, SIGNAL(released()), renderArea, SLOT(voronoiDiagram()));
    connect(voronoiButton, SIGNAL(voronoiDiagram()), this, SLOT(repaint())); //TODO

    QGridLayout *mainLayout = new QGridLayout;
    mainLayout->setColumnStretch(0, 1);
    mainLayout->setColumnStretch(3, 1);
    mainLayout->addWidget(renderArea, 0, 0, 1, 4);
    mainLayout->addWidget(saveImageButton, 2, 0);
    mainLayout->addWidget(resetButton, 2, 1);
    mainLayout->addWidget(previousPointsButton, 2, 2);
    mainLayout->addWidget(savePointsButton, 2, 3);
    mainLayout->addWidget(loadPointsButton, 3, 1);
    mainLayout->addWidget(exterieurContoursCheckBox, 3, 3, 1, 2);
    mainLayout->addWidget(linkPointsButton, 3, 2);
    mainLayout->addWidget(voronoiButton, 4, 1);
    setLayout(mainLayout);

    exterieurContoursCheckBox->setChecked(true);

    setWindowTitle(tr("Interface - Filtrage du diagramme de Voronoi"));
}

void Window::saveImage() //Permet d'enregister l'image
{
    QString path = QFileDialog::getSaveFileName(this, "Enregistrer une image", QString(), "Images (*.png *.gif *.jpg *.jpeg)");//Je saisis sous quel nom je veux enregistrer l'image
    QImage img_save(800, 800, QImage::Format_RGB32); // /!\Taille de l'image
    renderArea->render(&img_save);//je récupère le contenu de renderArea
    img_save.save(path);//J'enregistre le contenu dans le nom que j'ai choisi auparavant /!\ Ne pas oublier de mettre un suffixe (ex : .png)
}
