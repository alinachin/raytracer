//-------------------------------------------------------------------------------------------
//  University of Victoria Computer Science Department
//	FrameWork for OpenGL application under QT
//  Course title: Computer Graphics CSC305
//-------------------------------------------------------------------------------------------
//These two lines are header guiardians against multiple includes
#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <QGLWidget>
#include <QProgressBar>
#include "foundation.h"
#include <QtGui>
#include <QtOpenGL>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <map>
#include "myclasses.h"

//This is our OpenGL Component we built it on top of QGLWidget
class GLWidget : public QGLWidget
{
    Q_OBJECT

public:
    //Constructor for GLWidget
    GLWidget(QWidget *parent = 0);

    //Destructor for GLWidget
    ~GLWidget();

    void openImage(QString fileBuf);
    void saveImage( QString fileBuf);
    void makeImage(QString fileBuf);
    void renderImage();
    void about();
    void help();

    // other settings
    bool antialias;
    int aa_samples;
    int ray_bounces;
    float mat_density;

protected:
    //Initialize the OpenGL Graphics Engine
    void initializeGL();

    //All our painting stuff are here
    void paintGL();

    //When user resizes main window, the scrollArea will be resized and it will call this function from
    //its attached GLWidget
    void resizeGL(int w, int h);

    //Handle mouse press event in scrollArea
    void mousePressEvent(QMouseEvent * );
    void mouseReleaseEvent(QMouseEvent * );
    //Handle Mouse Move Event
    void mouseMoveEvent(QMouseEvent * );
    void wheelEvent(QWheelEvent * );  // for zoom


private:
    void clear();
    int renderWidth, renderHeight;
    void displayImage();
    QProgressBar* pbar;
    void prepareImageDisplay(QImage* myimage); // converts from Qt to opengl format
    QImage glimage, qtimage;  // paintGL will display the gl formatted image
    // keep the qtimage around for saving (one is a copy of the other

    // parsing files
    map<string,Material*> material_table;
    map<string,SphereData*> sphere_table;
    map<string,MeshData*> mesh_table;

    // scene objects
    MyScene scene;

/*
private slots:
    void callOpenImage();
    void callSaveImage();
    void callOpenTextFile();
    void callRender();
    void setAntialiasing(bool);
    */
};


#endif
