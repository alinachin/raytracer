//-------------------------------------------------------------------------------------------
//   Painting with Flowsnakes
// fsnake program modified to use open gl vertex arrays  Brian Wyvill October 2012
// added save/restore and postscript driver November 2012
// fixed memory management November 2012 delete works properly now.
// added progress bar to show memory usage.
//-------------------------------------------------------------------------------------------

#include "glwidget.h"
#include "myclasses.h"

GLWidget::GLWidget(QWidget *parent)
    : QGLWidget(parent)
{

}

GLWidget::~GLWidget()
{    

}

void GLWidget::clear()
{
     updateGL();
}

void GLWidget::initializeGL()
{
    //Background color will be white
    glClearColor(1.0, 1.0, 1.0, 1.0);
    glShadeModel( GL_FLAT );
    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();
    glPointSize(5);

    // todo: get from UI?
    antialias = false;
    aa_samples = 3;
    ray_bounces = 5;
    mat_density = 0.4;
}

void GLWidget::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    displayImage();
}

/* 2D */
void GLWidget::resizeGL( int w, int h )
{
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    glOrtho(0.0,GLdouble(w),0,GLdouble(h),-10.0,10.0);
    glFlush();
    glMatrixMode(GL_MODELVIEW);
    glViewport( 0, 0, (GLint)w, (GLint)h );
    cerr << "gl new size "<< w SEP h NL;
    renderWidth = w;
    renderHeight = h;
}

// no mouse events in this demo
void GLWidget::mousePressEvent( QMouseEvent * )
{
}

void GLWidget::mouseReleaseEvent( QMouseEvent *)
{
}

void GLWidget::mouseMoveEvent ( QMouseEvent * )
{
}

// wheel event
void GLWidget::wheelEvent(QWheelEvent *)
{
}

void GLWidget::openImage(QString fileBuf)
{     
    QImage myImage;
    myImage.load(fileBuf);
    prepareImageDisplay(&myImage);
}

void GLWidget::prepareImageDisplay(QImage* myimage)
{   
    glimage = QGLWidget::convertToGLFormat( *myimage );  // flipped 32bit RGBA stored as mi
    updateGL();    
}

void GLWidget::displayImage()
{
    if (glimage.width()==0) {
        cerr << "Null Image\n";
        return;
    } else {
        glRasterPos2i(0,0);
        glDrawPixels(glimage.width(), glimage.height(), GL_RGBA, GL_UNSIGNED_BYTE, glimage.bits());
        glFlush();
    }
}

void GLWidget::saveImage( QString fileBuf )
{
    // there is no conversion  back toQImage
    // hence the need to keep qtimage as well as glimage
    qtimage.save ( fileBuf );   // note it is the qtimage in the right format for saving
}

void GLWidget::renderImage()
{
    cerr << "render image from memory\n";

    QImage myimage(renderWidth, renderHeight, QImage::Format_RGB32);

    MyCamera arbitrarycamera; // TODO
    scene.density = mat_density;

    cerr << "beginning ray trace\n";
    qreal nx, ny;
    nx = (double)myimage.width();
    ny = (double)myimage.height();
    qreal b, t, l, r;
    b = -10;
    t = 10;
    l = nx / ny * b;
    r = -l;
    qreal d = -scene.camera.position.z();  // distance from camera to viewing plane

    int n;  // antialias sample rate (n squared)

    if (!antialias)  n = 1;
    else             n = aa_samples;

    int N = n*n;

    for (int j=0; j<myimage.height(); ++j)  {
        for (int i=0; i<myimage.width(); ++i) {
            // antialiasing
            RawRgb pc; // pixel colour (average of n-squared samples)
            for (int p=0; p<n; ++p)
                for (int q=0; q<n; ++q)  {
                    // make a ray
                    qreal u, v;
                    if (antialias)  {
                        u = l + (r-l)*(i+(p+rfloat())/n) / nx;
                        v = t - (t-b)*(j+(q+rfloat())/n) / ny;
                    }
                    else  {
                        u = l + (r-l)*(i+0.5) / nx;
                        v = t - (t-b)*(j+0.5) / ny;
                    }

                    QVector3D ray_e = scene.camera.position;
                    QVector3D ray_d(u, v, d);
                    //ray_d.normalize(); // don't normalize? (textbook 4.7)

                    // determine if the ray hits anything (if so, get surface normal)
                    qreal t0 = 0;
                    qreal t1 = BIGNUMBER;

                    /* start raycolor() */
                    pc = pc + scene.raycolor(ray_e, ray_d, t0, t1, ray_bounces);
                    /* end raycolor() */
                }

            pc = pc / N;
            pc.clamp();
            myimage.setPixel(i, j, pc.toQColor().rgb());
        }
        if (!(j%16))  {
            prepareImageDisplay(&myimage);
        }
    }
    qtimage=myimage.copy(0, 0,  myimage.width(), myimage.height()); // make copy for saving to file

    prepareImageDisplay(&myimage);
}

void GLWidget::makeImage(QString fileBuf)
{   
   std::ifstream myfile(fileBuf.toUtf8());
   if (!myfile.is_open())  {
       cerr << "unable to open file\n" << "previous scene data still saved\n";
       return;
   }

    cerr << "clearing material & data tables\n";
    for(map<string,Material*>::iterator it = material_table.begin(); it!=material_table.end(); ++it)
        delete it->second;
    for(map<string,SphereData*>::iterator it = sphere_table.begin(); it!=sphere_table.end(); ++it)
        delete it->second;
    for(map<string,MeshData*>::iterator it = mesh_table.begin(); it!=mesh_table.end(); ++it)
        delete it->second;
    material_table.clear();
    sphere_table.clear();
    mesh_table.clear();

    cerr << "clearing scene\n";
    for(list<MyLightSource*>::iterator it = scene.lights.begin(); it!=scene.lights.end(); ++it)
        delete *it;
    scene.lights.clear();
    scene = MyScene();

    list<MyGroup*> currentgroup;
    currentgroup.push_front(&(scene.surfaces));

    string s;
    while (myfile >> s)  {
        if (s.compare("SPHERES") == 0)  {
            string v;
            myfile >> v;
            if (v.compare("(") == 0)  {  // get list of spheres
                while ((myfile >> v)&&(v.compare(")") != 0)) {
                    double rad;
                    myfile >> rad;
                    sphere_table[v] = new SphereData(rad);
                    cerr << "added new sphere " << v << "\n";
                }
            }
        }
        else if (s.compare("MESHES") == 0)  {
            string v;
            myfile >> v;
            if (v.compare("(") == 0)  {  // get list of meshes
                while ((myfile >> v)&&(v.compare(")") != 0)) {
                    string objfile;
                    myfile >> objfile;
                    string currfile = fileBuf.toStdString();
                    string path = currfile.substr(0, currfile.find_last_of("/\\"));
                    objfile = path + "/" + objfile;
                    ifstream test(objfile.c_str());
                    if (test.is_open())
                        test.close();
                    else  {
                        cerr << "Unable to open " << objfile << "\n";
                        continue;
                    }
                    mesh_table[v] = new MeshData(objfile);
                }
            }
        }
        else if (s.compare("MATERIALS") == 0)  {
            /* invariants:
             * if refraction is enabled, specular reflection is disabled (handled in refraction code)
             * translucent objects (alpha != 1.0f) have a refractive index = 1.55 (arbitrary)
             * opaque objects are not refractive
             */
            string v;
            myfile >> v;
            if (v.compare("(") == 0)  {
                while ((myfile >> v)&&(v.compare(")") != 0)) {
                    //string varname = v;
                    float r, g, b, f;
                    myfile >> r >> g >> b >> f;
                    QColor surface = QColor::fromRgbF(r, g, b, f);
                    myfile >> r >> g >> b;
                    QColor spec = QColor::fromRgbF(r, g, b);
                    int specpow;
                    myfile >> specpow;

                    // new:
                    float reflect;
                    myfile >> reflect;
                    float refract; // 0 = no reflection, 1 = perfect reflection (scales to specular)
                    myfile >> refract; // index of refraction - set to 0 to disable refraction

                    // check invariants
                    if (surface.alphaF() == 1.0f)  refract = 0.0f;
                    else
                        if (refract == 0.0f)  refract = 1.55f; // sets refraction index if it was disabled
                    if (refract > 0.0f)  reflect = 0.0f;

                    material_table[v] = new Material(surface, spec, specpow, reflect, refract);
                }
            }
        }

        else if (s.compare("ambient") == 0)  {
            float r, g, b;
            myfile >> r >> g >> b;
            scene.ambientlight = QColor::fromRgbF(r, g, b);
        }
        else if (s.compare("camera") == 0)  {
            float x, y, z;
            myfile >> x >> y >> z;
            scene.camera = MyCamera(QVector3D(x, y, z));
        }
        else if (s.compare("lightsrc") == 0)  {
            string t;
            myfile >> t;
            float x, y, z;
            myfile >> x >> y >> z;
            QVector3D v(x, y, z);
            float r, g, b;
            myfile >> r >> g >> b;
            QColor colour = QColor::fromRgbF(r, g, b);
            if (t.compare("point")==0)
                scene.lights.push_back(new MyLightSource(v, colour));
            else
                scene.lights.push_back(new DirLightSource(v, colour));
        }

        else if (s.compare("{") == 0)  {  // start of a group
            QMatrix4x4 transform;
            // read transform(s)
            string v;
            myfile >> v;
            if (v.compare("(") == 0)   // list of transforms
                while ((myfile >> v)&&(v.compare(")") != 0))  {
                    if (v.compare("trans") == 0)  {
                        float x, y, z;
                        myfile >> x >> y >> z;
                        transform.translate(x, y, z);
                        cerr << "translate\n";
                    }
                    else if (v.compare("rot") == 0)  {
                        float degx, degy, degz;
                        myfile >> degx >> degy >> degz ;
                        if (degx != 0)
                            transform.rotate(degx, 1, 0, 0);
                        if (degy != 0)
                            transform.rotate(degy, 0, 1, 0);
                        if (degz != 0)
                            transform.rotate(degz, 0, 0, 1);
                        cerr << "rotate\n";
                    }
                    else if (v.compare("scale") == 0)  {
                        float x, y, z;
                        myfile >> x >> y >> z;
                        transform.scale(x, y, z);
                        cerr << "scale\n";
                    }
                }
            // make group & push onto stack & **add under current group**
            MyGroup *group = new MyGroup(transform);
            currentgroup.front()->add(group);
            currentgroup.push_front(group);
            cerr << "pushed a new group\n";
        }
        else if (s.compare("}") == 0)  {  // end of a group
            currentgroup.pop_front();
            // check if overall group has been popped?
            cerr << "popped a group\n";
        }

        else if ((s.compare("sphere") == 0)||(s.compare("mesh") == 0))  {
            QMatrix4x4 transform;
            SurfaceData *data;
            Material *material;

            string v;
            myfile >> v;
            if (v.compare("(") == 0)   // list of transforms
                while ((myfile >> v)&&(v.compare(")") != 0))  {
                    if (v.compare("trans") == 0)  {
                        float x, y, z;
                        myfile >> x >> y >> z;
                        transform.translate(x, y, z);
                    }
                    else if (v.compare("rot") == 0)  {
                        float degx, degy, degz;
                        myfile >> degx >> degy >> degz ;
                        if (degx != 0)
                            transform.rotate(degx, 1, 0, 0);
                        if (degy != 0)
                            transform.rotate(degy, 0, 1, 0);
                        if (degz != 0)
                            transform.rotate(degz, 0, 0, 1);
                    }
                    else if (v.compare("scale") == 0)  {
                        float x, y, z;
                        myfile >> x >> y >> z;
                        transform.scale(x, y, z);
                    }
                }
            myfile >> v;
            if (s.compare("sphere") == 0)
                if (sphere_table.count(v))  {
                    data = sphere_table[v];
                }
                else  {
                    cerr << "unable to find value assigned to " << v << "\n";
                    continue;
                }
            else
                if (mesh_table.count(v))  {
                    data = mesh_table[v];
                }
                else  {
                    cerr << "unable to find value assigned to " << v << "\n";
                    continue;
                }
            myfile >> v;
            if (material_table.count(v))  {
                material = material_table[v];
            }
            else  {
                cerr << "unable to find value assigned to " << v << "\n";
                continue;
            }
            MySurface *sphere = new MyInstance(transform, data, material);
            currentgroup.front()->add(sphere);
        }
    }
    myfile.close();

    renderImage();

}

void GLWidget::about()
{
    QString vnum;
    QString mess, notes;
    QString title="Images in Qt and Opengl ";

    vnum.setNum ( MYVERSION );
    mess="Qt OpenGl image demo Release Version: ";
    mess = mess+vnum;
    notes = "\n\n News: Every QImage is now on stack, there is no way for memory leak. -- Lucky";
    mess = mess+notes;
    QMessageBox::information( this, title, mess, QMessageBox::Ok );
}

void GLWidget::help()
{
    QString vnum;
    QString mess, notes;
    QString title="qtglimages";

    vnum.setNum ( MYVERSION);
    mess="Simple Image Handling in Qt/Opengl by Brian Wyvill Release Version: ";
    mess = mess+vnum;
    notes = "\n\n Save and Load images for display.  Also Make an image such as output from a ray tracer. ";
    mess = mess+notes;
    QMessageBox::information( this, title, mess, QMessageBox::Ok );
}

