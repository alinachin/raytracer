#ifndef MYCLASSES_H
#define MYCLASSES_H

#include "foundation.h"
#include <list>
#include <vector>
#include <cctype>
#include <QVector3D>
#include <QMatrix4x4>
#include <QColor>

#define BIGNUMBER 1000

float rfloat();
QVector3D reflect(QVector3D d, QVector3D n);
bool refract(QVector3D d, QVector3D n, float index_t, QVector3D *t);

struct RawRgb  {  // unbounded floats >= 0
    qreal r, g, b;
public:
    RawRgb();
    //RawRgb(const QColor&);
    void clamp();
    RawRgb operator+ (const RawRgb& param);
    RawRgb operator/ (const double param);
    RawRgb operator* (const RawRgb&);
    RawRgb operator* (const QColor& param);
    RawRgb operator* (const float&) const;
    QColor toQColor();
};
RawRgb operator* (const float&, const RawRgb&);

struct BoundingSphere {
    QVector3D centre;
    qreal radius;
};

struct Material {
    QColor surface; // diffuse colour + alpha
    QColor specular; // specular colour
    int specpower; // specular power
    // new:
    float reflect; // scale down reflection
    float refract_n; // index of refraction >1

    Material(QColor d, QColor s, int n, float ref, float index);
};


struct Hit {
    qreal t;
    QVector3D n;
    Material *m;
};

class SurfaceData {
    BoundingSphere bounds;
protected:
    void setboundingsphere(QVector3D p, qreal r);
public:
    bool checkboundingbox(QVector3D e, QVector3D d);
    virtual bool calchit(QVector3D, QVector3D, qreal, qreal, Hit*) { return false; }
};

class SphereData : public SurfaceData {
    double radius;
public:
    SphereData(double r);
    bool calchit(QVector3D e, QVector3D d, qreal t0, qreal t1, Hit *h);
};
class MeshData : public SurfaceData {
    std::vector<QVector3D> vertices;
    std::vector<QVector3D> normals;
    std::list<unsigned int*> faces;
    bool smoothshading;
public:
    MeshData(std::string fname);
    ~MeshData();
    bool calchit(QVector3D e, QVector3D d, qreal t0, qreal t1, Hit *h);
};


class MySurface {
    QMatrix4x4 transform;
protected:
    MySurface(QMatrix4x4 t);
    QVector3D transform_vector(QVector3D a);
    QVector3D transform_point(QVector3D p);
    QVector3D transform_normal(QVector3D n);
public:
    virtual bool hits(QVector3D, QVector3D, qreal, qreal, Hit*) { return false; }
};

class MyGroup : public MySurface {
    std::list<MySurface*> members;

public:
    void add(MySurface *s);
    MyGroup();
    MyGroup(QMatrix4x4 t);
    ~MyGroup();
    bool hits(QVector3D e, QVector3D d, qreal t0, qreal t1, Hit *h);
};


class MyInstance : public MySurface {
    SurfaceData *data;
    Material *material;
public:
    bool hits(QVector3D e, QVector3D d, qreal t0, qreal t1, Hit *h);
    MyInstance(QMatrix4x4 t, SurfaceData *d, Material *m);
};

// point light source
struct MyLightSource {
    QVector3D position;
    QColor colour;
    qreal a0, a1, a2;

    MyLightSource(QVector3D p, QColor c);
    virtual QVector3D L(QVector3D p)  { return (position - p); }
    qreal att(QVector3D p)  { return 1.0/(a0 + a1*L(p).length() + a2*L(p).lengthSquared()); }
};

struct DirLightSource : public MyLightSource {
    QVector3D direction;

    DirLightSource(QVector3D d, QColor c);
    QVector3D L(QVector3D p)  { return (-100 * direction); }
};

struct MyCamera {
    QVector3D position;
    QVector3D headup;
    QVector3D lookat;

    MyCamera(QVector3D p);
    MyCamera();
};



struct MyScene {
    QColor ambientlight;
    MyCamera camera;
    std::list<MyLightSource*> lights;
    MyGroup surfaces;
    float density;

    RawRgb raycolor(QVector3D ray_e, QVector3D ray_d, qreal t0, qreal t1, int rec);
    RawRgb shadowcolour(MyLightSource* light, QVector3D point);
    RawRgb calculate_k(qreal, qreal, QColor);
};

#endif // MYCLASSES_H

