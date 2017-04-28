#include "foundation.h"
#include "myclasses.h"

float rfloat()  {
    return ((float) rand()) / (float) RAND_MAX;
}

QVector3D reflect(QVector3D d, QVector3D n)  {
    return d.normalized() - 2 * QVector3D::dotProduct(d.normalized(), n) * n;
}

bool refract(QVector3D d, QVector3D n, float index_t, QVector3D *t)  {
    d.normalize();
    n.normalize();
    //float index = 1.0f;
    qreal c = QVector3D::dotProduct(d, n);
    float discr = 1.0f - (1.0f - c*c) / (index_t * index_t);
    if (discr < 0)  return false;

    QVector3D temp = (1.0/index_t) * (d - c * n);
    temp = temp - n * sqrt(discr);
    //cerr << "t=" << temp.x() << "," << temp.y() << "," << temp.z() << " ";
    *t = temp;  // copy refracted vector
    return true;
}

RawRgb::RawRgb()  { r=0; g=0; b=0; }
//RawRgb::RawRgb(const QColor& q)  { r=q.redF(); g=q.greenF(); b=q.blueF(); }

void RawRgb::clamp()  {
    r = min(1.0, r);
    g = min(1.0, g);
    b = min(1.0, b);
}
RawRgb RawRgb::operator+ (const RawRgb& param) {
    RawRgb temp;
    temp.r = r + param.r;
    temp.g = g + param.g;
    temp.b = b + param.b;
    return temp;
}
RawRgb RawRgb::operator/ (const double param)  {
    RawRgb temp;
    temp.r = r / param;
    temp.g = g / param;
    temp.b = b / param;
    return temp;
}
RawRgb RawRgb::operator* (const RawRgb& param)  {
    RawRgb temp;
    temp.r = r * param.r;
    temp.g = g * param.g;
    temp.b = b * param.b;
    return temp;
}
RawRgb RawRgb::operator* (const QColor& param)  {
    RawRgb temp;
    temp.r = r * param.redF();
    temp.g = g * param.greenF();
    temp.b = b * param.blueF();
    return temp;
}

RawRgb RawRgb::operator* (const float& param) const {
    RawRgb temp;
    temp.r = r * param;
    temp.g = g * param;
    temp.b = b * param;
    return temp;
}

RawRgb operator* (const float& param1, const RawRgb& param2) {
    return param2 * param1;
}

QColor RawRgb::toQColor()  {
    return QColor::fromRgbF(min(1.0, r), min(1.0, g), min(1.0, b));
}


Material::Material(QColor d, QColor s, int n, float r, float index) : surface(d), specular(s),
    specpower(n), reflect(r), refract_n(index) { }

void SurfaceData::setboundingsphere(QVector3D p, qreal r)  {
    bounds.centre = p;
    bounds.radius = r;
}

bool SurfaceData::checkboundingbox(QVector3D e, QVector3D d)  {
    return (bounds.centre.distanceToLine(e, d.normalized()) <= bounds.radius);
}

SphereData::SphereData(double r) : radius(r) {
    setboundingsphere(QVector3D(), r);
}
bool SphereData::calchit(QVector3D e, QVector3D d, qreal t0, qreal t1, Hit *h)  {
    bool hit = false;
    qreal A = d.lengthSquared();
    // using instancing, position is always (0,0,0)
    QVector3D e_min_c = e;
    qreal B = QVector3D::dotProduct(d, e_min_c);
    qreal C = e_min_c.lengthSquared() - radius*radius;

    qreal discr = B*B - A*C;
    if (discr == 0)  {
        qreal t = -B / A;
        h->t = t;
        h->n = (e + h->t * d).normalized();
        hit = true;
    }
    else if (discr > 0)  {
        //qreal raydist = (d-e).length();
        qreal raydist = 0.001; // epsilon? for debugging shadows

        qreal t2 = (-B - (qreal)sqrt((double)discr)) / A;
        qreal t3 = (-B + (qreal)sqrt((double)discr)) / A;
        if (t2 > raydist) {
            h->t = t2;
            h->n = (e + h->t * d).normalized();
            hit = true;
        }
        else if (t3 > raydist)  {
            h->t = t3;
            h->n = (e + h->t * d).normalized();
            hit = true;
        }
    }
    return hit;
}

MeshData::MeshData(std::string fname)  {
    std::ifstream myfile(fname.c_str());
    if (myfile.is_open())  {
        std::string s;
        while (myfile >> s)  {
            if (s.compare("v") == 0)  {
                float x, y, z;
                myfile >> x >> y >> z;
                QVector3D v(x, y, z);
                vertices.push_back(v);
            }
            else if (s.compare("vn") == 0)  {
                float x, y, z;
                myfile >> x >> y >> z;
                QVector3D vn(x, y, z);
                normals.push_back(vn);
            }
            else if (s.compare("f") == 0)  {
                unsigned *f = new unsigned[3];
                myfile >> f[0] >> f[1] >> f[2];
                --f[0]; --f[1]; --f[2];
                faces.push_back(f);
            }
            else  {
                cerr << "Error reading word: " << s << "\n";
                continue;
            }
        }
        myfile.close();
    }

    // if vector normals provided, use smooth shading
    smoothshading = !(normals.empty());

    // make bounding sphere (assume there's at least 1 face)
    if (!faces.empty())  {
        cerr << "making bounding sphere for mesh\n";
        // simple algorithm adapted from http://www.mvps.org/directx/articles/using_bounding_spheres.htm
        QVector3D centre;
        for(unsigned int i=0; i<vertices.size(); ++i)  {
            centre += vertices[i];
        }
        centre /= (qreal)vertices.size();

        QVector3D j = vertices[0];
        qreal distsq = (centre - j).lengthSquared();
        for(unsigned int i=1; i<vertices.size(); ++i)  {
            qreal distsq1 = (centre - vertices[i]).lengthSquared();
            if (distsq1 > distsq)
                distsq = distsq1;
        }
        setboundingsphere(centre, (qreal)sqrt(distsq));
    }
}
MeshData::~MeshData()  {
    vertices.clear();
    normals.clear();
    while (!faces.empty())  delete faces.front(), faces.pop_front();
    faces.clear();
}

bool MeshData::calchit(QVector3D e, QVector3D d, qreal t0, qreal t1, Hit *h)  {
    bool hit = false;

    std::list<unsigned int*>::iterator i;
    for (i=faces.begin(); i!=faces.end(); ++i)  {
        unsigned *face = *i;

        QVector3D a = vertices[face[0]];
        QVector3D b = vertices[face[1]];
        QVector3D c = vertices[face[2]];
        QVector3D B = a - b;
        QVector3D C = a - c;
        QVector3D E = a - e;
        qreal ei_hf = C.y()*d.z() - d.y()*C.z();
        qreal gf_di = d.x()*C.z() - C.x()*d.z();
        qreal dh_eg = C.x()*d.y() - d.x()*C.y();
        qreal ak_jb = B.x()*E.y() - E.x()*B.y();
        qreal jc_al = E.x()*B.z() - B.x()*E.z();
        qreal bl_kc = B.y()*E.z() - E.y()*B.z();

        qreal M = B.x()*ei_hf + B.y()*gf_di + B.z()*dh_eg;

        qreal t = -(C.z()*ak_jb + C.y()*jc_al + C.x()*bl_kc) / M;
        if ((t < t0)||(t > t1))
            continue;
        qreal gamma = (d.z()*ak_jb + d.y()*jc_al + d.x()*bl_kc) / M;
        if ((gamma < 0)||(gamma > 1))
            continue;
        qreal beta = (E.x()*ei_hf + E.y()*gf_di + E.z()*dh_eg) / M;
        if ((beta < 0)||(beta > (1-gamma)))
            continue;

        // we have a hit, but not sure if it's the closest one
        hit = true;
        h->t = t;
        // if smoothshading (normals stored), interpolate normals
        if (smoothshading)  {
            QVector3D na, nb, nc;  // the vertex normals provided
            na = normals[face[0]];
            nb = normals[face[1]];
            nc = normals[face[2]]; // remember p = a + beta(b-a) + gamma(c-a)
            h->n = na + beta*(nb-na) + gamma*(nc-na);  // Phong shading - interpolate normals
        }
        else
            h->n = QVector3D::normal(a, b, c);  // not BxC - use original vertex positions
        t1 = t;

    }
    return hit;
}

MySurface::MySurface(QMatrix4x4 t) : transform(t) { }
QVector3D MySurface::transform_vector(QVector3D a)  { return transform.inverted().mapVector(a); }
QVector3D MySurface::transform_point(QVector3D p)  { return transform.inverted().map(p); }
QVector3D MySurface::transform_normal(QVector3D n)  { return transform.inverted().transposed().mapVector(n); }

void MyGroup::add(MySurface *s) { members.push_back(s); }
MyGroup::MyGroup() : MySurface(QMatrix4x4()) { }
MyGroup::MyGroup(QMatrix4x4 t) : MySurface(t) { }
MyGroup::~MyGroup()
{
    while (!members.empty())  delete members.front(), members.pop_front();
    members.clear();
}

bool MyGroup::hits(QVector3D e, QVector3D d, qreal t0, qreal t1, Hit *h)  {
    bool hit = false;
    Hit hitobject;

    std::list<MySurface*>::iterator o;
    for (o=members.begin(); o!=members.end(); ++o)
        if ((*o)->hits(transform_point(e), transform_vector(d), t0, t1, &hitobject)) {
            if ((hitobject.t < t0) || (hitobject.t > t1)) continue;
            hit = true;

            // copy hitobject into h
            h->t = hitobject.t;
            h->n = transform_normal(hitobject.n);
            h->m = hitobject.m;

            t1 = h->t;
        }

    return hit;
}


bool MyInstance::hits(QVector3D e, QVector3D d, qreal t0, qreal t1, Hit *h)  {
    if (data->checkboundingbox(transform_point(e), transform_vector(d)))
        if (data->calchit(transform_point(e), transform_vector(d), t0, t1, h))  {
            h->m = material;
            h->n = transform_normal(h->n);
            return true;
        }
    return false;
}

MyInstance::MyInstance(QMatrix4x4 t, SurfaceData *d, Material *m) : MySurface(t), data(d), material(m) { }


MyLightSource::MyLightSource(QVector3D p, QColor c) : position(p), colour(c) {
    a0 = a2 = 0;
    a1 = 0.176;  // linear distance attenuation
}

DirLightSource::DirLightSource(QVector3D d, QColor c) : MyLightSource(QVector3D(), c), direction(d) {
    a0 = 1;  // constant attenuation
    a1 = a2 = 0;
}

MyCamera::MyCamera(QVector3D p) : position(p) {
    headup = QVector3D();
    lookat = QVector3D();
}
MyCamera::MyCamera()  { MyCamera(QVector3D()); }


RawRgb MyScene::raycolor(QVector3D ray_e, QVector3D ray_d, qreal t0, qreal t1, int rec)  {
    RawRgb c;
    if (rec < 0)  return c; // recursion test

    Hit hit;
    if (surfaces.hits(ray_e, ray_d, t0, t1, &hit))  {
        QVector3D p = ray_e + hit.t * ray_d;
        QVector3D N = hit.n.normalized();

        QColor diffuse = hit.m->surface;
        QColor specular = hit.m->specular;
        int specn = hit.m->specpower;

        std::list<MyLightSource*>::iterator light;
        for (light=lights.begin(); light!=lights.end(); ++light)  {
            QColor colour = (*light)->colour;
            RawRgb c_i;

            // translucent shadows - no refraction
            RawRgb k = shadowcolour(*light, p);
            if ((k.r)||(k.g)||(k.b))  {
                QVector3D L = (*light)->L(p);
                L.normalize();

                float max_diffuse = max(0.0f, QVector3D::dotProduct(N, L));

                c_i.r += diffuse.redF() * colour.redF() * max_diffuse;
                c_i.g += diffuse.greenF() * colour.greenF() * max_diffuse;
                c_i.b += diffuse.blueF() * colour.blueF() * max_diffuse;

                if (specular != QColor(0,0,0))  {
                    QVector3D H = L + (-ray_d).normalized();  // L + V (the vector from surface to eye)
                    H.normalize();

                    float max_specular = max(0.0f, (float)pow(QVector3D::dotProduct(N, H), specn));

                    c_i.r += specular.redF() * colour.redF() * max_specular;
                    c_i.g += specular.greenF() * colour.greenF() * max_specular;
                    c_i.b += specular.blueF() * colour.blueF() * max_specular;
                }

                c = c + ((*light)->att(p) * k * c_i);  // multiply by shadow colour
            }
        }
        c.r += diffuse.redF() * ambientlight.redF();
        c.g += diffuse.greenF() * ambientlight.greenF();
        c.b += diffuse.blueF() * ambientlight.blueF();

        // simple translucency - pick a different function?
        c = diffuse.alphaF() * c;

        // simple reflection
        if (hit.m->reflect)  {
            QVector3D r = reflect(ray_d, N);
            // colour = direct + reflected
            c = c + hit.m->reflect * (raycolor(p, r, EPSILON, t1, rec-1) * specular);
        }

        // refraction
        float refract_n = hit.m->refract_n;
        if (refract_n > 0)  {
            qreal crit;
            bool TIR = false;
            QVector3D t;  // transmitted ray
            QVector3D r = reflect(ray_d, N);  // reflected ray
            RawRgb k;  // light-attenuation constant
            qreal temp = QVector3D::dotProduct(ray_d.normalized(), N);
            if (temp < 0)  { // outside the material
                refract(ray_d, N, refract_n, &t);
                crit = -temp;
                //cerr << "d.n=" << temp << " ";
                k.r = k.g = k.b = 1;
            }
            else  { // inside the material
                RawRgb a;
                // usually a lightened version of the material colour
                QColor temp = diffuse;
                temp.setHsv(diffuse.hue(), 10, 250); // arbitrary value
                a.r = 1-temp.redF();
                a.g = 1-temp.greenF();
                a.b = 1-temp.blueF();

                //qreal density = 0.4;  // stored in MyScene
                qreal att = density * hit.t;
                k.r = exp(-a.r * att);
                k.g = exp(-a.g * att);
                k.b = exp(-a.b * att);

                // inside object = change sign of normal vector
                if (refract(ray_d, -N, 1.0f/refract_n, &t))  {
                    crit = QVector3D::dotProduct(t.normalized(), N);
                }
                else  {
                    TIR = true;
                }
            }

            // attenuate for partially-opaque materials (e.g. plastic)
            k = (1.0f - diffuse.alphaF()) * k;

            if (TIR)  {
                // total internal reflection - R=0
                c = c + k * raycolor(p, r, EPSILON, t1, rec-1);
            }
            else  {
                // return recursive calls with refracted & reflected rays
                float R0 = (refract_n-1)*(refract_n-1) / (refract_n+1)*(refract_n+1);
                float R = R0 + (1 - R0)*pow(1-crit, 5); // approximation of Fresnel equations

                c = c + k * (R * raycolor(p, r, EPSILON, t1, rec-1) + (1-R) * raycolor(p, t, EPSILON, t1, rec-1));
            }
        }

        c.clamp();
    }
    return c;
}

RawRgb MyScene::shadowcolour(MyLightSource *light, QVector3D p)  {
    RawRgb k;
    k.r = k.g = k.b = 1;
    QVector3D e = p;
    QVector3D d = light->L(p);
    qreal t0, t1;
    t0 = 0;
    int inside = 0;
    const float MINRGB = 0.0039f; // about 1/256
    Hit hit;
    QColor lastmaterial;

    while ((k.r > MINRGB)||(k.g > MINRGB)||(k.b > MINRGB))  {
        if (surfaces.hits(e, d, t0+EPSILON, 1, &hit))  {
            if (hit.m->surface.alphaF() == 1.0)  {
                k.r = k.g = k.b = 0;
                return k;
            }
            t1 = hit.t;
            // check if hit outside
            if (QVector3D::dotProduct(d.normalized(), hit.n.normalized()) < 0)  {
                if (inside)
                    k = k * calculate_k(t0, t1, lastmaterial);
                inside++;
            }
            else  {
                k = k * calculate_k(t0, t1, hit.m->surface);
                inside--;
            }
            lastmaterial = hit.m->surface;
            t0 = t1;
        }
        else  {
            if (inside)
                k = k * calculate_k(t0, 1, lastmaterial);
            return k;
        }
    }
    return k;
}

RawRgb MyScene::calculate_k(qreal t0, qreal t1, QColor material)  {
    RawRgb k;
    float att = density * (t1-t0);
    k.r = exp((material.redF()-1) * att);
    k.g = exp((material.greenF()-1) * att);
    k.b = exp((material.blueF()-1) * att);
    k = (1.0-material.alphaF()) * k;
    return k;
}
