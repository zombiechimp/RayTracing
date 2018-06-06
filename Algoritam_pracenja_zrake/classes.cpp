//
//  classes.cpp
//  Algoritam_pracenja_zrake
//
//  Created by Tea Jakić on 01/06/2018.
//  Copyright © 2018 Tea Jakić. All rights reserved.
//

#include <stdio.h>
#include <string.h>
#include <istream>
#include <sstream>
#include <fstream>
#include "Vector.cpp"
#include <vector>
#include "Matrix.cpp"


using namespace std;
class Object;

class Face3DInd{
    int v1_;
    int v2_;
    int v3_;
    
public:
    Face3DInd(const int v1, const int v2, const int v3):
    v1_(v1),
    v2_(v2),
    v3_(v3)
    {}
    
    int v1(){
        return v1_;
    }
    
    int v2(){
        return v2_;
    }
    
    int v3(){
        return v3_;
    }
    
    void set_v1(int v1){
        v1_ = v1;
    }
    
    void set_v2(int v2){
        v2_ = v2;
    }
    
    void set_v3(int v3){
        v3_ = v3;
    }
    
    
};



class Vertex3D{
    
public:
    double x_;
    double y_;
    double z_;
    double h_;
    
    Vertex3D(){
        x_ = 0;
        y_ = 0;
        z_ = 0;
        h_ = 1;
        {}
    }
    Vertex3D(const double x, const double y, const double z):
    x_(x),
    y_(y),
    z_(z),
    h_(0)
    {}
    Vertex3D(const double x, const double y, const double z, const double h):
    x_(x),
    y_(y),
    z_(z),
    h_(h)
    {}
    double x() {
        return x_;
    }
    double y() {
        return y_;
    }
    double z(){
        return z_;
    }
    double h(){
        return h_;
    }
    
    void set_x(double x){
        x_ = x;
    }
    
    void set_y(double y){
        y_ = y;
    }
    
    void set_z(double z){
        z_ = z;
    }
    
    void set_h(double h){
        h_ = h;
    }
    
    std::string toString(){
        std::string result;
        
        for(int i = 0; i < 4; i++){
            
        }
        return result;
    }
};



class Colour{
public:
    double r_;
    double g_;
    double b_;
    
    
    Colour(double r, double g, double b):
    r_(r),
    g_(g),
    b_(b)
    {}
    
    Colour(){
        Colour(0,0,0);
    }
    
};


class Source{
public:
    Vertex3D position_;
    Colour colour_;
    
    Source(double x, double y, double z, double r, double g, double b):
        position_(Vertex3D(x,y,z)),
        colour_(Colour(r,g,b))
    {}
    Source(double x, double y, double z):
    position_(Vertex3D(x,y,z)),
    colour_(Colour(1,1,1))
    {}
    
    
    
};

class Coefficients{
    double A_;
    double B_;
    double C_;
    double D_;
    
public:
    Coefficients(const double A, const double B, const double C, const double D):
    A_(A),
    B_(B),
    C_(C),
    D_(D)
    {}
    double A(){
        return A_;
    }
    
    double B(){
        return B_;
    }
    
    double C(){
        return C_;
    }
    
    double D(){
        return D_;
    }
    
    void set_A(double A) {
        A_ = A ;
    }
    
    void set_B(double B) {
        B_ = B;
    }
    
    void set_C(double C){
        C_ = C;
    }
    
    void set_D(double D){
        D_ = D;
    }
    
    
    
};


class Intersection{
    
public:
    Object& object_;
    double lambda_;
    bool front_;
    Vertex3D point_;
    Intersection(Object& object, double lambda, bool front, Vertex3D point):
    object_(object),
    lambda_(lambda),
    front_(front),
    point_(point)
    {}
    
    Intersection(Object &object):
        object_(object),
        lambda_(INFINITY),
        front_(false),
        point_(Vertex3D({0,0,0,0}))
    {}

};



class Object{
public:
    

    
    virtual void updateIntersection(Intersection& inters, Vector3d start, Vector3d d) = 0;
    virtual Vector3d getNormalInPoint(Vertex3D point) = 0;
};




class Orb: public Object{
    Vertex3D center_;
    double radius_;
    Colour ambiance_;
    Colour difusing_;
    Colour reflecting_;
    double n_;
    double kref_;
    
public:
    Orb(const double x, const double y, const double z, const double r, const double ar, const double ag, const double ab, const double dr, const double dg, const double db, const double rr, const double rg, const double rb, const double n, const double kref):
    center_(Vertex3D(x, y, z)),
    radius_(r),
    ambiance_(Colour(ar, ag, ab)),
    difusing_(Colour(dr, dg, db)),
    reflecting_(Colour(rr, rg, rb)),
    n_(n),
    kref_(kref)
    {}
    
    void set_center(Vertex3D center){
        center_ = center;
    }
    
    void set_radius(double radius){
        radius_ = radius;
    }
    
    void set_ambiance(Colour ambiance){
        ambiance_= ambiance;
    }
    
    void set_difusing(Colour difusing){
        difusing_ = difusing;
    }
    
    void set_reflecting(Colour reflecting){
        reflecting_ = reflecting;
    }
    
    void set_n(double n){
        n_ = n;
    }
    
    void set_kref(double kref){
        kref_ = kref;
    }
    
    Vertex3D center(){
        return center_;
    }
    
    double radius(){
        return radius_;
    }
    
    Colour ambiance(){
        return ambiance_;
    }
    
    Colour difusing(){
        return difusing_;
    }
    
    Colour reflecting(){
        return reflecting_;
    }
    
    double n(){
        return n_;
    }
    
    double kref(){
        return kref_;
    }
    
    void updateIntersection(Intersection& inters, Vector3d start, Vector3d d) override{
        Vector3d v_center({center().x(), center().y(), center().z()});
        
//        double a = d.scalar_product(d);
        double b = 2*d.scalar_product(v_center - start);
        double c = (v_center-start).scalar_product(v_center-start) - radius() * radius();
        
        double disc = b*b - 4.0*c;
        
        if(disc < 0) return;
        disc = sqrt(disc);
        double lambda1 = (-b-disc)/(2.0);
        double lambda2 = (-b+disc)/(2.0);
        
        double lambda;
        
        
        if((lambda1 < 0 && lambda2 < 0 )|| lambda1 == 0 || lambda2 == 0) return;
        else if(lambda1 < 0) lambda = lambda2;
        else if(lambda2 < 0) lambda = lambda1;
        else lambda = (lambda1<lambda2) ? lambda1 : lambda2;
       
        
        if(!(lambda < inters.lambda_)) return;
        
        Vertex3D point(start[0] + d[0] * lambda, start[1] + d[1] * lambda, start[2] + d[2] * lambda);
        inters.object_ = *this;
        inters.front_ = true;
        inters.lambda_ = lambda;
        inters.point_ = point;
        
        return;
    }
    
    
    Vector3d getNormalInPoint(Vertex3D point) override{
        Vector3d position({point.x(),point.y(),point.z()});
        Vector3d centerx({center().x(), center().y(), center().z()});
        return (position-centerx).normalize();
    }
    
    
};



class Krpica: public Object{
    Vertex3D center_;
    Vector3d v1_;
    Vector3d v2_;
    double width_;
    double height_;
    Colour ambiance1_;
    Colour ambiance2_;
    Colour difusing1_;
    Colour difusing2_;
    Colour reflecting1_;
    Colour reflecting2_;
    double n1_;
    double n2_;
    double kref1_;
    double kref2_;
    
public:
    
    Krpica(double x, double y, double z, double v1x, double v1y, double v1z, double v2x, double v2y, double v2z, double width, double height, double ar, const double ag, const double ab, const double dr, const double dg, const double db, const double rr, const double rg, const double rb, const double n, const double kref, double ar2, const double ag2, const double ab2, const double dr2, const double dg2, const double db2, const double rr2, const double rg2, const double rb2, const double n2, const double kref2):
    center_(Vertex3D(x, y, z)),
    v1_(Vector3d(std::array<double,3>{v1x, v1y, v1z})),
    v2_(Vector3d(std::array<double,3>{v2x, v2y, v2z})),
    width_(width),
    height_(height),
    ambiance1_(Colour(ar, ag, ab)),
    difusing1_(Colour(dr, dg, db)),
    reflecting1_(Colour(rr, rg, rb)),
    n1_(n),
    kref1_(kref),
    ambiance2_(Colour(ar2, ag2, ab2)),
    difusing2_(Colour(dr2, dg2, db2)),
    reflecting2_(Colour(rr2, rg2, rb2)),
    n2_(n2),
    kref2_(kref2)
    {}
    
    void set_center(Vertex3D center){
        center_ = center;
    }
    
    void set_v1(Vector3d v1){
        v1_= v1;
    }
    
    void set_v2(Vector3d v2){
        v2_=v2;
    }
    
    void set_width(double width){
        width_ = width;
    }
    
    void set_height(double height){
        height_ = height;
    }
    
    Vertex3D center(){
        return center_;
    }
    
    Vector3d v1(){
        return v1_;
    }
    
    Vector3d v2(){
        return v2_;
    }
    
    double width(){
        return width_;
    }
    
    double height(){
        return height_;
    }
    
    Colour ambiance1(){
        return ambiance1_;
    }
    
    Colour difusing1(){
        return difusing1_;
    }
    
    Colour reflecting1(){
        return reflecting1_;
    }
    
    double n1(){
        return n1_;
    }
    
    double kref1(){
        return kref1_;
    }
    
    Colour ambiance2(){
        return ambiance2_;
    }
    
    Colour difusing2(){
        return difusing2_;
    }
    
    Colour reflecting2(){
        return reflecting2_;
    }
    
    double n2(){
        return n2_;
    }
    
    double kref2(){
        return kref2_;
    }
    
    void updateIntersection(Intersection& inters, Vector3d start, Vector3d d) override{
        
        
        Matrix<double> system1(std::vector<std::vector<double>>{std::vector<double>{v1_[0], v2_[0], -d[0]},
            std::vector<double>{v1_[1], v2_[1], -d[1]}, std::vector<double>{v1_[2], v2_[2], -d[2]}});
        
        Matrix<double> system2(std::vector<std::vector<double>>{std::vector<double>{start[0] - center().x()}, std::vector<double>{start[1] - center().y()}, std::vector<double>{start[2] - center().z()}});
        double detSys1 = system1.determinant();
        
  
                               
        Matrix<double> x1(std::vector<std::vector<double>>{{system2[0][0], system1[0][1], system1[0][2]}, {system2[1][0], system1[1][1], system1[1][2]}, {system2[2][0], system1[2][1], system1[2][2]}});
        
        Matrix<double> y1(std::vector<std::vector<double>>{{system1[0][0], system2[0][0], system1[0][2]}, {system1[1][0], system2[1][0], system1[1][2]}, {system1[2][0], system2[2][0], system1[2][2]}});
        
        Matrix<double> z1(std::vector<std::vector<double>>{{system1[0][0], system1[0][1], system2[0][0]}, {system1[1][0], system1[1][1], system2[1][0]}, {system1[2][0], system1[2][1], system2[2][0]}});
        
        
        double lambda = x1.determinant()/detSys1;
        double mi = y1.determinant()/detSys1;
        double eta = z1.determinant()/detSys1;
        
        if((-width()/2.0 <= lambda <= width()/2.0) && (-height()/2.0 <= mi <= height()/2.0) && eta >= 0 && (inters.lambda_ > lambda)){
            inters.object_ = *this;
            inters.lambda_ = lambda;
            inters.front_ = true;
            Vertex3D point(center().x() + lambda * v1()[0] + mi * v2()[0], center().y() + lambda * v1()[1] + mi * v2()[1], center().z() + lambda * v1()[2] + mi * v2()[2]);
            inters.point_ = point;

            
        }
        
        
        return;
    }
    
 
    Vector3d getNormalInPoint(Vertex3D point) override{
        
        return (v1()*v2()).normalize();
    }
    
    
};
