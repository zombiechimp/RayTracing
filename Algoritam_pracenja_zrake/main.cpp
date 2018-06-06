//
//  main.cpp
//  Algoritam_pracenja_zrake
//
//  Created by Tea Jakić on 01/06/2018.
//  Copyright © 2018 Tea Jakić. All rights reserved.
//

#include "classes.cpp"
#include <vector>
#include <GLUT/glut.h>
#include <float.h>
#include <string.h>
#define BasicAngle 18.4349488
#define BasicR 3.16227766
#define increment 1
#define radToC 57.295
#define PI 3.14159265

using namespace std;

double angle = BasicAngle;
double r = BasicR;
GLuint width = 800;
GLuint height = 600;

Vertex3D eye;               /*ociste*/
Vector3d view;
Vector3d viewUp;
Colour gl_ambient_light;
double  dist;           /*udaljenost ravnine prikaza od ocista*/
double angle_h;             /*horizontalni kut gledanja*/
double angle_v;             /*vertikalni kut gledanja*/
double l_;
double r_;
double b_;
double t_;
Vertex3D G;
Vector3d xAxis;
Vector3d yAxis;

std::vector<Source> sources;
std::vector<Orb> orbs;
std::vector<Krpica> krpice;


Colour pixel_cols[600][800];

void reshape(int width, int height);
void display();

char test[600][800];



class Ray{
public:
    Vertex3D start_;
    Vector3d d_;
    
    Ray(double x, double y, double z){
//        start_ = Vertex3D(x, y, z);
        start_ = eye;
        Vector3d end({x, y, z});

        Vector3d eyez({eye.x(), eye.y(), eye.z()});

        d_= eyez - end;
        if(d_.magnitude()) d_ = d_.normalize();
    }
    
};


void computeKS(){
    double v_norm = sqrt(view[0]*view[0] + view[1]*view[1] + view[2]*view[2]);
    G.x_ = eye.x() + dist*(view[0]/view.magnitude());
    G.y_ = eye.y() + dist*(view[1]/view.magnitude());
    G.z_ = eye.z() + dist*(view[2]/view.magnitude());
    
    l_ = r_ = dist*tan(angle_h/(2*radToC));
    t_ = b_ = dist*tan(angle_v/(2*radToC));
    
    xAxis = view * viewUp;
    xAxis.normalize();
    yAxis = xAxis * viewUp;
    yAxis.normalize();
    
    
    
}


Colour setColour(Object& object, Ray ray, Intersection inters){

//    Vertex3D point(inters.lambda_ * ray.d_[0] + ray.start_.x(), inters.lambda_ * ray.d_[1] + ray.start_.y(), inters.lambda_ * ray.d_[2] + ray.start_.z());
//    Colour black(0,0,0);
//
//    for(int i = 0; i < sources.size(); i++){
//
//    }
    return Colour(1,1,1);
}



Ray ray_from_eye_to_point(int x, int y){
    
    return Ray(x, y, -dist);
    
}




Colour next(Ray ray){
    Intersection inters(orbs[0]);
    Vector3d start({ray.start_.x(), ray.start_.y(), ray.start_.z()});
    for(int i = 0; i < orbs.size(); i++){
        orbs[i].updateIntersection(inters, start, ray.d_);
        
    }
    for(int i = 0; i < krpice.size(); i++){
        krpice[i].updateIntersection(inters, start, ray.d_);
    }
    if(inters.lambda_ != INFINITY)
        return setColour(inters.object_, ray, inters);
    return Colour(0,0,0);
}



void raytrace(){
    Colour colour(1,1,1);
    
    for(int y = 0; y < height; y ++){
        for (int x = 0; x < width; x++){

            Ray ray = ray_from_eye_to_point(x, y);
            pixel_cols[y][x] = next(ray);
            if(colour.r_ == 1){
//                cout<<x<<" "<<y<<endl;
            }


            
        }
    }
}

int main(int argc, char * argv[]) {
    vector<Vertex3D> vertices;
    vector<Face3DInd> faces;
    
    string line;
    
    ifstream file;
    system("pwd");
    file.open( (argv[1]));
    
    if(!file.is_open()){
        printf("Nema filea! Ne moze ga se otvoriti! \n");
        return 0;
    }
    
    for (string line; getline(file, line);)
    {
        stringstream sstream(line);
        string entity_type = "\0";
        sstream >> entity_type;
        if(entity_type == "e"){
                double x = 0, y = 0, z = 0;
                if (!(sstream >> x >> y >> z))
                {
                    throw domain_error("vertex format invalid");
                }
                eye = Vertex3D(x, y, z);
            }
        else if(entity_type == "v") {
                
                double x = 0, y = 0, z = 0;
                if (!(sstream >> x >> y >> z))
                {
                    throw domain_error("vertex format invalid");
                }
                
                view = Vector3d({x, y, z});
                

            }
        else if(entity_type == "vu") {
                double x = 0, y = 0, z = 0;
                
                if (!(sstream >> x>> y >> z))
                {
                    throw domain_error("face format invalid");
                }
                
                viewUp = Vector3d({x, y, z});
                

            }
                
            else if(entity_type == "h"){
                if (!(sstream >> dist))
                {
                    throw domain_error("face format invalid");
                }

            }
            else if(entity_type ==  "xa"){
                if (!(sstream >> angle_h))
                {
                    throw domain_error("face format invalid");
                }

            }
            else if(entity_type == "ya"){
                if (!(sstream >> angle_v))
                {
                    throw domain_error("face format invalid");
                }
 
            }
            else if(entity_type == "ga"){
                double r = 0, g= 0, b = 0;
                if (!(sstream >> r >> g >> b))
                {
                    throw domain_error("vertex format invalid");
                }
                gl_ambient_light = Colour(r, g, b);
  
            }
            else if(entity_type == "i"){
                double x = 0, y = 0, z = 0, r = 0, g = 0, b = 0;
                if (!(sstream >> x >> y >> z >> r >> g >> b))
                {
                    throw domain_error("vertex format invalid");
                }
                sources.emplace_back(Source(x, y, z, r, g, b));

            }
            else if(entity_type == "o"){
                string entity_type = "\0";
                sstream >> entity_type;
                if(entity_type == "s"){
                
                    double x = 0, y = 0, z = 0, r = 0, ar = 0, ag = 0, ab = 0, dr = 0, dg = 0, db = 0, rr = 0, rg = 0, rb = 0, n = 0, kref = 0;
                    if (!(sstream >> x >> y >> z >> r >> ar >> ag >> ab >> dr >> dg >> db >> rr >> rg >> rb >> n >> kref))
                    {
                        throw domain_error("vertex format invalid");
                    }
                
                    orbs.emplace_back(Orb(x, y, z, r, ar, ag, ab, dr, dg, db, rr, rg, rb, n, kref));
                }
                else if(entity_type == "p"){
                    double x = 0, y = 0, z = 0, v1x = 0, v1y = 0, v1z = 0, v2x = 0, v2y = 0, v2z = 0, w = 0, h = 0,ar = 0, ag = 0, ab = 0, dr = 0, dg = 0, db = 0, rr = 0, rg = 0, rb = 0, n = 0, kref = 0, ar2 = 0, ag2 = 0, ab2 = 0, dr2 = 0, dg2 = 0, db2 = 0, rr2 = 0, rg2 = 0, rb2 = 0, n2 = 0, kref2 = 0;
                    if (!(sstream >> x >> y >> z >> v1x >> v1y >> v1z >> v2x >> v2y >> v2z >> w >> h >> ar >> ag >> ab >> dr >> dg >> db >> rr >> rg >> rb >> n >> kref >>  ar2 >> ag2 >> ab2 >> dr2 >> dg2 >> db2 >> rr2 >> rg2 >> rb2 >> n2 >> kref2)){
                        throw domain_error("vertex format invalid");
                    }
                    krpice.emplace_back(Krpica(x, y, z, v1x, v1y, v1z, v2x, v2y, v2z, w, h, ar, ag, ab, dr, dg, db, rr, rg, rb, n, kref,ar2, ag2, ab2, dr2, dg2, db2, rr2, rg2, rb2, n2, kref2));

            }
        }
    }
    computeKS();

    raytrace();

    
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
    glutInitWindowSize(800, 600);
    glutInitWindowPosition(0, 0);
    glutCreateWindow("Pls work");
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutReportErrors();

    glutMainLoop();
    
    return 0;
}

void display() {
    
    
 
//    gluPerspective(angle_v/radToC, (float)width/height, 0.5, 8.0);
//    gluLookAt (eye.x(), eye.y(), eye.z(), G.x(), G.y(), G.z(), viewUp[0], viewUp[1], viewUp[2]);


  
    glPointSize(1.0);
    for(int y = 0; y < height; y++){
        for(int x = 0; x < width; x++){
            glBegin(GL_POINTS);
            glColor3f(pixel_cols[x][y].r_, pixel_cols[x][y].g_, pixel_cols[x][y].b_);
            glVertex2i(x, y);
            glEnd();
            
        }

    }

    

    glutSwapBuffers();

    glFlush();
}

void reshape(int w, int h) {
    
    width = w;
    height = h;
    
    glMatrixMode (GL_PROJECTION);        // aktivirana matrica projekcije
    glLoadIdentity ();
//    gluPerspective(angle_v/radToC, (float)width/height, 1, dist); // kut pogleda, x/y, prednja i straznja ravnina odsjecanja
    glViewport(0, 0, width, height);
    gluOrtho2D(0, width-1, height-1, 0);
    glMatrixMode (GL_MODELVIEW);         // aktivirana matrica modela
    glLoadIdentity ();
//    gluLookAt (eye.x(), eye.y(), eye.z(), G.x(), G.y(), G.z(), viewUp[0], viewUp[1], viewUp[2]);
}
    

