#include <cmath>
#include <iostream>
#include <fstream>
#include <iterator>
#include <string>
#include <iomanip>
#include "de2d.h"

double Omega(const Vector2D & rM, const Vector2D & rN, const Vector2D nN)
{
    Vector2D rNM = rM - rN; 
    Vector2D gradPhi = (0.5/M_PI/pow(rNM.module(),2.0))*rNM;
    return dotProd(gradPhi,nN);
}

double dl(const Segm2D & S)
{
   return (S.B-S.A).module();
}

Vector2D normal(const Segm2D & S)
{
    Vector2D r12 = S.B - S.A;
    Vector2D n = Vector2D( -r12.y, r12.x );
    return (1.0/n.module())*n;
}

double phi0(const Vector2D & rM, const double & q1, const Vector2D & r1)
{
    return .5*q1*M_1_PI*log((rM-r1).module());
}

Vector2D v0(const Vector2D & rM, const double & q1, const Vector2D & r1)
{
    Vector2D r = rM - r1;
    return .5*q1*M_1_PI*r/pow(r.module(),2.0);
}

Vector2D vVortex(const Vector2D & rM, const Vector2D & A)
{
    if(rM==A) return Vector2D(0.0,0.0);
    double r = (rM-A).module();
    double dPsi2x = 0.5/M_PI*(rM.x-A.x)/pow(r,2.0);
    double dPsi2y = 0.5/M_PI*(rM.y-A.y)/pow(r,2.0);
    return Vector2D(dPsi2y,-dPsi2x);
}

double Theta(const Vector2D & rM, const Vector2D & A, const double & VEps){
   if(rM==A) return 0.0;
   Vector2D rAM =  rM - A ;
   double r = rAM.module();
   double koeff = 1.0;
   if(r<VEps) 
   koeff = (63*pow(r/VEps,5.0)-90*pow(r/VEps,7.0)+35*pow(r/VEps,9.0))/8.0;
   return koeff;
}

Vector2D vSegmTheta(const Vector2D & rM, const Segm2D & s, const double & VEps)
{
  return vVortexTheta(rM,s.A,VEps)-vVortexTheta(rM,s.B,VEps);
}

Vector2D vVortexTheta(const Vector2D & rM, const Vector2D & A, const double & VEps){
   if(rM==A) return Vector2D(0.0,0.0);
   double r = (rM-A).module();
   double dPsi2x = 0.5/M_PI*(rM.x-A.x)/pow(r,2.0);
   double dPsi2y = 0.5/M_PI*(rM.y-A.y)/pow(r,2.0);
   double koeff = Theta(rM,A,VEps);
   return koeff*Vector2D(dPsi2y,-dPsi2x);
}

void save_border(const std::vector<Segm2D> & f, const std::string & s)
{
    std::ofstream file(s.c_str());
    size_t n = f.size();
    std::copy
    (
        f.begin(),
        f.end(),
        std::ostream_iterator<Segm2D>(file,"\n")
    );
}

void load_border(std::vector<Segm2D> & f, const std::string & s)
{
    std::ifstream file(s.c_str());
    std::istream_iterator<Segm2D> is(file), eof;
    std::copy(is, eof, std::back_inserter(f) );
}
