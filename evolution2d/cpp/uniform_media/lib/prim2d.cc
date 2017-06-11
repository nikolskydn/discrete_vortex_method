#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>
#include "prim2d.h"

Vector2D::Vector2D(const double & _x, const double & _y )
    :x(_x), y(_y) {}

Vector2D Vector2D::operator -() const 
{ 
    return Vector2D(-x,-y); 
}

Vector2D& Vector2D::operator +=(const Vector2D & rV) 
{
    x += rV.x; 
    y += rV.y;  
    return *this;
}

Vector2D & Vector2D::operator=(const Vector2D & A)
{
    x=A.x;
    y=A.y; 
    return *this; 
}

double Vector2D::module() const 
{ 
    return sqrt(x*x+y*y); 
}

Vector2D operator+ (const Vector2D & A,  const Vector2D & B)
{
   return Vector2D(A.x+B.x, A.y+B.y);
}

Vector2D operator- (const Vector2D & A,  const Vector2D & B)
{
   return Vector2D(A.x-B.x, A.y-B.y);
}

bool operator== (const Vector2D & A, const Vector2D & B)
{
   return (fabs(A.x-B.x) < std::numeric_limits<double>::epsilon()) && (fabs(A.y-B.y)<std::numeric_limits<double>::epsilon());
}

Vector2D operator* (const double & d, const Vector2D & A)
{
   return Vector2D(d*A.x,d*A.y);
}

Vector2D operator* (const Vector2D & A, const double & d)
{
   return Vector2D(d*A.x,d*A.y);
}

Vector2D operator/ (const Vector2D & A, const double & d)
{
   return Vector2D(A.x/d,A.y/d);
}

std::ostream & operator<<(std::ostream & s, const Vector2D & P)
{
    s  << std::fixed << std::showpoint<< std::setprecision(10)
       << P.x << ' ' << P.y  << std::endl;  return s;
}



std::istream & operator>>(std::istream & s, Vector2D & P)
{
    s >>  P.x  >> P.y  ;  
    return s;
}
   

double dotProd(const Vector2D & A, const Vector2D & B) 
{
  return A.x*B.x+A.y*B.y;
}

double crossProd(const Vector2D & A, const Vector2D & B) 
{
   return A.x*B.y-A.y*B.x;
}

Segm2D::Segm2D(const Vector2D & _A, const Vector2D & _B) 
    : A(_A), B(_B) {}

Segm2D & Segm2D::operator=(const Segm2D & S)
{
    A=S.A; 
    B=S.B; 
    return *this;
}

std::ostream & operator<<(std::ostream & s, const Segm2D & F){
   s << F.A  << F.B   << std::endl  ; 
   return s;
}
   
   
std::istream & operator>>(std::istream & s, Segm2D & F) {
  s >> F.A >> F.B  ;
}
