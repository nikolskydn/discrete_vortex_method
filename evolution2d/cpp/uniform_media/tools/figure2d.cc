#include <fstream>
#include <algorithm>
#include <iomanip>  
#include <iterator> 
#include <cmath>
#include <string>

#include "figure2d.h"

#include <iostream>

El2D::El2D(const size_t _n, const double & _a, const double & _b, std::vector<Segm2D>  & f)
{
    n=_n;  a=_a; b=_b; 
    for(size_t m=0; m<n-1; ++m) {
        double xA, yA, xB, yB, zB, theta1, theta2;
        theta2 = m*2.0*M_PI/static_cast<double>(n-1);
        theta1 = (m+1)*2.0*M_PI/static_cast<double>(n-1);
        xA=a*cos(theta1);  
        yA=b*sin(theta1);
        xB=a*cos(theta2);  
        yB=b*sin(theta2);
        Vector2D A(xA,yA), B(xB,yB) ;
        if(  ! (A==B)  )
        {
           Segm2D F(A,B);
           f.push_back(F);
        }
   }//for_m
}

/*
Line2D::Line2D(const size_t _n, const double & _y, const double & _l, std::vector<Segm2D>  & f)
{
    n=_n;  y=_y; l=_l;
    for(size_t m=0; m<n-1; ++m) {
        Vector2D A , B ;
        A.y = B.y = y;
        A.x = -l + 2.0*l*m/static_cast<double>(n-1.0);
        B.x = -l + 2.0*l*(m+1.0)/static_cast<double>(n-1.0); 
        Segm2D F(A,B);
        f.push_back(F);
   }//for_m
}


LineTheta2D::LineTheta2D(const size_t _n, const double _y, const double _th, std::vector<Segm2D>  & f)
{
    n=_n;  y=_y; th=_th*M_PI/180; double h=1.0; double theta1, theta2;
    for(size_t m=0; m<n-1; ++m) {
        Vector2D A , B ;
        A.y = B.y = y;
        theta1=-th+2.0*th*m/static_cast<double>(n-1.0);
        theta2=-th+2.0*th*(m+1.0)/static_cast<double>(n-1.0);
        A.x = h * tan(theta1) ; 
        B.x = h * tan(theta2); 
        Segm2D F(A,B);
        f.push_back(F);
   }//for_m
}
*/
