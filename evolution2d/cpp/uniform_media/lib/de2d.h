#ifndef discrete_elements_2d_h_
#define discrete_elements_2d_h_

#include <vector>
#include "prim2d.h"

double Omega(const Vector2D & rM, const Vector2D & rN, const Vector2D nN);
double dl(const Segm2D &);
Vector2D normal(const Segm2D &);
double phi0(const Vector2D & rM, const double & q1, const Vector2D & r1);
Vector2D v0(const Vector2D & rM, const double & q1, const Vector2D & r1);
Vector2D vSegm(const Vector2D & rM, const Segm2D & sigma);
Vector2D vVortex(const Vector2D & rM, const Vector2D & A);
double Theta (const Vector2D & rM, const Vector2D & A, const double & VEps);
Vector2D vSegmTheta(const Vector2D & rM, const Segm2D & sigma, const double & VEps);
Vector2D vVortexTheta(const Vector2D & rM, const Vector2D & A, const double & VEps);

void save_border(const std::vector<Segm2D> & , const std::string &);
void load_border(std::vector<Segm2D> & , const std::string & );

#endif

