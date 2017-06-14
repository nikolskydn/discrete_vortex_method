
#ifndef _FIGURE3D_H_
#define _FIGURE3D_H_

#include <vector>
#include <string>

#include "../lib/prim2d.h"


class El2D{
  size_t n;
  double a,b;
public:
   El2D(const size_t _n, const double & _a, const double & _b,  std::vector<Segm2D>  & f);
};


/*
class Line2D{
  size_t n;
  double y,l;
public:
  Line2D(const size_t _n, const double _y, const double _l, std::vector<Segm2D>  & f);
};

class LineTheta2D{
  size_t n;
  double y,th;
public:
  LineTheta2D(const size_t _n, const double _y, const double _th, std::vector<Segm2D>  & f);
};
*/
#endif
