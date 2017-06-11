#ifndef prim_2d_h_
#define prim_2d_h_

struct Vector2D {
   double x, y;
   explicit Vector2D(const double & x = .0, const double & y = .0);
   Vector2D & operator=(const Vector2D &);
   Vector2D  operator -() const;
   Vector2D &  operator +=(const Vector2D &);
   double module() const;
};

struct Segm2D {
   Vector2D A,B;
   Segm2D(const Vector2D & _A = Vector2D(.0,.0), const Vector2D & _B = Vector2D(.0,.0));
   Segm2D & operator=(const Segm2D &);
};

Vector2D operator+ (const Vector2D &, const Vector2D &);
Vector2D operator- (const Vector2D &, const Vector2D &);
Vector2D operator* (const double &, const Vector2D &);
Vector2D operator* (const Vector2D & , const double &);
Vector2D operator/ (const Vector2D & , const double &);
bool operator== (const Vector2D &, const Vector2D &);
double dotProd(const Vector2D &, const Vector2D &);
double crossProd(const Vector2D &, const Vector2D &);
std::ostream & operator<<(std::ostream & , const Vector2D & );
std::istream & operator>>(std::istream &, Vector2D &);

std::ostream & operator<<(std::ostream & , const Segm2D & );
std::istream & operator>>(std::istream & , Segm2D & ) ;

#endif

