#include <vector>
#include <string>

#include "figure2d.h"
#include "../lib/de2d.h"


int main( int argc, char* argv[] )
{
    size_t n = 401;
    double a = 1.0;
    double b = 1.0;
    std::string fileName = "b0";

    std::vector<Segm2D> border;
    El2D e(n,a,b,border);
    save_border(border, fileName);
}

