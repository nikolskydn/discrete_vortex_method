#include <boost/numeric/ublas/vector.hpp> 
#include <boost/numeric/ublas/matrix.hpp> 
#include <boost/math/tools/solve.hpp> 
#include <iostream>
#include <iomanip>
#include <regex>
#include <omp.h>
#include "../lib/prim2d.h"
#include "../lib/de2d.h"

#define ub  boost::numeric::ublas

#define DEBUG 2

int main(int argc, char* argv[]){

    if( !argv[1] ) 
    {
        std::cout << "set the input file name\n";
        exit(1);
    }
    std::string fileName0(argv[1]);

    std::vector<Segm2D> sigma_t;
    sigma_t.reserve(1000);

    load_border(sigma_t, fileName0);
    size_t n = sigma_t.size();
    #if DEBUG > 1
       std::cout << "\033[31;1mn = " << n << "\033[0m\n";
    #endif

    bool isMove = true;
    int cycle = 0;
    double t = 0.0;
    double dt = 0.01;
    const double T = 5.;
    const double lambda = 0.5;
    const Vector2D rM1 = Vector2D(2.5,0.0);
    const double q1 = -M_PI;
    const double T0 = 1.0;
    const double VEps = 0.1;
    const size_t max_cycle = 20;
    ub::identity_matrix<double> I (n);

    ub::matrix<double> A(n,n);
    ub::vector<double> f(n);
    ub::vector<double> g(n), gp(n);
    double err, eps = 1e-5;;

    for(;cycle<max_cycle&&isMove;++cycle)
    { 

        // ***********   stage 1 *****************// 


        double ts1 = omp_get_wtime();
        #pragma omp parallel for
        for(size_t i=0; i<n; ++i)
        {
            Vector2D rM = (sigma_t[i].A + sigma_t[i].B)/2.0;
            for(size_t j=0; j<n; ++j)
            {
                if(i!=j)
                {
                    Vector2D rN = 0.5*(sigma_t[j].A + sigma_t[j].B );
                    Vector2D nN = normal(sigma_t[j]);
                    A(i,j)=-2.0*lambda*Omega(rM,rN,nN)*dl(sigma_t[j]);
                }
                else A(i,j)=1.0;
            }
            f(i)=2.0*lambda*phi0(rM,q1,rM1); 
        }
        double te1 = omp_get_wtime();
        std::cout << "\t\033[33mstage 1: "
                  << std::fixed << std::showpoint<< std::setprecision(10)
                  << te1-ts1 << " sec.\033[0m";

        // ***********   stage 2 *****************// 
        double ts2 = omp_get_wtime();
        //g = boost::math::tools::solve(A,f);
        err = 1e+10; 

        while ( err > eps )
        {
            #pragma omp parallel for
            for( size_t i=0; i<n; ++i )
                gp(i) = g(i);
        
            #pragma omp parallel for 
            for( size_t i=0; i<n; ++i )
            {
                g(i) = f(i);
                for( size_t j=0; j<n; ++j )
                    g(i) += ( I(i,j)-A(i,j) ) * gp(j) ;
            }
            err=fabs(g(0)-gp(0));

            #pragma omp parallel for reduction(max:err)
            for( size_t i=1; i<n; ++i )
            {
                double tmp = fabs( g(i)-gp(i) );
                if( tmp > err ) 
                    err = tmp;
            }
        } // while

        double te2 = omp_get_wtime();
        std::cout << "\t\033[34mstage 2: " 
                  << std::fixed << std::showpoint<< std::setprecision(5)
                  << te2-ts2 << " sec.\033[0m";


        // ***********   stage 3 *****************// 

        std::vector<Vector2D> vA(n), vB(n);

        double ts3 = omp_get_wtime();
        #pragma omp parallel for
        for(size_t i=0; i<n; ++i){
            Vector2D rMA = sigma_t[i].A;
            vA[i] = v0(rMA,q1,rM1);
            for(size_t j=0; j<n; ++j)
                vA[i] += g(j)*vSegmTheta(rMA,sigma_t[j],VEps);
        }// for_i

        #pragma omp parallel for
        for(size_t i=0; i<n; ++i){
            sigma_t[i].A += vA[i]*dt;
            sigma_t[(i+1)%n].B = sigma_t[i].A;
 
        }// for_i

        double te3 = omp_get_wtime();
        std::cout << "\t\033[35mstage 3: " 
                  << std::fixed << std::showpoint<< std::setprecision(10)
                  << te3-ts3 << " sec. \033[0m\n";

        t+=dt/T0;

        #if DEBUG > 1 
        for( size_t k=1; k<5; ++k)
            if(cycle==(k*(max_cycle-1)/5)) 
            {
                std::stringstream ss;
                ss << k;
                std::string suffix;
                ss >> suffix; 

                std::string outFileName = std::regex_replace(
                    fileName0, 
                    std::regex("0"), 
                    suffix
                );
                std::cout << "\033[31;1mcycle=" << std::setw(4) << cycle 
                          << "\tsave " << outFileName 
                          << "\t\tt = " << t 
                          << "\033[0m\n"; 
                save_border(sigma_t, outFileName); 
            }
        #endif 

        // set your conditions for the variable isMove here:
        // if ( u cond. )  isMove = false;
   
    } // cycle

    std::cout <<  "\033[32;1msimulation time T = " << T << "\033[0m\n";
    std::string outFileName = std::regex_replace(
        fileName0, 
        std::regex("0"), 
        "_end"
    );
    save_border(sigma_t, outFileName); 

}//main
