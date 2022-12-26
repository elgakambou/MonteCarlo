// vim: set sw=4 ts=4 sts=4:

#include <iostream>
#include <ctime>
#include "MonteCarlo.hpp"
#include "BSCall.hpp"
#include "Exo1.hpp"

using namespace std;
int main()
{

     PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
     pnl_rng_sseed(rng, std::time(NULL));

    // tp2
    // exo 1 
    Exo1 exo1(2.0, 0.25, 0.02, 100.0, 0.95, rng);
    double price, std_err; ;
    exo1.mc_classique(price, std_err);
    exo1.mc_anthitetique_sur_M(price, std_err);
    exo1.mc_anthitetique_sur_M_sur_2(price, std_err);
    exo1.mc_var_controle1( price,  std_err);
    exo1.mc_var_controle2( price,  std_err);

    // exo 2 : IMPORTANCE SAMPLING
    cout << " \n\n exo 2 : IMPORTANCE SAMPLING \n\n" <<endl;
     double prix, stddev;
    BSCall product(2., 0.2, 0.03, 100., 120.);
    MonteCarlo pricer(product, 50000);
    pricer.mc(prix, stddev, rng);

    //test de MonteCarlo exo 2
    double prixMonteCarlo, stddevMonteCarlo;
    double approxTheta = 0;
    int n = 50000;
    PnlVect* lambdas = pnl_vect_create(n);
    double gamma = 0.5;
    pricer.is(lambdas, gamma, n, rng, approxTheta);
    pricer.mcis(prixMonteCarlo, stddevMonteCarlo, approxTheta, rng);
    //std::cout << "thetaAttendu : " << log(120.0/100) * 5 << "\n";
    std::cout << "thetaopt : " << approxTheta << "\n";
    std::cout << "prix MonteCarlo : " << prixMonteCarlo << " (IC = " << stddevMonteCarlo * 1.96 << ")\n";
    // FILE* f = fopen("output.txt", "w");
    // pnl_vect_fprint(f, lambdas);
    // fclose(f);
   

    pnl_rng_free(&rng);
    exit(0);
}
