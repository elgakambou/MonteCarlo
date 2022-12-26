// vim: set sw=4 ts=4 sts=4:

#include <iostream>
#include <ctime>
#include "DeltaCall.hpp"
#include "pnl/pnl_finance.h"

int main()
{
    double prix, stddev;

    PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng, std::time(NULL));
    double spot = 100.0;
    double sigma = 0.2;
    double r = 0.05;
    double T = 2.0;
    double K = 110.;
    double epsilon = 0.0668740304976422; // PAS trop petit sinon on ne capte plus les variations  mais les bruits MonteCarlo  , eps > n ^ (-0.25)

    DeltaCall deltacall(spot, sigma, r, T, K);
    double pnl_price =  0.0;
    double pnl_delta = 0.0;

    std::cout << pnl_cf_call_bs(spot, K, T, r, 0.0, sigma, &pnl_price, &pnl_delta)<< "\n";
    std::cout << "Prix à trouver : " <<  pnl_price << "  delta à trouver : " << pnl_delta << "\n";

    deltacall.delta_df(prix, stddev, epsilon, rng);
    deltacall.delta_likelihood(prix, stddev, rng);
    deltacall.delta_pathwise(prix, stddev, rng);
    pnl_rng_free(&rng);
    exit(0);
}