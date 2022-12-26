// vim: set sw=4 ts=4 sts=4:

#include <iostream>
#include <ctime>
#include "MonteCarlo.hpp"
#include "BSBarrier.hpp"
#include "ImportanceSampling.hpp"

int main()
{
    double prix, stddev;
    BSBarrier product(2., 0.2, 0.05, 100., 110., 80, 24);
    MonteCarlo pricer(product, 50000);
    PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng, std::time(NULL));

    double price, stderr;
    std::cout << "Prix à trouver :" << 11.2 << "\n";
    ImportanceSampling impSampl (rng, product, 50000);
    pricer.mc(prix, stddev, rng);
    impSampl.importance_sampling_mc(price, stderr);
    std::cout << "\n On converge plus vite par importance sampling, on a besoin de moins de <<samples>> pour atteindre le resultat\n";

    pnl_rng_free(&rng);
    exit(0);
}
