// vim: set sw=4 ts=4 sts=4:
#include <iostream>
#include <cmath>
#include "MonteCarlo.hpp"

MonteCarlo::MonteCarlo(BSCall &product, size_t samples)
    : m_product(product), m_samples(samples)
{}


void MonteCarlo::mc(double &prix, double &stddev, PnlRng *rng)
{
    double sum = 0.;
    double var = 0.;
    for (size_t i = 0; i < m_samples; i++) {
        double g = pnl_rng_normal(rng);
        double flow = m_product.payoff(g);
        sum += flow;
        var += flow * flow;
    }
    prix = sum / m_samples;
    var = var / m_samples - prix * prix;
    stddev = std::sqrt(var / m_samples);
    std::cout << "prix : " << prix << " (half-IC = " << stddev * 1.96 << ")\n";

}

void MonteCarlo::is(PnlVect *lambda, double gamma, int n, PnlRng *rng, double &lambdaApprox) {

    double beta  = 0.75;
    double theta0 = 0;
    double theta = theta0;
    double alpha = 0;
    

    for (int i = 1; i <= n; i++) {  
        double g = pnl_rng_normal (rng);
        double thetaUnDemi = theta - m_product.gamma(gamma, i) * m_product.computeU(theta, g );

        if ( pow( thetaUnDemi, 2) <= log(alpha + 1) ) {
            theta = thetaUnDemi;         
        }
        else {
            theta = theta0;
            alpha +=1;
        }

        if (i == n) { // c'est la derniere valeur
            lambdaApprox = theta;
        }

        pnl_vect_set(lambda, i-1, theta);          
    }
}

void MonteCarlo::mcis(double &prix, double &stddev, double lambda, PnlRng *rng) {
    int M = this->m_samples;

    double somme = 0;
    double somme2 = 0;

    for (int i = 0; i < M; i++){
        double g = pnl_rng_normal (rng);
        double val = m_product.payoff(g, lambda) * m_product.weight_minus( g, lambda);
        somme += val;
        somme2 += val*val;
    }
    prix = somme/M;
    stddev = somme2/M - prix*prix;
    stddev = sqrt(stddev/M);

    //commentaire : le thetaOpt ne change pas la moyenne mais reduire la variance d'ou l'interet de la chose 
}
