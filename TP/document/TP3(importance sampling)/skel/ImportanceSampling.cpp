#include "ImportanceSampling.hpp"
#include <iostream>


ImportanceSampling::ImportanceSampling  ( PnlRng * rng, BSBarrier product, int samples) :
    myRng(rng), myProduct(product), m_samples(samples) {}
        

void ImportanceSampling :: computeIncrementBrownien(PnlMat * GMat) {
    int size = myProduct.m_dates;
    PnlVect * G = pnl_vect_create(size);
    for (int i = 0; i < n; i ++) {
    pnl_vect_rng_normal(G, size, myRng);    
    pnl_mat_set_row(GMat, G, i);
    }
    pnl_vect_free(&G);

}

double ImportanceSampling:: vp(double lambda, PnlMat * Gmat, PnlVect * fcarre) {
    
    double output = 0.0;
    PnlVect * G =  pnl_vect_create(myProduct.m_dates); 
    for (int i = 0; i < n ; i ++) {
        pnl_mat_get_row(G, Gmat, i);
        double WT = myProduct.compute_BT(G);
        double d_exp_plus = myProduct.d_weight_plus(G , lambda);
        double f2 = pnl_vect_get(fcarre, i);
        output += d_exp_plus * f2;
    }
    pnl_vect_free(&G);
    return output/n;
   }

double ImportanceSampling:: vpp(double lambda, PnlMat * Gmat, PnlVect * fcarre) {
    
    double output = 0.0;
    PnlVect * G =  pnl_vect_create(myProduct.m_dates); 
    for (int i = 0; i < n ; i ++) {
        pnl_mat_get_row(G, Gmat, i);
        double WT = myProduct.compute_BT(G);
        double exp_plus = myProduct.weight_plus(G , lambda);
        double f2 = pnl_vect_get(fcarre, i);
        output += (  myProduct.m_maturity +   pow( lambda*myProduct.m_maturity - WT, 2) )  * exp_plus * f2;
    }
    pnl_vect_free(&G);
    return output/n;
   }



void ImportanceSampling :: computeLambaOpt(double &lambdaOpt) {
    int ncol = myProduct.m_dates;
    PnlMat * Gmat = pnl_mat_create(n, myProduct.m_dates);

    // calcul des browniens et stochage 
    computeIncrementBrownien(Gmat);

    PnlVect * fcarre = pnl_vect_create(n);
    PnlVect * path = pnl_vect_new();
    PnlVect * G = pnl_vect_create(ncol);

    // calcul de fcarre
    for (int i = 0; i < n ; i++) {
        pnl_mat_get_row(G, Gmat, i);
        myProduct.asset(path, G);
        double payoff = myProduct.payoff(path);
        pnl_vect_set(fcarre, i, payoff*payoff);
    }

    // calcul de lambda
    double lambda = 0.0;
    double lambda1 = 0.0;
    double e = 5.0;
    //int N = 500;
    while (e > 0.0001)
    {
        lambda1  += -vp(lambda, Gmat, fcarre)/vpp(lambda, Gmat, fcarre);
        e = lambda1 - lambda;
        lambda = lambda1;
    }
    



    // for (int i = 0; i < N; i++) {
    //     lambda  += -vp(lambda, Gmat, fcarre)/vpp(lambda, Gmat, fcarre);
    // }
    lambdaOpt = lambda;

    pnl_vect_free(&fcarre);
    pnl_vect_free(&path);
    pnl_vect_free(&G);
    pnl_mat_free(&Gmat);


}

void ImportanceSampling::importance_sampling_mc (double &prix, double &stddev) {
    double lambdaOpt = 0.0;
    int ncol = myProduct.m_dates;
    computeLambaOpt(lambdaOpt);

    PnlVect * path = pnl_vect_new();
    PnlVect * G = pnl_vect_create(ncol);

    double somme = 0.0;
    double somme2 = 0.0;
    double tmp;
    for (int i = 0; i < m_samples; i++) {
        pnl_vect_rng_normal(G, ncol, myRng); 
        myProduct.asset(path, G, lambdaOpt);
        tmp = myProduct.payoff(path) * myProduct.weight_minus(G , lambdaOpt);
        somme += tmp;
        somme2 += tmp*tmp;
    }
    prix = somme / m_samples;
    double var = somme2 / m_samples - prix * prix;
    stddev = std::sqrt(var / m_samples);
    std::cout << "Imp Sampling  price : " << prix << " (half-IC = " << stddev * 1.96 << ")  avec lambdaOpt : " << lambdaOpt << "\n";

    pnl_vect_free(&path);
    pnl_vect_free(&G);
}
