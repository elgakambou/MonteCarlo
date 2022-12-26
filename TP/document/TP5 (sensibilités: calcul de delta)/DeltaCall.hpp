#pragma once
#include "pnl/pnl_vector.h"
#include "pnl/pnl_random.h"

class DeltaCall   {

    public : 
        double m_S0;
        double m_sigma;
        double m_r;
        double m_T;
        double m_K;

    DeltaCall (double spot, double sigma, double r, double T, double K);

    double payoffCall(double S0,double sigma, double r, double T, double K, double g);
    double derive_payoffCall(double S0,double sigma, double r, double T, double K, double g);
    void delta_df(double&price, double &std_err, double epsilon, PnlRng* rng);
    void delta_likelihood (double&price, double &std_err, PnlRng* rng);
    void delta_pathwise (double&price, double &std_err, PnlRng* rng);
    private :
        int M = 50000;

};