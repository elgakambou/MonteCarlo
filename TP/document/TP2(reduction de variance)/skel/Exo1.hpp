// vim: set sw=4 ts=4 sts=4:
#pragma once
#include "pnl/pnl_random.h"
using namespace std;



class Exo1 {
public:
    Exo1(double maturity, double volatility, double interest_rate, double spot, double strike, PnlRng *rng);
    double m_maturity;
    double m_volatility;
    double m_interest_rate;
    double m_spot;
    double m_strike;
    PnlRng *m_rng;

    int M = 50000;
    double computeSclassique( int nbValues);
    double computeSantitetique1(int nbValues);// retourne le payoffanti1 f(Xi)
    double computeSantitetique2(int nbValues); // retourne le payoffanti2 f(-Xi)
    // tout depend de cov(f(Xi), f(-Xi)), plus f est lineaire mieux c'est 
    double  compute_var_controle1(int nbValues);
    double  compute_var_controle2(int nbValues, double& maxPath);

    void mc_classique(double &price, double &std_err);
    void mc_anthitetique_sur_M(double &price, double &std_err);
    void mc_anthitetique_sur_M_sur_2(double &price, double &std_err);
    void mc_var_controle1(double &price, double &std_err);
    void mc_var_controle2(double &price, double &std_err);



    // /// Calculer une valeur de S_T
    // double asset(double G, double theta = 0.);
    // /// Calculer psi(G)
    // double payoff(double G, double theta = 0.);
    // /// Calculer exp(- theta G + theta^2/2)
    // double weight_plus(double G, double theta);
    // /// Calculer exp(- theta G - theta^2/2)
    // double weight_minus(double G, double theta);
    // // Calculer la dérivée de weight_plus
    // double d_weight_plus(double G, double theta);
    // // moi
    // double computeU(double param1,  double param2 );
    // double gamma(double theGamma, int n, double beta = 0.75 );

};
