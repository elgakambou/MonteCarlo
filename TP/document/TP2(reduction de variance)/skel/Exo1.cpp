
#include <iostream>
#include <cmath>
#include <algorithm>
#include "Exo1.hpp"

using namespace std;

Exo1::Exo1(double maturity, double volatility, double interest_rate, double spot, double strike, PnlRng *rng) {
    m_maturity = maturity;
    m_volatility = volatility;
    m_interest_rate = interest_rate;
    m_spot = spot;
    m_strike = strike;

    m_rng = rng;
    cout << "Spot " << m_spot << " vol " << m_volatility << " r " << m_interest_rate << " strike " << m_strike << endl;
}

// vector de taille nbval + 1
// remplir le vecteur S et donner le payoff de celui ci
double  Exo1::computeSclassique( int nbValues) {
    double maxS = m_spot;
    double step = m_maturity/nbValues;
    double sqrt_step = sqrt(step);

    double lastVal = m_spot;

    double bs_drift = m_interest_rate - m_volatility * m_volatility / 2;
    
    for (int i=1; i <= nbValues; i++) {
        double g = pnl_rng_normal(m_rng);
        lastVal *= exp(bs_drift * step + m_volatility * sqrt_step * g);
        if (lastVal > maxS) {
            maxS = lastVal;
        }
    }
    double diff = m_strike*maxS - lastVal;
    return diff > 0 ? diff : 0.0;
    //return lastVal - m_strike > 0 ? (lastVal - m_strike)*exp(-m_interest_rate) : 0.0;
}


void Exo1::mc_classique(double &price, double &std_err) {
    int n = 24;
    double somme = 0.0;
    double somme2 = 0.0;
    for (int i = 0; i< M; i ++) { 
        double payoff = computeSclassique(n) ;
        somme += payoff;
        somme2 += payoff*payoff;
    }
    price = somme/M;
    double var = (somme2 / M) - price*price;
    std_err = sqrt(var / M);
    cout << "exo1 : MC classique  price : " << price << " et IC : " << std_err << endl;

}



// ANTITHETIQUE
// simulation d'une trajectoir de S
// avec >S1 et S2 des variables antithetique, c'est à dire de meme loi (symetrie de N(0, p))
// et de covariance negative pour reduire la variance 
// E(S1) = E(S2) = E(f(X))
// voir aussi dans mon doc MC
double Exo1::computeSantitetique1(int nbValues) { 
    double maxS1 = m_spot;

    double step = m_maturity/nbValues;
    double sqrt_step = sqrt(step);

    double lastVal1 = m_spot;

    double bs_drift = m_interest_rate - m_volatility * m_volatility / 2;
    
    for (int i=1; i <= nbValues; i++) {
        double g = pnl_rng_normal(m_rng);

        lastVal1 *= exp(bs_drift * step + m_volatility * sqrt_step * g);

        if (lastVal1 > maxS1) {
            maxS1 = lastVal1;
        }
    }
    double diff1 = m_strike*maxS1 - lastVal1;
    diff1 = diff1 > 0 ? diff1 : 0.0;
    return diff1; // 
    //return lastVal1 - m_strike > 0 ? (lastVal1 - m_strike)*exp(-m_interest_rate) : 0.0; // pour un call ou la cov est plus negative car le payoff est lineaire
}

double Exo1::computeSantitetique2(int nbValues) {
    double maxS2 = m_spot;

    double step = m_maturity/nbValues;
    double sqrt_step = sqrt(step);

    double lastVal2 = m_spot;

    double bs_drift = m_interest_rate - m_volatility * m_volatility / 2;
    
    for (int i=1; i <= nbValues; i++) {
        double g = pnl_rng_normal(m_rng);

        lastVal2 *= exp(bs_drift * step + m_volatility * sqrt_step * (-g));

        if (lastVal2 > maxS2) {
            maxS2 = lastVal2;
        }
    }
    double diff2 = m_strike*maxS2 - lastVal2;
    diff2 = diff2 > 0 ? diff2 : 0.0;
   return diff2; // 
    //return lastVal2 - m_strike > 0 ? (lastVal2 - m_strike)*exp(-m_interest_rate) : 0.0; // pour un call ou la cov est plus negative car le payoff est lineaire
}


void Exo1::mc_anthitetique_sur_M(double &price, double &std_err) {
    int n = 24;
    int M = 50000;
    double somme = 0.0;
    double somme2 = 0.0;
    for (int i = 0; i< M; i ++) {  
        double payoff1 = computeSantitetique1(n) ;
        double payoff2 = computeSantitetique2(n) ;
        double sum_moy = (payoff1+payoff2)*0.5;

        somme += sum_moy;
        somme2 += sum_moy*sum_moy;

    }
    price =  somme/M; // E(X)
    double var =  (somme2 / M) - (price)*(price);
    std_err = sqrt(var/ M);
    cout << "exo1 : MC antithetique sur M  price : " << price << " et IC : " << std_err << endl;

}

void Exo1::mc_anthitetique_sur_M_sur_2(double &price, double &std_err) {
    int n = 24;
    double somme = 0.0;
    double somme2 = 0.0;
    for (int i = 0; i< M/2; i ++) {  // ON PRENDS QUE LA MOITIE
    // on peut faire un classique mc,   avec M
        double payoff1 = computeSantitetique1(n) ;
        double payoff2 = computeSantitetique2(n) ;
        double sum_moy = (payoff1+payoff2);

        somme += sum_moy;
        somme2 += sum_moy*sum_moy*0.5; // ATTENTION ICI ou a mettre en sortir de boucle

    }
    price =  somme/M; // E(X)
    double var =  (somme2 / M) - (price)*(price);
    std_err = sqrt((var) / M);
    cout << "exo1 : MC antithetique sur M/2  price : " << price << " et IC : " << std_err << endl;

}


// VARIABLE DE CONTROLE H(Y) la variable ce controle tp var(f(X) - h(Y)) << var (f(X) )
// E = E(h(Y)) (connu) + E(f(X) - h(Y)) 
// S = E(h(Y)) + 1/n * sum [ E(f(Xi) - h(Yi)) ] converge vers E
// var (S) = 1/ var(f(X)) + [ 1/n * var(h(Y)) - cov (f(X), h(Y))]
// autre ex : C - P = S - K => C = (S - K) + P


// h(Y) = ST - S0*e(rT)
double  Exo1::compute_var_controle1(int nbValues) { 
    double maxS = m_spot;
    double step = m_maturity/nbValues;
    double sqrt_step = sqrt(step);

    //pnl_vect_set(vector_S, 0, m_spot);
    double lastVal = m_spot;

    double bs_drift = m_interest_rate - m_volatility * m_volatility / 2;
    
    for (int i=1; i <= nbValues; i++) {
        double g = pnl_rng_normal(m_rng);
        lastVal *= exp(bs_drift * step + m_volatility * sqrt_step * g);
        //pnl_vect_set(vector_S, i, lastVal);
        if (lastVal > maxS) {
            maxS = lastVal;
        }
    }
    double fxi = m_strike*maxS - lastVal;
    double hxi = lastVal - m_spot* exp(m_interest_rate*m_maturity);
    return fxi-hxi;

}


void Exo1::mc_var_controle1(double &price, double &std_err) {
    int n = 24;
    double somme = 0.0;
    double somme2 = 0.0;
    for (int i = 0; i< M; i ++) { 
        double payoff = compute_var_controle1(n) ;
        somme += payoff;
        somme2 += payoff*payoff;
    }

    // h(Y) = ST - S0exp(rT) mais pa tres interessant
    double esp_h_Y = 0; //(connu) car E(ST) = S0*exp(rT)
    double moyenne_fxi_moin_hyi = somme / M;
    price = esp_h_Y + moyenne_fxi_moin_hyi;

    double var = (somme2 / M) - moyenne_fxi_moin_hyi*moyenne_fxi_moin_hyi;
    std_err = sqrt(var / M);
    cout << "exo1 : MC var_controle1 ( h(Y) = St - S0*e(rT) )  price : " << price << " et IC : " << std_err << endl;

}



// parite call put
// h(Y) = S0 - K car C0 = (S0 - K) + P0 pour les MBG
// mais ici h(Y) = K(maxSt)*e(-rT) - S0) et donc f(x) - h(y) = (ST - KmaxSt)+
double  Exo1::compute_var_controle2(int nbValues, double& maxPath) { 
    double maxS = m_spot;
    double step = m_maturity/nbValues;
    double sqrt_step = sqrt(step);
    double E_var_controle = 0.0;
    double lastVal = m_spot;

    double bs_drift = m_interest_rate - m_volatility * m_volatility / 2;
    
    for (int i=1; i <= nbValues; i++) {
        double g = pnl_rng_normal(m_rng);
        lastVal *= exp(bs_drift * step + m_volatility * sqrt_step * g);
        if (lastVal > maxS) {
            maxS = lastVal;
            E_var_controle += maxS;
        }
    }
    maxPath = maxS;
    return lastVal > m_strike*maxS ? lastVal - m_strike*maxS : 0;
}




void Exo1::mc_var_controle2(double &price, double &std_err) {
    int n = 24;
    double somme = 0.0;
    double somme2 = 0.0;
    double E_max = 0.0;
    for (int i = 0; i< M; i ++) { 
        double maxS;
        double payoff = compute_var_controle2(n, maxS);
        E_max += maxS;
        somme += payoff;
        somme2 += payoff*payoff;
    }
    E_max /=M;
    double esp_h_Y = m_strike*E_max*exp(-m_interest_rate * m_maturity) - m_spot; 
    double moyenne_fxi_moin_hyi = somme / M;
    price = esp_h_Y + moyenne_fxi_moin_hyi;

    // var(In) = var(sum (f(x) -h(x))) = E(f(X) - h(X)^2) - E(f(X) - h(X))^2
    double var = (somme2 / M) - moyenne_fxi_moin_hyi*moyenne_fxi_moin_hyi; 
    std_err = sqrt(var / M);
    cout << "exo1 : MC var_controle2 ( h(Y) = K(maxSt)*e(-rT) - S0)  price : " << price << " et IC : " << std_err << endl;

}



// remarque 
// le calcul d'une integrale :
            // sur [0, 1] peut etre vu comme esperance de f(U)
            // sur R  / // / //                           f(N)
    // ainsi on peut appliquer nos methides de reduction à savoir car :
    // X =  N(0, d) est symetrique don meme loi que -X
    // Y = 1- U a meme loi que U