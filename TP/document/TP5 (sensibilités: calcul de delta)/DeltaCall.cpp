#include "DeltaCall.hpp"
#include <iostream>
#include <cmath>

    DeltaCall::DeltaCall (double spot, double sigma, double r, double T, double K) :
    m_S0(spot), m_sigma(sigma), m_r(r), m_T(T), m_K(K) {}

    double DeltaCall:: payoffCall(double S0,double sigma, double r, double T, double K, double g) {

            double ST = S0*exp((r- sigma*sigma/2)*T + sigma*sqrt(T)*g);
            return ST > K ? (ST - K)  : 0.0;

    }

     double DeltaCall:: derive_payoffCall(double S0,double sigma, double r, double T, double K, double g) {

            double ST_prime_S0 = exp((r- sigma*sigma/2)*T + sigma*sqrt(T)*g);
            double ST = S0*ST_prime_S0;
            return ST > K ? ST_prime_S0  : 0.0;
    }


    // simple,  toujours possible mais biaisé à cause de epsilon
    void DeltaCall :: delta_df(double&price, double &std_err, double epsilon, PnlRng* rng) {

        double g;
        double somme_plus_eps = 0.0;
        double somme_minus_eps = 0.0;

        for (int i = 0; i < M; i++) {
            g = pnl_rng_normal(rng);
            // on garde la meme trajectoire de brownien
            somme_plus_eps += payoffCall(m_S0 + epsilon, m_sigma, m_r, m_T, m_K, g);
            somme_minus_eps += payoffCall(m_S0 - epsilon, m_sigma, m_r, m_T, m_K, g);
        }

        somme_plus_eps /= M;
        somme_minus_eps /=M;
        double delta  = exp(-m_r*m_T) * (somme_plus_eps - somme_minus_eps) / (2*epsilon); // actualsation
        std::cout << "delta_df delta : " << delta << "\n";
    }

    // ici on derive la densité de Y_theta
    // precis et probleme de derivabilite comme pathwise
    // densité a connaitre
    void DeltaCall :: delta_likelihood (double&price, double &std_err, PnlRng* rng) { // voir cours pour la formaul dans le cas BS
        
        double somme = 0.0;
        double WT ;
        double g ;
        for (int i = 0; i < M; i++) {
            g = pnl_rng_normal(rng);
            // on garde la meme trajectoire de brownien
            WT = sqrt(m_T)*g;
            somme += payoffCall(m_S0, m_sigma, m_r, m_T, m_K, g) * WT;
        }
        somme /=M; 
        somme /= m_sigma*m_T*m_S0;
        somme *= exp(-m_r*m_T); // actualisation
        std::cout << "delta_likelihood delta : " << somme << "\n";
    }
    // ici on derive la fonction payoff
    // non interessant si on a pas de formule fermé et quand le payoff n'est pas derivable 
    void DeltaCall :: delta_pathwise (double&price, double &std_err, PnlRng* rng) { // plus simple et non biaisé

        double somme = 0.0;
        double g ;
        for (int i = 0; i < M; i++) {
            g = pnl_rng_normal(rng);
            somme += derive_payoffCall(m_S0, m_sigma, m_r, m_T, m_K, g); // derive du payoff par rapport a theta, ici thetta = S0
        }
        somme /=M; 
        somme *= exp(-m_r*m_T); // actualisation
        std::cout << "delta_likelihood delta : " << somme << "\n";
    }