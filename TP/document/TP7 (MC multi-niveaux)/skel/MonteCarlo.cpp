#include "MonteCarlo.hpp"
#include "iostream"

MonteCarlo::MonteCarlo(Option *opt, Model *mod, PnlRng *rng)
{
  m_opt = opt;
  m_mod = mod;
  m_rng = rng;
}




void MonteCarlo::run(double &prix, double &std_dev, int nSamples, int nTimeSteps)
{
    double somme = 0.;
    double somme2 = 0.;
    double payoff = 0.0;

    PnlVect * path = pnl_vect_create (nTimeSteps + 1);

    // martrice des N(0,1)
    PnlMat *G = pnl_mat_create (nTimeSteps, 2);
    for (int i = 0; i < nSamples; i ++) {
      pnl_mat_rng_normal(G, nTimeSteps, 2, m_rng);
      m_mod->simul(path, m_opt->m_maturity, nTimeSteps, G);
      payoff = m_opt->payoff(path);
      somme += payoff;
      somme2 += payoff*payoff;
    }

    prix = (somme / nSamples) * exp(-m_opt->m_maturity * m_mod->m_r);
    std_dev = (somme2/nSamples) - prix*prix;
    std_dev = sqrt(std_dev/nSamples)* exp(- m_mod->m_r * m_opt->m_maturity);
    pnl_vect_free(&path);
    pnl_mat_free(&G);
}

void MonteCarlo::mse(double &prix, double &std_dev, double &mse, long long nSamples, int nTimeSteps){
    run(prix, std_dev, nSamples, nTimeSteps);
    mse = std_dev * std_dev + (prix - 3.847906) * (prix - 3.847906);

};