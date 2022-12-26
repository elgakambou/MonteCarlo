#include "MultiLevelMonteCarlo.hpp"

MultiLevelMonteCarlo::MultiLevelMonteCarlo(Option *opt, Model *mod, PnlRng *rng)
{
  m_opt = opt;
  m_mod = mod;
  m_rng = rng;
}

long long MultiLevelMonteCarlo::nSamples(int level, int m, int L)
{
    return L * pnl_pow_i(m, 2 * L - level);
}


void MultiLevelMonteCarlo::collapse(PnlMat *Gcrude, PnlMat *Gfine, int m)
{
    long m_l_moins_un = Gfine->m / m;
    double sqrt_m = sqrt(m);
    pnl_mat_resize(Gcrude, m_l_moins_un, 2);
    for (int i=0; i<m_l_moins_un; i++){
        double g1 = 0, g2 = 0;
        for (int j=0; j<m; j++){
            int index = i*m + j;
            g1 += pnl_mat_get(Gfine, index, 0);
            g2 += pnl_mat_get(Gfine, index, 1);
        }
        pnl_mat_set(Gcrude, i, 0, g1/sqrt_m);
        pnl_mat_set(Gcrude, i, 1, g2/sqrt_m);
    }
}

void MultiLevelMonteCarlo::run(double &prix, double &std_dev, int m, int L)
{
    PnlVect *pathCrude = pnl_vect_new();
    PnlVect *pathFine = pnl_vect_new();
    PnlMat *Gcrude = pnl_mat_new();
    PnlMat *Gfine = pnl_mat_new();
    double sum = 0., var = 0., payoff;

    // Treat Level 0
    long N0 = nSamples(0, m, L);
    for (int i=0; i< N0; i++){
        pnl_mat_rng_normal(Gfine, 1, 2, m_rng);
        m_mod->simul(pathFine, m_opt->m_maturity, 1, Gfine);
        payoff = m_opt->payoff(pathFine);
        sum += payoff;
        var += payoff * payoff;
    }
    sum /= N0;
    var = (var / N0 - sum * sum) / N0;
    prix = sum;
    std_dev = sqrt(var);

    // Loop on all the levels > 0
    for (int l = 1; l < L; l++) {
        long long Nl = nSamples(l, m, L);
        int ml = pnl_pow_i(m, l);
        sum = 0; var = 0;
        for (int i = 0; i < Nl; i++) {
            pnl_mat_rng_normal(Gfine, ml, 2, m_rng);

            // Simulation of the fine model
            m_mod->simul(pathFine, m_opt->m_maturity, ml, Gfine);

            // Simulation of the crude model
            collapse(Gcrude, Gfine, m);
            m_mod->simul(pathCrude, m_opt->m_maturity, ml / m, Gcrude);

            // Payoffs of both models
            payoff = m_opt->payoff(pathFine);
            payoff -= m_opt->payoff(pathCrude);
            // adding both payoffs to the sum
            sum += payoff;
            var += payoff * payoff;
        }

        sum /= Nl;
        var = (var / Nl - sum * sum) / Nl;
        prix += sum;
        std_dev += sqrt(var/Nl);
    }

    // Actualisation de nos prix et variance
    prix *= exp(- m_mod->m_r * m_opt->m_maturity);
    std_dev *= exp(- m_mod->m_r * m_opt->m_maturity);

    pnl_mat_free(&Gcrude);
    pnl_mat_free(&Gfine);
    pnl_vect_free(&pathCrude);
    pnl_vect_free(&pathFine);
}