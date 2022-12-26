#include "Heston.hpp"
#include "iostream"
#include <cmath>
#include <pnl/pnl_mathtools.h>

HestonModel::HestonModel(double r, double spot, double initVol, double sigma, double kappa, double theta, double rho)
  : Model(2, r)
  , m_spot(spot)
  , m_initVol(initVol)
  , m_sigma(sigma)
  , m_kappa(kappa)
  , m_theta(theta)
  , m_rho(rho)
{}

// g1 et g2 sont iid N(0,1)
void HestonModel::compute_W1_W2_correled (double& w1, double& w2, double rho, double g1, double g2) {
  w1 = g1;
  w2 =  rho* g1 + sqrt(1-rho*rho) * g2 ;
}



void HestonModel::simul(PnlVect *path, double maturity, double nTimeSteps, PnlMat *G)
{
  pnl_vect_resize(path, nTimeSteps + 1);

  pnl_vect_set(path, 0, m_spot);
  double current_S = m_spot;
  double vol = m_initVol > 0 ? m_initVol : 0.;
  double step = maturity / nTimeSteps;

  double dw1, dw2;
  for (int i = 1; i < nTimeSteps + 1  ; i++) {
    compute_W1_W2_correled(dw1, dw2, m_rho, pnl_mat_get(G, i-1, 0), pnl_mat_get(G, i-1, 1));

    //std::cout << dw1 << " " << dw2 << " " <<  current_S <<"\n";
    //std::cout << "oui" << "\n";
    dw1 = dw1*sqrt(step);
    dw2 = dw1*sqrt(step);
    current_S = current_S + m_r*current_S*step + sqrt(vol) * current_S *  dw1 ;
    vol = vol +  m_kappa * ( m_theta - vol) * step + m_sigma * sqrt(vol) * dw2 ;
    vol  = vol > 0 ? vol : 0.0;
    pnl_vect_set(path, i, current_S);
  }
  // pnl_vect_resize(path, nTimeSteps + 1);
  //   double dt = maturity / nTimeSteps;
  //   double sqrt_dt = sqrt(dt);
  //   double sqrt_one_minus_square_rho = sqrt(1- m_rho * m_rho);
  //   double path_k = m_spot;
  //   double v_k = m_initVol;

  //   LET(path, 0) = path_k;
  //   for (int k=1; k<=nTimeSteps; k++){
  //       double w2 = sqrt_one_minus_square_rho * pnl_mat_get(G, k-1, 1) + m_rho * pnl_mat_get(G, k-1, 0);
  //       path_k += dt * m_r * path_k + sqrt(MAX(v_k, 0.)) * path_k * sqrt_dt * pnl_mat_get(G, k-1, 0);
  //       v_k += dt * m_kappa * (m_theta - v_k) + m_sigma * sqrt(MAX(v_k, 0.)) * sqrt_dt * w2; // TODO : Reponse : vk après sinon vk n'est plus mesurable
  //       LET(path, k) = path_k;
  //   }
}