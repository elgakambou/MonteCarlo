#pragma once

#include "pnl/pnl_vector.h"
#include "pnl/pnl_random.h"
#include "BSBarrier.hpp"


class ImportanceSampling {

    public :    
        BSBarrier myProduct;
        PnlRng* myRng;
        int m_samples;
        int n = 500;
        ImportanceSampling (PnlRng* rng, BSBarrier product, int samples);
        void computeIncrementBrownien(PnlMat * G);
        void computeLambaOpt(double &lambdaOpt);
        double vp(double lambda, PnlMat * Gmat, PnlVect * fcarre);
        double vpp(double lambda, PnlMat * Gmat, PnlVect * fcarre); 
        void importance_sampling_mc (double &prix, double &stddev);
};