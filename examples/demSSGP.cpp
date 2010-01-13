#include <iostream>

#include "itpp/itbase.h"
#include "itpp/itstat.h"
#include "io/csvstream.h"

#include "covarianceFunctions/SumCovarianceFunction.h"
#include "covarianceFunctions/GaussianCF.h"
#include "covarianceFunctions/WhiteNoiseCF.h"

#include "optimisation/SCGModelTrainer.h"

#include "likelihoodModels/GaussianLikelihood.h"

#include "gaussianProcesses/SequentialGP.h"
#include "gaussianProcesses/GaussianProcess.h"

#include "GraphPlotter/GraphPlotter.h"

using namespace itpp;

int main() {
    
    vec Xtrn, Xtst, Ytrn, Ytst;
    vec gpmean, gpvar, ssgpmean, ssgpvar;
    
    GraphPlotter gplot = GraphPlotter();
    
    // Generate some data from a GP
    double range  = 2.0;
    double sill   = 2.0;
    double nugget = 0.05;
    

    // Covariance function: Gaussian + Nugget
    GaussianCF   g1(range, sill);            
    WhiteNoiseCF g2(nugget);
    SumCovarianceFunction gCmp(g1);
    gCmp.addCovarianceFunction(g2);
    
    // Generate some data from a GP
    int n = 200;
    Xtst = 15.0*(randu(n)-0.5);
    sort(Xtst);
    mat Xtstmat = Xtst;
    mat K = zeros(Xtst.length(), Xtst.length());
    gCmp.computeSymmetric(K, Xtstmat);
    Ytst = (chol(K)).transpose()*randn(n);
    
    ivec itrn = to_ivec(linspace(0,Xtst.length()-1,40));
    cout << itrn << endl;
    Xtrn = Xtst(itrn);
    Ytrn = Ytst(itrn);
    mat Xtrnmat = Xtrn;
    
    // gplot.plotPoints(Xtst, Ytst, "true function", LINE, RED);
        
    
    // SSGP
    int n_active = 2;
    
    SequentialGP ssgp(1, 1, n_active, Xtrnmat, Ytrn, gCmp);
            
    // Gaussian observation likelihood
    GaussianLikelihood gaussLik(nugget);
    // ssgp.computePosterior(gaussLik);
          
    for (int i=0; i<8; i++) {
        cout << "SSGP with " << n_active << " active points" << endl;
        
        ivec iActive;
        if (n_active == Xtrn.length())
            iActive = to_ivec(linspace(0,Xtrn.length()-1,n_active));
        else {
            iActive = to_ivec(linspace(0,Xtrn.length()-1,n_active+2));
            iActive.del(iActive.length()-1);
            iActive.del(0);    
        }
        
                
        // ssgp.resetPosterior();
        ssgp.computePosteriorFixedActiveSet(gaussLik, iActive);
        
        // ssgp.computePosterior(gaussLik);
        ssgp.makePredictions(ssgpmean, ssgpvar, Xtst, g1);
        
        gplot.clearPlot();
        
        // Plot SSGP mean and error bars
        gplot.plotPoints(Xtst, ssgpmean, "ssgp mean", LINE, BLUE); // SSGP
        gplot.plotPoints(Xtst, ssgpmean + 2.0*sqrt(ssgpvar), "error bar", LINE, CYAN);
        gplot.plotPoints(Xtst, ssgpmean - 2.0*sqrt(ssgpvar), "error bar", LINE, CYAN);
        
        // Plot observations
        gplot.plotPoints(Xtrn, Ytrn, "training set", CROSS, RED);  
                
        // Plot active points
        vec activeX = (ssgp.getActiveSetLocations()).get_col(0);
        vec activeY = Ytrn(ssgp.getActiveSetIndices());
        gplot.plotPoints(activeX, activeY, "active points", CIRCLE, BLUE);
        
        mat Xtrnmatgp;
        
        iActive = ssgp.getActiveSetIndices();
        Xtrnmatgp = Xtrn(iActive);
        vec Ytrngp = Ytrn(iActive);
        
        GaussianProcess gp(1, 1, Xtrnmatgp, Ytrngp, gCmp);
        
        vec gpmean, gpvar;
        gp.makePredictions(gpmean, gpvar, Xtst, g1);
        
        // Plot GP mean and error bars
        gplot.plotPoints(Xtst, gpmean, "gp mean", LINE, GREEN); // SSGP
        gplot.plotPoints(Xtst, gpmean + 2.0*sqrt(gpvar), "gp error bar", LINE, GREEN);
        gplot.plotPoints(Xtst, gpmean - 2.0*sqrt(gpvar), "gp error bar", LINE, GREEN);
                
        
        n_active += 1;
        ssgp.setActiveSetSize(n_active);
        cout << "Press a key to continue" << endl;
        getchar();
    }
    
    
    return 0;
    
}
