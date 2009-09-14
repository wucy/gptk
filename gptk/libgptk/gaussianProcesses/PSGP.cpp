#include "PSGP.h"

/**
 * Constructor
 * 
 * @params:
 * 
 * X             Matrix of inputs (locations)
 * Y             Vector outputs (observations)
 * nActivePoints Maximum number of active points 
 * 
 */
PSGP::PSGP(mat& X, vec& Y, CovarianceFunction& cf, int nActivePoints, int _iterChanging, int _iterFixed) 
: ForwardModel(X.cols(), 1), Locations(X), Observations(Y), covFunc(cf)
{
    assert(Locations.rows() == Observations.size());

    resetPosterior();
    
    maxActiveSet = nActivePoints;
    epsilonTolerance = 1e-6;
    gammaTolerance = 1e-3;
    
    
    momentProjection = true;
    
    iterChanging = _iterChanging;
    iterFixed = _iterFixed;

    nObs = Locations.rows();

    likelihoodType = Approximate;
}


/**
 * Destructor
 */
PSGP::~PSGP()
{
}


/**
 * Compute posterior with a single likelihood model
 */
void PSGP::computePosterior(const LikelihoodType& noiseModel)
{
    bool fixActiveSet = false;

    // Cycle several times through the data, first allowing the active
    // set to change (for iterChanging iterations) and then fixing it
    // (for iterFixed iterations)
    for(int cycle = 1; cycle <= (iterChanging + iterFixed); cycle++)
    {
        if(cycle > iterChanging) fixActiveSet = true;
        
        // Present observations in a random order
        ivec randObsIndex = itppext::randperm(nObs);
        
        // DEBUG START
        // ivec randObsIndex = to_ivec(linspace(0,nObs-1,nObs));
        // DEBUG END
        
        for(int i=0; i<nObs; i++)	
        {
            cout << "\rProcessing observation: " << i+1 << "/" << nObs  << flush;
            processObservationEP(randObsIndex(i), noiseModel, fixActiveSet);
        }
        cout << endl;
    }
}


/**
 * Compute posterior with a different likelihood model for each observation
 */
void PSGP::computePosterior(const ivec& modelIndex, const Vec<LikelihoodType *> noiseModel)
{
    // Check if we have an index and a model per observation
    assert(nObs == modelIndex.length()); 
    assert(nObs == noiseModel.size());

    bool fixActiveSet = false;

    for(int cycle = 1; cycle <= (iterChanging + iterFixed); cycle++)
    {
        if(cycle > iterChanging) fixActiveSet = true;

        // Present observations in a random order
        ivec randObsIndex = itppext::randperm(nObs);

        for(int iObs=0; iObs<nObs; iObs++)	
        {
            int iModel = modelIndex(randObsIndex(iObs));
            cout << "\rProcessing observation: " << iObs+1 << "/" << nObs;
            cout << " using likelihood model " << iModel << flush;
            
            assert(iModel < noiseModel.length() && iModel >= 0);
            
            processObservationEP(iObs, *noiseModel(iModel), fixActiveSet);
        }
        cout << endl;
    }
}


/**
 * This is the core method, implementing the sparse EP algorithm
 * 
 * This algorithm follows Appendix G in in Lehel Csato's PhD thesis 
 * "Gaussian Processes - Iterative Sparse Approximations", 2002,
 * NCRG, Aston University.
 * 
 * Further reference numbers correspond to sections/equations in the same
 * document. 
 */
void PSGP::processObservationEP(const int iObs, const LikelihoodType &noiseModel, const bool fixActiveSet) 
{
    double sigmaLoc;                // Auto-covariance of location
    vec k = zeros(sizeActiveSet);   // Covariance between location and active set
    
    double r, q;                    // Online coefficients
    double cavityMean, cavityVar;   // Cavity mean and variance
    double gamma;                   // Mean of ??
    vec    eHat;                    // Variance of ??
    
    // Retrieve location and observation at index iObs
    vec loc    = Locations.get_row(iObs);
    double obs = Observations(iObs);
    
    // Remove previous contribution of observation iObs
    // Appendix G.(a)
    EP_removePreviousContribution(iObs);
    
    // Intermediate computations: covariances, cavity mean and variance
    // Appendix G.(c)
    EP_updateIntermediateComputations(cavityMean, cavityVar, sigmaLoc, k, gamma, eHat, loc);
    
    // Compute the updated q and r online coefficients for the specified 
    // likelihood model, using the updated alpha and C.
    // Appendix G.(b), Sec. 2.4, Eq. 2.42, 3.28
    double logEvidence = noiseModel.updateCoefficients(q, r, obs, cavityMean, cavityVar);
    
    stabiliseCoefficients(q, r, cavityMean, cavityVar, 1e3, 1e-6);
    
    // Update EP parameters 
    // Appendix G.(d)
    meanEP(iObs) = cavityMean - q/r;
    varEP(iObs) = -r / (1.0 + r*cavityVar);
        
    // Scaling factor (depends on update type, sparse or full)
    double eta;
    
    vec s;
    if (sizeActiveSet>0) {
        s = C*k;
    }
    
    // Perform full or sparse update depending on geometry (gamma)
    // Appendix G.(e)
    if (gamma >= gammaTolerance*sigmaLoc && !fixActiveSet)
    {
        cout << " (full update)";
        // Full update - add location to active set
        addActivePoint(iObs);

        // EP scaling factor
        eta = 1.0;
        
        // Update S matrix (Appendix D.2)
        // vec z = zeros(S.rows());
        // S.append_col(z);
        // double snew = 1.0/(r*eta) + dot(k,s) + sigmaLoc;
        // S.append_row( concat(z, 1.0/snew));
        
        // Increase size of s
        s.set_size(sizeActiveSet, true);
        s(sizeActiveSet-1) = 1.0; 
        
        // e is the unit vector for dimension sizeActiveSet
        vec e = zeros(nObs);
        e(iObs) = 1.0;
        
        // Update P matrix
        P.append_col(e);
        
        // Update KB matrix - add auto and cross covariances
        KB.append_col(k);
        KB.append_row( concat(k, sigmaLoc) );
        
        // update Q matrix
        eHat.set_size(sizeActiveSet, true);
        eHat(sizeActiveSet - 1) = -1.0;
        Q.set_size(sizeActiveSet, sizeActiveSet, true);
        Q += outer_product(eHat, eHat) / gamma;
    }
    else
    {
        cout << " (sparse update)";
        // Sparse update
        s += eHat;
        P.set_row(iObs, eHat);
        
        eta = 1.0 / (1.0 + gamma*r);
    }
    
    // update EP parameters (Eq. 4.18)
    double ratio = q / r;
    logZ(iObs)= logEvidence + (log(2.0 * pi) - log(abs(r)) - (q * ratio)) / 2.0;
    meanEP(iObs) = cavityMean - ratio;
    varEP(iObs) = -r / (1.0 + (r * cavityVar));
    
    // Update GP parameters
    // Appendix G.(f)
    alpha += eta * s * q;
    C     += r * eta * outer_product(s, s, false);

    // Appendix G.(g)
    if (sizeActiveSet > maxActiveSet) 
    {
        // Compute scores for each active point
        vec scores = scoreActivePoints(FullKL);
        
        // Remove active point with lowest score, and resize/update
        // alpha, C, Q and P
        int removalCandidate = min_index(scores);
        deleteActivePoint(removalCandidate);
    }
    
    // Remove unneeded points based on geometry
    EP_removeCollapsedPoints();
    

}


/**
 * Substract contribution of observation iObs from previous iterations 
 * 
 * Appendix G.(a), Eq. 4.19 and 4.20
 */
void PSGP::EP_removePreviousContribution(int iObs) 
{
    if (varEP(iObs) > LAMBDA_TOLERANCE)   
    {
        vec p = P.get_row(iObs);
        vec Kp = KB * p;

        // Update alpha and C
        vec h  = C * Kp + p;
        double nu = varEP(iObs) / ( 1.0 - varEP(iObs) * dot(Kp,h) );
        alpha += h * nu * ( dot(alpha, Kp) - meanEP(iObs) );
        C += nu * outer_product(h,h);
    }
}


/**
 * Compute cavity mean and variance for specified input/location x
 * 
 * Appendix G.(c), Eq. 3.3
 */
void PSGP::EP_updateIntermediateComputations(double &cavityMean, double &cavityVar, double &sigmaLoc,
                                              vec &k, double &gamma, vec &eHat, vec loc) 
{
    assert(k.length() == sizeActiveSet);
    
    covFunc.computeSymmetric(sigmaLoc, loc);           // Auto-variance of location
    
    if(sizeActiveSet == 0)
    {
        cavityVar = sigmaLoc;
        cavityMean = 0.0;
        eHat = zeros(0);
        gamma = sigmaLoc;
    } 
    else 
    {
        covFunc.computeSymmetric(sigmaLoc, loc);           // Auto-variance of location
        covFunc.computeCovariance(k, ActiveSet, loc);      // cov(location, active set)
        cavityVar = sigmaLoc + dot(k, C * k);
        cavityMean = dot(k, alpha);
        eHat = backslash(KB, k);
        gamma = sigmaLoc - dot(k, eHat);
    }
}


/**
 * Add a given location to active set
 */
void PSGP::addActivePoint(int iObs)
{
    // Increase size of active set
    sizeActiveSet++; 
   
    // Append index of observation to indexes in active set
    idxActiveSet.set_size(sizeActiveSet, true);
    idxActiveSet(sizeActiveSet - 1) = iObs;

    // Increase storage size for active set and store new observation
    // ActiveSet.set_size(sizeActiveSet, Locations.cols(), true);
    // ActiveSet.set_row(sizeActiveSet - 1, Locations.get_row(iObs));
    ActiveSet.append_row(Locations.get_row(iObs));

    // Increase size of C and alpha
    alpha.set_size(sizeActiveSet, true);
    C.set_size(sizeActiveSet, sizeActiveSet, true);
}


/**
 * Delete an active point from the active set
 *
 * @params:
 *  iObs    The index of the active point to be deleted
 */
void PSGP::deleteActivePoint(int iObs)
{

    // Elements for iObs (correspond to the * superscripts in Csato)
    double alpha_i = alpha(iObs);
    double c_i     = C(iObs, iObs);
    double q_i     = Q(iObs, iObs);
    vec    P_i     = P.get_col(iObs);

    // Covariance between iObs and other active points
    vec    C_i     = C.get_row(iObs);
    vec    Q_i     = Q.get_row(iObs);
    C_i.del(iObs);  // Delete cov(iObs, iObs), we only 
    Q_i.del(iObs);  // want the cross terms
    
    // Updated elements without iObs (correspond to the "r" superscripts in Csato)
    alpha.del(iObs);
    C.del_col(iObs);
    C.del_row(iObs);
    Q.del_col(iObs);
    Q.del_row(iObs);
    P.del_col(iObs);
    
    // Update new (reduced) elements
    alpha -= ( alpha_i / (c_i + q_i) ) * (Q_i + C_i);    // Eq. 3.19
    mat QQq = outer_product( Q_i, Q_i ) / q_i;
    C += QQq - outer_product( Q_i+C_i, Q_i+C_i ) / ( q_i+c_i );
    Q -= QQq;
    P -= outer_product( P_i, Q_i) / q_i;
    
    // Update S (Appendix D.2)
    // Not used - too computationally expensive
    // See workaround in scoreActivePoints
    /*
    S = KB - KB*S*KB;
    S.del_col(iObs);
    S.del_row(iObs);
    S = Q - Q*S*Q;
    */
    
    // Update active set
    KB.del_row(iObs);
    KB.del_col(iObs);
    
    ActiveSet.del_row(iObs);
    idxActiveSet.del(iObs);
    sizeActiveSet--;
}

/**
 * addOne - split into parts
 *
 */
void PSGP::EP_removeCollapsedPoints()
{
    while(sizeActiveSet > 0)
    {
        vec scores = scoreActivePoints(Geometric);
        int removalCandidate = min_index(scores);
        
        if(scores(removalCandidate) >= (gammaTolerance / 1000.0))
        {
            break;
        }
        deleteActivePoint(removalCandidate);
    }    
}


/**
 * Score active points according to scoring method
 */
vec PSGP::scoreActivePoints(ScoringMethod sm)
{
    vec diagC, diagS, term1, term2, term3;
    vec diagInvGram = diag(Q);
        
    switch(sm)
    {
    
    case Geometric : 
        return(1.0 / diagInvGram);
        break;
    
    case MeanComponent : 
        return(elem_div(elem_mult(alpha,alpha), diag(C) + diagInvGram));
        break;

    case FullKL : // Lehel: Eq. 3.23 
        diagC = diag(C);
        
        // Computation of the S matrix is too expensive - use P instead
        // diagS = -diag(S);
        // 
        // We are trying to compute the expression below:
        // diagS = diag((P.transpose() * diag(varEP)) * P);
        //
        // This form is faster:
        diagS = zeros(P.cols());
        for (int i=0; i<P.cols(); i++) {
          diagS(i) = elem_mult_sum(varEP, elem_mult(P.get_col(i),P.get_col(i)));
        }
        // cout << "|diagS2 - diagS| = " << sum(abs(diagS2 - diagS)) << endl;
        // getchar();
        
        term1 = elem_div(elem_mult(alpha,alpha), diagC + diagInvGram);
        term2 = elem_div(diagS, diagInvGram);
        term3 = log(1.0 + elem_div(diagC , diagInvGram));
        return (term1 + term2 + term3);
        break;
    
    default : 
        cerr << "Unknown scoring method" << endl;
        break;
    }
    return zeros(sizeActiveSet);
}


/**
 * Check update coefficients to ensure stability
 */
void PSGP::stabiliseCoefficients(double& q, double& r, double cavityMean, double cavityVar, double upperTolerance, double lowerTolerance)
{
    double sqrtPt = sqrt(cavityVar);
    double tu = -sqrtPt * r * sqrtPt; 
    bool mod = false;
    if(tu > upperTolerance)
    {
        tu = upperTolerance;
        cout << "Shrinking lambda" << endl;
        mod = true;
    }

    if(tu < lowerTolerance)
    {
        tu = lowerTolerance;
        cout << "(REV) Shrinking lambda" << endl;
        mod = true;
    }

    if(mod) {
        r = - (tu / sqrtPt) / tu;
        r = r + itpp::eps;
        r = r + r;
    }
}


/**
 * Make predictions
 */
void PSGP::makePredictions(vec& Mean, vec& Variance, const mat& Xpred, CovarianceFunction& cf) const
{
    assert(Mean.length() == Variance.length());
    assert(Xpred.rows() == Mean.length());

    // Predictive mean
    mat ktest(Xpred.rows(),sizeActiveSet); 
    cf.computeCovariance(ktest, Xpred, ActiveSet);
    Mean = ktest*alpha;
    
    // Predictive variance
    vec kstar(Xpred.rows());
    covFunc.computeDiagonal(kstar, Xpred);
    Variance = kstar + sum(elem_mult((ktest * C), ktest), 2);
}


/**
 * Same as above, but using the current (stored) covariance function to make
 * the predictions.
 **/
void PSGP::makePredictions(vec& Mean, vec& Variance, const mat& Xpred) const
{
    makePredictions(Mean, Variance, Xpred, covFunc);
}


/**
 * Simulate from PSGP
 */
vec PSGP::simulate(const mat& Xpred, bool approx) const
{
    mat cov, vCov, kxbv(Xpred.rows(), sizeActiveSet);
    vec dCov, samp;

    covFunc.computeCovariance(kxbv, Xpred, ActiveSet);

    if(approx)
    {
        cov = Q + C;
        vCov.set_size(sizeActiveSet, sizeActiveSet);
        dCov.set_size(sizeActiveSet);
        eig_sym(cov, dCov, vCov);
        vCov = kxbv * vCov;
        samp = randn(sizeActiveSet);
    }
    else
    {
        mat kxx(Xpred.rows(), Xpred.rows());
        covFunc.computeSymmetric(kxx, Xpred);
        cov = kxx + ((kxbv * C) * kxbv.transpose());
        eig_sym(cov, dCov, vCov);
        samp = randn(Xpred.rows());
    }

    dCov = sqrt(abs(dCov));

    vec a1 = kxbv * alpha;
    mat a2 = vCov * diag(dCov);

    return(a1 + (a2 * samp));
}


/**
 * Get covariance function parameters
 */
vec PSGP::getParametersVector() const
{
    return covFunc.getParameters();
}


/**
 * Set covariance function parameters
 */
void PSGP::setParametersVector(const vec p)
{
    //	cout << "Set parameters vector"<< endl;
    covFunc.setParameters(p);
}


/**
 * Recompute posterior parameters
 */
void PSGP::recomputePosterior()
{
    mat KBold = KB;
    mat Kplus(Observations.length(), sizeActiveSet);
    covFunc.computeSymmetric(KB, ActiveSet);
    covFunc.computeCovariance(Kplus, Locations, ActiveSet);
    
    // P should be transpose(inv(KB)*Kplus), not inv(KB)*Kplus (size is wrong otherwise)   
    mat Ptrans(P.cols(), P.rows());
    backslash(KB, Kplus.transpose(), Ptrans); // output is rightmost parameter
    P = Ptrans.transpose();

    mat varEPdiag = diag(varEP);
    mat projLam = P.transpose() * varEPdiag;
    mat UU = projLam * P;
    mat CC = UU * KB + eye(sizeActiveSet);
    alpha = backslash(CC, projLam * meanEP);
    C = -backslash(CC, UU);
    Q = computeInverseFromCholesky(KB);
}

/**
 * RB: Note sure what the point of this is. I thint I added it at some point,
 *     but really ought to get rid of it.  
 */
void PSGP::updateModel() 
{
    recomputePosterior();
}

/**
 * Reset posterior representation
 */
void PSGP::resetPosterior()
{
    C.set_size(0, 0);
    KB.set_size(0, 0);
    Q.set_size(0, 0);
    // S.set_size(0,0);
    alpha.set_size(0);
    ActiveSet.set_size(0, getInputDimensions());
    idxActiveSet.set_size(0);
    P = zeros(Observations.length(), 0);
    varEP = zeros(Observations.length());
    meanEP = zeros(Observations.length());
    logZ = zeros(Observations.length());
    sizeActiveSet = 0;
}

/**
 * Objective function
 * 
 * RB: Can we possibly replace the switch() with an OO version? I.e.
 * have an Evidence class which is extended by EvidenceFull, 
 * EvidenceApproximate and EvidenceUpperBound? This would also solve the issue
 * of an unknown evidence model.
 */
double PSGP::objective() const
{
    double evidence;

    switch(likelihoodType)
    {
    case FullEvid : 
        evidence = compEvidence();
        break;
    
    
    case Approximate : 
        evidence = compEvidenceApproximate();
        break;
    
    case UpperBound :	
        evidence = compEvidenceUpperBound();
        break;
        
    default:	
        // RB: This really ought to throw an exception
        cerr << "Error in PSGP::objective: Unknown likelihood type." << endl;
        return 0.0;
    }
    return evidence;
}

/**
 * Gradient of the objective function
 *
 * RB: Can we possibly replace the switch() with an OO version? I.e.
 * have an Evidence class which is extended by EvidenceFull, 
 * EvidenceApproximate and EvidenceUpperBound? Would also solve the issue
 * of an unknown evidence model.
 
 */
vec PSGP::gradient() const
{
    vec g;

    switch(likelihoodType)
    {
    case FullEvid : 
        g = gradientEvidence();
        break;
        
    case Approximate : 
        g = gradientEvidenceApproximate();
        break;
        
    case UpperBound :	
        g = gradientEvidenceUpperBound();
        break;
        
    default : 
        // RB: This really ought to throw an exception
        g.set_size(covFunc.getNumberParameters());
    }

    return g;
}

/**
 * Full evidence for current covariance function
 */
double PSGP::compEvidence() const
{

    cvec es;
    mat KB_new(sizeActiveSet, sizeActiveSet);

    covFunc.computeSymmetric(KB_new, ActiveSet);

    double evid = sum(log(varEP));

    evid -= sum(elem_mult(pow(meanEP, 2.0), varEP));
    evid += 2.0 * sum(logZ);
    evid -= varEP.length() * log(2.0 * pi);
    mat Klp = P.transpose() * diag(varEP);
    mat Ksm = (Klp * P) * KB_new + eye(sizeActiveSet);
    vec Kall = Klp * meanEP;	
    mat Kinv= backslash(Ksm.transpose(), KB_new.transpose());

    evid += dot(Kall, Kinv.transpose() * Kall);

    if(eig(Ksm, es))
    {
        evid -= sum(log(real(es)));
    }
    else
    {
        cerr << "Problem computing evidence: eig_sum()" << endl;
    }

    return -evid / 2.0;
}

/**
 * Approximate evidence for current covariance function
 */
double PSGP::compEvidenceApproximate() const
{
    mat cholSigma(sizeActiveSet, sizeActiveSet);
    mat Sigma(sizeActiveSet, sizeActiveSet);
    
    covFunc.computeSymmetric(Sigma, ActiveSet);
    mat invSigma = computeInverseFromCholesky(Sigma);
    vec obsActiveSet = Observations(idxActiveSet);
    
    vec alpha = invSigma * obsActiveSet;

    double like1 = sum(log(diag(chol(Sigma))));
    double like2 = 0.5 * dot(obsActiveSet, alpha);

    return like1 + like2 + 0.5 * sizeActiveSet * log(2 * pi);
}

/**
 * Upper bound on the evidence for current covariance function
 */
double PSGP::compEvidenceUpperBound() const
{
    mat KB_new(sizeActiveSet, sizeActiveSet);
    covFunc.computeSymmetric(KB_new, ActiveSet);
    
    // cout << "KB_new = " << KB_new << endl;
    
    // cout << "KB=" << endl << KB << endl;
    // cout << "KB_new=" << endl << KB_new << endl;
    
    /*
	double cond = itppext::cond(KB_new);
	while (cond < 1e-6) {
	    KB_new = KB_new + 1e-6*eye(KB_new.rows());
	    cout << "Ill-conditionned matrix (cond = " << cond << ")" << endl;
	    cond = itppext::cond(KB_new);
	}
     */

    mat U(KB_new.rows(), KB_new.cols());
    if (!chol(KB_new, U)) { 
        cout << "Error in Cholesky decomposition of KB_new" << endl;
    }
    
    // cout << "sum(log(eig(KBnew))) = " << 2.0*sum(log(diag(U))) << endl;
    double like1 = 2.0 * (sum(log(diag(U))));
    double like2 = trace((eye(sizeActiveSet) + 
            (KB * (C + outer_product(alpha, alpha)))) * backslash(KB_new, KB));
    // cout << "lik1 = " << like1 << endl;
    // cout << "lik2 = " << like2 << endl;
    
    like2 = trace((eye(sizeActiveSet) + 
                (KB * (C + outer_product(alpha, alpha)))) * backslash(U,backslash(U.transpose(), KB)) );
    // cout << "lik2 = " << like2 << endl;
    
    return like1 + like2;
}

/**
 * Gradient of full evidence
 */
vec PSGP::gradientEvidence() const
{
    vec grads = zeros(covFunc.getNumberParameters());
    return grads;

}

/**
 * Gradient of approximate evidence
 */
vec PSGP::gradientEvidenceApproximate() const
{
    vec grads(covFunc.getNumberParameters());

    mat cholSigma(sizeActiveSet, sizeActiveSet);
    mat Sigma(sizeActiveSet, sizeActiveSet);

    covFunc.computeSymmetric(Sigma, ActiveSet);
    cholSigma = computeCholesky(Sigma);
    mat invSigma = computeInverseFromCholesky(Sigma);
    vec obsActiveSet = Observations(idxActiveSet);
    vec alpha = invSigma * obsActiveSet;

    mat W = (invSigma - outer_product(alpha, alpha, false));

    mat partialDeriv(sizeActiveSet, sizeActiveSet);

    for(int i = 0; i < covFunc.getNumberParameters(); i++)
    {
        covFunc.getParameterPartialDerivative(partialDeriv, i, ActiveSet);
        grads(i) = elem_mult_sum(W, partialDeriv) / 2.0;
    }
    return grads; 
}

/**
 * Gradient of upper bound on evidence
 */
vec PSGP::gradientEvidenceUpperBound() const
{

    vec grads(covFunc.getNumberParameters());

    mat W = eye(sizeActiveSet);
    mat KB_new(sizeActiveSet, sizeActiveSet);
    covFunc.computeSymmetric(KB_new, ActiveSet);

    // RB: This gives the correct gradient for the length scale
    mat partialDeriv(sizeActiveSet, sizeActiveSet);
    mat U = backslash(KB_new,KB);

    W +=  KB * (C + outer_product(alpha, alpha));

    for(int i = 0; i < covFunc.getNumberParameters(); i++)
    {
        covFunc.getParameterPartialDerivative(partialDeriv, i, ActiveSet);
        mat V1 = backslash(KB_new,partialDeriv);
        mat V2 = W*backslash(KB_new,partialDeriv*U);

        grads(i) = trace(V1-V2);
    }

    return grads;
}

/**
 * Set the likelihood type
 */
void PSGP::setLikelihoodType(LikelihoodCalculation lc)
{
    likelihoodType = lc;
}


/**
 * Display current parameters
 */
void PSGP::displayModelParameters() const
{
    cout << "Summary Sequential Gaussian Process" << endl;
    cout << "  Kernel Matrix size         : " << KB.rows() << " x " << KB.cols() << endl;
    cout << "  Inverse Kernel Matrix size : " << Q.rows() << " x " << Q.cols() << endl;
    cout << "  alpha size                 : " << alpha.size() << endl;
    cout << "  C size                     : " << C.rows() << " x " << C.cols() << endl;
    cout << "  Projection matrix size     : " << P.rows() << " x " << P.cols() << endl;
    cout << "  Lambda                     : " << varEP.size() << endl;
    cout << "  projection alpha           : " << meanEP.size() << endl;
    cout << "  log evidence vector        : " << logZ.size() << endl;
    cout << "  ----------------------------" << endl;
    cout << "  Predicion locations        : " << Locations.rows() << " x " << Locations.cols() << endl;
    cout << "  Observations               : " << Observations.size() << endl;
    cout << "  Active set size            : " << ActiveSet.rows() << " (max = " << maxActiveSet << ")" << endl;
    cout << "  Epsilon tolerance          : " << epsilonTolerance << endl;
    cout << "  Iterations Changing/Fixed  : " << iterChanging << "/" << iterFixed << endl;
    cout << "  Moment projection          : " << momentProjection << endl;
}


/**
 * Compute cholesky decomposition of matrix
 * 
 * RB: This belongs in a different class/library
 * TODO: Move to ITPPExt
 */
mat PSGP::computeCholesky(const mat& iM) const 
{
    mat M = iM;
    assert(M.rows() == M.cols());

    const double ampl = 1.0e-10;
    const int maxAttempts = 10;

    mat cholFactor(M.rows(), M.cols());

    int l = 0;
    bool success = chol(M, cholFactor);
    if(success)
    {
        return cholFactor;
    }
    else
    {
        double noiseFactor = abs(ampl * (trace(M) / double(M.rows())));
        while(!success)
        {
            M = M + (noiseFactor * eye(M.rows()));

            if(l > maxAttempts)
            {
                cerr << "Unable to compute cholesky decomposition" << endl;
                break;
            }
            l++;
            noiseFactor = noiseFactor * 10;
            success = chol(M, cholFactor);
        }
        cout << "Matrix not positive definite.  After " << l << " attempts, " << noiseFactor << " added to the diagonal" << endl;
    }
    return cholFactor;

}

/**
 * Compute inverse of a square matrix using Cholesky decomposition
 * 
 * TODO: Move to ITPPExt
 */
mat PSGP::computeInverseFromCholesky(const mat& C) const
{
    mat cholFactor = computeCholesky(C);
    mat invChol = backslash(cholFactor, eye(cholFactor.rows()));
    return invChol * invChol.transpose();
}

