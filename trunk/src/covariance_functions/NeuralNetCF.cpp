#include "NeuralNetCF.h"

/*
 * Constructor - pass in the length scale, variance and offset
 */
NeuralNetCF::NeuralNetCF(double ls, double var, double offst)
: CovarianceFunction("Neural network covariance function", 3), lengthScale(parameters[0]),
variance(parameters[1]), offset(parameters[2])
{
    // Make sure parameters are positive 
    assert( ls > 0.0 && var > 0.0 && offset > 0.0);  
    
    lengthScale = ls;   // Length scale
    variance = var;
    offset = offst;
    
    // Parameter names
    parametersNames[0] = "length scale";
    parametersNames[1] = "variance";
    parametersNames[2] = "offset length scale";
}


/**
 * Destructor
 */ 
NeuralNetCF::~NeuralNetCF()
{
}

/**
 * Covariance between two points A and B
 */
inline double NeuralNetCF::computeElement(const vec& A, const vec& B) const
{
    vec Ahat = A/lengthScale;
    vec Bhat = B/lengthScale;
        
    double gamma = 1.0 / sqr(offset);
    
    return variance * asin( (gamma+dot(Ahat,Bhat)) / 
                            sqrt( (1.0+gamma+dot(Ahat,Ahat)) * (1.0+gamma+dot(Bhat,Bhat)) )
                          );
}

/** 
 * Gradient of cov(X) w.r.t. given parameter number
 */
void NeuralNetCF::covarianceGradient(mat& grad, const int parameterNumber, const mat& X) const
{
    assert(parameterNumber < getNumberParameters());
    assert(parameterNumber >= 0);

    Transform* t = getTransform(parameterNumber);
    double gradientModifier = t->gradientTransform(getParameter(parameterNumber));

    switch(parameterNumber)
    {
        case 0 :
        {
            double gamma = 1.0 / sqr(offset);
            mat Xhat = X / lengthScale;
            
            vec Xi, Xj;
            for (int i=0; i<Xhat.rows(); i++)
            {
                Xi = Xhat.get_row(i);
                double vi = 1.0 + gamma + dot(Xi, Xi);
                
                for (int j=0; j<i; j++) 
                {
                    Xj = Xhat.get_row(j);
                    
                    double u  = gamma + dot(Xi,Xj);
                    double vj = 1.0 + gamma + dot(Xj, Xj);
                    // double vij = vi*vj;
                    // PD(i,j) = u * (vi+vj) / (vij * sqrt(vij - u*u) );
                    
                    double vij = sqrt(vi*vj);
                    grad(i,j) = (2*gamma*vij - (gamma+1)*u*(vi+vj)/vij) / (vij*sqrt(vij*vij-u*u));
                    grad(j,i) = grad(i,j);
                }
                // PD(i,i) = 2.0*(vi-1.0) / (vi*sqrt(2*vi-1.0));
                grad(i,i) = 2.0*(gamma+1.0-vi) / (vi*sqrt(2*vi-1.0));
            }
            grad *= variance * gradientModifier / lengthScale;
            break;
        }

        case 1 :
        {
            covariance(grad, X);
            grad *= (gradientModifier / variance);
            break;
        }

        case 2:
            double gamma = 1.0 / sqr(offset);
            mat Xhat = X / lengthScale;

            vec Xi, Xj;
            for (int i=0; i<Xhat.rows(); i++) 
            {
                Xi = Xhat.get_row(i);
                double vi = 1.0 + gamma + dot(Xi, Xi);

                for (int j=0; j<i; j++) 
                {
                    Xj = Xhat.get_row(j);

                    double u  = gamma + dot(Xi,Xj);
                    double vj = 1.0 + gamma + dot(Xj, Xj);
                    double vij = vi*vj;

                    grad(i,j) = (1.0 - 0.5*u*(vi+vj)/vij) / sqrt(vij-u*u);
                    grad(j,i) = grad(i,j);
                }
                grad(i,i) = 1.0 / ( vi * sqrt(2*vi-1.0) );
            }
            grad *= - 2.0 * variance * gradientModifier / pow(offset,3.0);
            break;
    }

}
