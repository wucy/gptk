\name{learnParameters}
\alias{learnParameters}
\title{Projected Sequential Gaussian Process}
\description{
  \code{learnParameters} performs maximum likelihood parameter estimation in the PSGP framework.
}
\synopsis{ 
  learnParameters(object) 
}
\usage{
  learnParameters(object)
}
\arguments{
  \item{object}{ a list object of Intamap type. Most arguments necessary for
    interpolation are passed through this object. 
    See \link[intamap]{00Introduction} for further description of the necessary 
    content of this variable.
  }
}

\section{Warning}{It is advised to use the Intamap wrapper \code{\link[=intamap]{estimateParameters}} rather than calling this method directly.}

\details{
  The Projected Spatial Gaussian Process (PSGP) framework provides  
  an approximation to the full Gaussian process in which the observations 
  are projected sequentially onto an optimal subset of 'active' observations. Spatial 
  interpolation is done using a mixture of covariance kernels (exponential and Matern 
  5/2).

  The function \code{learnParameters} is an internal function for estimating the
  parameters of the covariance function given the data, using a maximum likelihood
  approach. A valid Intamap \code{object} must be passed in.
  
  PSGP is able to also take the measurement characteristics (i.e. errors) into
  account using possibly many error models. For each error model, assumed Gaussian, a 
  variance (and optionally, a bias) can be specified. The vector 
  \code{object$observations$oevar} contains the variances for the error models (one 
  value per error model). Similarly, if specified, the vector  
  \code{object$observations$oebias} must contains the bias values for each model. Which 
  error model is used for a given observation is determined by the 
  \code{object$observations$oeid} vector of indices, which specifies the index of the 
  model to be used for each observation. In the current release, the biases are effectively ignored (the following release will aim to support them).
  
  Predictions can be done with \code{\link{makePrediction}} although it is advised to 
  use the intamap wrapper function \code{\link[=intamap]{spatialPredict}} instead, with 
  an \code{object} of type \code{psgp}.
}


\references{ 
    Csato and Opper, 2002; Ingram et al., 2008
    \url{http://www.intamap.org/}
}
\author{Ben Ingram, Remi Barillec}
\seealso{
  \code{\link{makePrediction}},
  \code{\link{learnParameters}},
  \code{\link[intamap]{estimateParameters}},
  \code{\link[intamap]{createIntamapObject}}
}
\examples{
# This example skips some steps that might be necessary for more complicated
# tasks, such as estimateParameters and pre- and postProcessing of the data
data(meuse)
coordinates(meuse) = ~x+y
meuse$value = log(meuse$zinc)
data(meuse.grid)
gridded(meuse.grid) = ~x+y
proj4string(meuse) = CRS("+init=epsg:28992")
proj4string(meuse.grid) = CRS("+init=epsg:28992")

# set up intamap object:
obj = createIntamapObject(
	observations = meuse,
	predictionLocations = meuse.grid,
	targetCRS = "+init=epsg:3035",
	class = "psgp"
)

# do interpolation step:
obj = conformProjections(obj)
obj = estimateParameters(obj)  # obj updated with variogram
}
\keyword{spatial}


