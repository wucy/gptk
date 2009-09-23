\name{learnParameters}
\alias{learnParameters}
\title{Projected Sequential Gaussian Process}
\description{\code{learnParameters} performs maximum likelihood parameter estimation in the PSGP framework.
}
\synopsis{ learnParameters(object) }
\usage{
learnParameters(object)
}
\arguments{
  \item{object}{ a list object of Intamap type. Most arguments necessary for interpolation
  are passed through this object. See \link[intamap]{00Introduction} for further 
  description of the necessary content of this variable}
}
 

\details{
The function \code{learnParameters} is a function for estimating variogram parameters with 
Projected Spatial Gaussian Process (PSGP) methods (Csato and Opper, 2002; Ingram et al., 2008)
through a maximum likelihood estimation. Predictions can be done with 
\code{\link{makePrediction}} or \code{\link[intamap]{spatialPredict}} with
an \code{object} of type \code{psgp}.
These methods are able to also take the measurement characteristics into account,
in this function implemented as the element \code{obsChar} in \code{object}.

Most of the method is implemented in C++, relying on the external library IT++
(\url{http://itpp.sourceforge.net}), which is a C++ library composed of
classes and functions for linear algebra (matrices and vectors). 

Instead of calling this function directly, a user is advised to call the wrapper 
function \code{\link[intamap]{estimateParameters}} with an \code{object} of class \code{psgp}.
}


\references{ 

\url{http://www.intamap.org/}

}
\author{Ben Ingram}
\seealso{
\code{\link[intamap]{estimateParameters}},\code{\link[intamap]{makePrediction}}
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
