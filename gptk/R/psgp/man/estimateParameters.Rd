\name{estimateParameters.psgp}
\alias{estimateParameters.psgp}
\alias{estimateParameters}
\title{Parameter estimation using a Projected Sequential Gaussian Process (PSGP)}
\description{
  \code{estimateParameters} overloads the \code{\link[=intamap]{estimateParameters}}
  routine from the Intamap package for interpolation using the PSGP method.
}
\synopsis{ 
  estimateParameters(object) 
}
\usage{
  estimateParameters(object)
}
\arguments{
  \item{object}{ a list object of Intamap type. Most arguments necessary for
    interpolation are passed through this object. 
    See \link[intamap]{00Introduction} for further description of the necessary 
    content of this variable.
  }
}

\details{
  See \code{\link{psgp-package}} and \code{\link{learnParameters}} for 
  further detail.
}

\author{Remi Barillec}
\seealso{
  \code{\link{learnParameters}},
  \code{\link[intamap]{estimateParameters}},
  \code{\link{makePrediction}},
  \code{\link[intamap]{createIntamapObject}}
}
\examples{
# load our favourite dataset
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
  class = "psgp"    # Use PSGP for parameter estimation/interpolation
)

# do interpolation step:
obj = conformProjections(obj)
obj = estimateParameters(obj)
}
\keyword{spatial}


