learnParameters <- function(object) {

save(object, file="Output.rdata")

#since we are using the full covariance matrix here, we need
#to limit the number of observations
inputs = object$observations
pred = object$predictionLocations

# JOS: the line below used to be higher and tried to access nrow(m), which
#      did not exist
#      Is it not rather arbitrary to choose the 1000 first observations?
#      Would it make more sense to sample from the observations if you want
#      to limit the number?
#maxObs = min(nrow(inputs), 1000)
#inputs = inputs[1:maxObs,]

# put data into an easy parseable format for the backend C++ code
x = coordinates(inputs)

# JOS: Allow for other names for the dependent variable
depVar = as.character(object$formulaString[[2]])
y = as.numeric(inputs@data[[depVar]])

# put obsError and sensorID into vectors
observationError = as.integer(inputs$oeid)
sensorModel = as.integer(inputs$sensor)

# vector of strings of meta data
metaData = object$obsChar

# error variance vector
# JOS: Introduced a check to see if var is exisiting in the object
#      Giving 0 otherwise
#      Is e supposed to include indices to error models, or only variances?      
#      Would you need to be able to pass uncObject as well, as we discussed 
#      in Wageningen
if ("var" %in% names(inputs)) {
  e = as.numeric(inputs@data$var)
} else {
  e = rep(0,length(y))
}


vario = array()
# some defaults
#vario[1] = 2; # variogram model
#vario[2] = 1; # range
#vario[3] = 1; # sill
#vario[4] = 0.1; # nugget

# JOS: Would it not be better to use automap here? 
#      Your default seems rather arbitrary...
#      The function only searches among the models you have defined in 
#      vmodels 
# vmodels = c("Exp","Gau")
vmodels = c("Exp");
varioModel = autofitVariogram(object$formulaString,
                  object$observations,model=vmodels)$var_model
# vario[1] = varioModel$model[2]
vario[1] = 2    # Exponential model
vario[2] = varioModel[[3]][2]
vario[3] = varioModel[[2]][2]
vario[4] = varioModel[[2]][1]

# call the C code
r <- .Call("estParam", x, y, e, vario, observationError, sensorModel, metaData,
	PACKAGE = "psgp")

# there must be a better way of doing this?
# does R have switch/case?
# JOS: switch exists (?switch), but I think this might be even better:

# PSGP covariance parameters (log transformed)
# These are stored in r from 6th element onwards
# We store them in the variogram model - this is a hack really, we just
# use the variogram model data frame for storage!

# Create dummy variogram model - this allows us to store 2x9 values
# as the first column is for text
vmodel = vgm(0, "Exp", 1, 0)
# Use N/A to indicate that something is wrong with this variogram model
# This should make sure it is not used or returns an error if it is.
vmodel[1:2,1] = NA    
vmodel[1,2:9] = r[1:8]
vmodel[2,2:9] = r[9:16]

object$variogramModel = vmodel
object
}

