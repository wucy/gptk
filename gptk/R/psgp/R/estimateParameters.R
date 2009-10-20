estimateParameters.psgp = function(object,...) {

  origObs = object$observations
  
  rotated = FALSE
  if (object$params$doAnisotropy) {
      object = estimateAnisotropy(object) 
      #rotate Data
    if (object$anisPar$doRotation && all(as.character(object$formulaString[[3]])=="1"))
      object$observations=rotateAnisotropicData(object$observations,object$anisPar)
    rotated = TRUE
  }  

  #if (require(astonGeostats)) 
  object = learnParameters(object)

  if (rotated) 
    object$observations = origObs
  object
}