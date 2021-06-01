#' getAtomicWeight
#'
#' This function return the standard atomic weight and its uncertainty for a given element.
#' The standard atomic weights are calculated using the reported interval values for  normal and naturally  occurring materials
#' according to the 2019 Table of Standard Atomic Weights published at \url{https://ciaaw.org/atomic-weights.htm} (1), based on
#' IUPAC Technical Report: Atomic weights of the elements 2013(2). The Commission on Isotopic Abundances and Atomic Weights (CIAAW)
#' defines a normal material as any terrestrial material
#' that "is a reasonably possible source for this element or its compounds in commerce, for industry or science;
#' the material is not itself studied for some extraordinary anomaly and its isotopic composition has not been
#' modified significantly in a geologically brief period"(3).
#'
#'  This function takes the more simple aproximation for standard uncertaninty calculation, based on the
#'  reported interval with the asumption of a rectangular distribution. Decisional uncertainty is returned when is reported(3).
#'  More detailed uncertainty estimations considering isotopic abundances for an specific application requires more complex models such as the reported in (2-4)
#'
#' @param element A chemical symbol for the element
#' @return A data.frame with the standard atomic weight (AW) and its standard uncertainty (u_AW)
#' @examples
#' getAtomicWeight("H")
#' @seealso [getMolarMass()]
#' @references
#' (1) CIAAW, Standard atomic weights [online], \url{https://ciaaw.org/atomic-weights.htm}. Retrieved: 2021-05-05
#'
#' (2) Juris Meija, et. al., Atomic weights of the elements 2013 (IUPAC Technical Report), Pure Appl. Chem. 2016; 88(3): 265–291. \url{https://doi.org/10.1515/pac-2015-0305}
#'
#' (3) Possolo, A.,  van der Veen, A.,   Meija, J.,  Hibbert, D.B., Interpreting and propagating the uncertainty of the standard atomic weights,  Pure Appl. Chem. 2018; 90(2): 395–424.\url{https://doi.org/10.1515/pac-2016-0402}
#'
#' (4) Juris Meija and Zoltan Mester. Atomic weight uncertainty calculation from isotopic composition of the elements, Metrologia 45 (2008) 459–463.\url{https://doi.org/10.1088/0026-1394/45/4/012}
#'
#' @export
getAtomicWeight<-function(element){
  atom<-atomicWeightsData[atomicWeightsData$Atom==element,]
  if( length(atom$Atom)>0){
    if( is.na(atom$u)){
      data.frame(AW=mean(c(atom$Lower_weight,atom$Upper_weight)),u_AW=(atom$Upper_weight-atom$Lower_weight)/(2*sqrt(3)))
    }else{
      data.frame(AW=atom$Lower_weight,u_AW=atom$u)
    }
  }
}

#' getMolarMass
#'
#' This function return the standard  molar mass and its standard uncertainty for a given molecular formula.
#' The standard atomic weights used are obtained from [getAtomicWeight()]
#'
#'  The function uses the GUM aproach based on the molecular formula(1).
#'
#' @param atomsType A string vector with chemical symbols for the elements in the molecular formula
#' @param atomsNumber A int vector with the number of atoms for each element in the molecular formula (by default =1)
#' @return A data.frame with the molar mass (MM) and its standard uncertainty (u_MM)
#' @examples
#' getMolarMass(c("C","H"),c(1,4))
#' @seealso [getAtomicWeights()]
#' @references
#' (1) Possolo, A.,  van der Veen, A.,   Meija, J.,  Hibbert, D.B., Interpreting and propagating the uncertainty of the standard atomic weights,  Pure Appl. Chem. 2018; 90(2): 395–424.\url{https://doi.org/10.1515/pac-2016-0402}
#' @export
getMolarMass<-function(atomsType,atomsNumber=1){
  atomsNumber_=atomsNumber
  if(length(atomsNumber)==1){atomsNumber_=rep(atomsNumber,length(atomsType))}
  aws=sapply(atomsType,getAtomicWeight)
  MM=sum(unlist(aws[1,])*atomsNumber_)
  u_MM=sqrt(sum((unlist(aws[2,])*atomsNumber_)^2))
  data.frame(MM=MM,u_MM=u_MM)
}


