get_DAmethod = function(method){
  if(method %in% c('HD-EnKF', 'localization')){
    DAinfo = list(
      taper = GC,
      interval = c((log(p) / n) ^ (-1 / 2) / 5, (log(p) / n) ^ (-1 / 2) * 5),
      num_eigen = p,
      distMat = distMat
    )
  }else{
    DAinfo = NULL
  }
  
  return(DAinfo)
}