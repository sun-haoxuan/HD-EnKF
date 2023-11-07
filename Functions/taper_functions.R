banding = function(z){
  if(z <= 1){
    return(1)
  }else{
    return(0)
  }
}

GC = function(z){
  z = z * 2
  if(z >= 0 & z <= 1){
    return(1 - 5 * z ^ 2 / 3 + 5 * z ^ 3 / 8 + z ^ 4 / 2 - z ^ 5 / 4)
  }else if(z > 1 & z <= 2){
    return(- 2 / 3 / z + 4 - 5 * z + 5 * z ^ 2 / 3 + 5 * z ^ 3 / 8 - z ^ 4 / 2 + z ^ 5 / 12)
  }else{
    return(0)
  }
}