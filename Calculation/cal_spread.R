cal_spread = function(x.en.temp){
  return(sqrt(mean((x.en.temp - apply(x.en.temp, 1, mean)) ^ 2)))
}