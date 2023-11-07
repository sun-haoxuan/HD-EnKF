RK4 = function(x, t, h, p, Fc){
  L96 = function(x, t){
    index1 = c(2:p, 1)
    index2 = c(p - 1, p, 1:(p - 2))
    index3 = c(p, 1:(p - 1))
    d = (x[index1] - x[index2]) * x[index3] - x[1:p] + Fc
    return(d)
  }
  k1 = h * L96(x, t)
  k2 = h * L96(x + 0.5 * k1, t + 0.5 * h)
  k3 = h * L96(x + 0.5 * k2, t + 0.5 * h)
  k4 = h * L96(x + k3, t + h)
  return(x + (k1 + 2 * k2 + 2 * k3 + k4) / 6)
}