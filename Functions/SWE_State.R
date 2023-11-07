SWE_State = function(opt){
  
  ## parameter ####
  opt$k = opt$k.true
  
  ## Generate true state
  nx = opt$nx # grid size of x
  ny = opt$ny # grid size of y
  L = opt$L
  D = opt$D
  g = opt$g # gravitational constant
  f = opt$f
  k = opt$k
  dt = opt$dt # hardwired timestep
  h0 = opt$h0
  h1 = opt$h1
  h2 = opt$h2
  
  ## Generate observations
  id = opt$id # rows with observation
  sigma.obs = opt$sigma.obs # observation error
  rho = opt$rho ## observation
  vf = opt$vf # for Wrong R
  
  ## Assimilation Step
  S = opt$S
  S.o.start = opt$S.o.start ## Start of the observation
  S.o.end = opt$S.o.end # assimilate 48h, then forecast 24h
  S.o.seq = opt$S.o.seq # assimilate every 3h
  
  ngrid <- nx * ny
  dx <- L / nx
  dy <- D / (ny - 1)
  nobs <- (S.o.end - S.o.start) %/% S.o.seq + 1 # forecast 1h, then assimilate
  ndim <- ngrid * 3
  nid <- length(id)
  xobs <- ny * nid * 3 # observe h only
  
  ## Generate true state ####
  x <- array(rep(0:(nx - 1), each = ny), dim = c(ny, nx)) / nx * L
  y <- array(rep(0:(ny - 1), times = nx), dim = c(ny, nx)) / (ny - 1) * D
  
  H0 <- matrix(1, nrow = ny + 2, ncol = nx + 2)
  U0 <- matrix(0, nrow = ny + 2, ncol = nx + 2)
  V0 <- matrix(0, nrow = ny + 2, ncol = nx + 2)
  Hx <- matrix(0, nrow = ny + 2, ncol = nx + 1)
  Ux <- matrix(0, nrow = ny + 2, ncol = nx + 1)
  Vx <- matrix(0, nrow = ny + 2, ncol = nx + 1)
  Hy <- matrix(0, nrow = ny + 1, ncol = nx + 2)
  Uy <- matrix(0, nrow = ny + 1, ncol = nx + 2)
  Vy <- matrix(0, nrow = ny + 1, ncol = nx + 2)
  
  H0[2:(ny + 1), 2:(nx + 1)] <- h0 + h1 * tanh(9 * (D / 2 - y) / (2 * D)) + h2 / cosh(9 * (D / 2 - y) / D)^2 * sin(2 * pi / L * x)
  U0[2:(ny + 1), 2:(nx + 1)] <- -g / f * (-h1 * (1 - tanh(9 * (D / 2 - y) / (2 * D))^2) * 9 / (2 * D) + 2 * h2 / cosh(9 * (D / 2 - y) / D)^2 * tanh(9 * (D / 2 - y) / D) * sin(2 * pi / L * x) * 9 / D)
  V0[2:(ny + 1), 2:(nx + 1)] <- g / f * h2 / cosh(9 * (D / 2 - y) / D)^2 * cos(2 * pi / L * x) * 2 * pi / L
  
  x.t <- matrix(NA, nrow = ndim, ncol = S)  # true state
  x.t[, 1] <- c(as.vector(U0[2:(ny + 1), 2:(nx + 1)]), as.vector(V0[2:(ny + 1), 2:(nx + 1)]), as.vector(H0[2:(ny + 1), 2:(nx + 1)]))
  for (t in 2:S) {
    x.t[, t] <- SWE_Forecast(x.t[, t - 1], opt)
  }
  
  ## Observation Matrix ####
  H <- matrix(0, nrow = xobs, ncol = ndim)
  for (i in 1:nid) {
    jy <- (id[i] - 1) * ny + 1:ny
    jx <- ((i - 1) * ny + 1):(i * ny)
    for (ii in 1:ny) {
      H[jx[ii], jy[ii]] <- 1
      H[jx[ii] + ny * nid, jy[ii] + ngrid] <- 1
      H[jx[ii] + ny * nid * 2, jy[ii] + ngrid * 2] <- 1
    }
  }
  
  ## Observation error covariance matrix ####
  aa1 <- matrix(NA, nrow = 2, ncol = ny * nid)
  for (i in 1:nid) {
    iy <- ((i - 1) * ny + 1):(i * ny)
    for (ii in 1:ny) {
      aa1[1, iy[ii]] <- ii
      aa1[2, iy[ii]] <- id[i]
    }
  }
  aa2 <- diag(1, ny * nid)
  for (i in 1:(ny * nid - 1)) {
    for (j in (i + 1):(ny * nid)) {
      dd <- sqrt((aa1[1, i] - aa1[1, j])^2 + (aa1[2, i] - aa1[2, j])^2)
      aa2[i, j] <- rho^dd
      aa2[j, i] <- aa2[i, j]
    }
  }
  R <- matrix(0, nrow = xobs, ncol = xobs)
  R[1:(ny * nid), 1:(ny * nid)] <- aa2 * vf * sigma.obs[1]^2
  R[(ny * nid + 1):(ny * nid * 2), (ny * nid + 1):(ny * nid * 2)] <- aa2 * vf * sigma.obs[2]^2
  R[(ny * nid * 2 + 1):(ny * nid * 3), (ny * nid * 2 + 1):(ny * nid * 3)] <- aa2 * vf * sigma.obs[3]^2
  
  ## Observation State ####
  t.o <- rep(0, times = S)
  for (i in seq(S.o.start, S.o.end, by = S.o.seq)){
    t.o[i] <- 1
  }
  so <- H %*% x.t[, t.o == 1]
  ww <- matrix(NA, nrow = ny * nid, ncol = nobs)
  for(i in 1:nobs){
    ww[, i] <- rnorm(ny * nid, 0, 1)
  }
  bb1 <- svd(R[1:(ny * nid), 1:(ny * nid)] / vf)
  sqrR_true1 <- bb1$u %*% diag(sqrt(bb1$d), ny * nid) # generate observation (true R)
  bb2 <- svd(R[(ny * nid + 1):(ny * nid * 2), (ny * nid + 1):(ny * nid * 2)] / vf)
  sqrR_true2 <- bb2$u %*% diag(sqrt(bb2$d), ny * nid) # generate observation (true R)
  bb3 <- svd(R[(ny * nid * 2 + 1):(ny * nid * 3), (ny * nid * 2 + 1):(ny * nid * 3)] / vf)
  sqrR_true3 <- bb3$u %*% diag(sqrt(bb3$d), ny * nid) # generate observation (true R)
  obs_err1 <- sqrR_true1 %*% ww
  obs_err2 <- sqrR_true2 %*% ww
  obs_err3 <- sqrR_true3 %*% ww
  y.o <- matrix(NA, nrow = xobs, ncol = nobs)
  y.o[1:(ny * nid), ] <- so[1:(ny * nid), ] + obs_err1
  y.o[(ny * nid + 1):(ny * nid * 2), ] <- so[(ny * nid + 1):(ny * nid * 2), ] + obs_err2
  y.o[(ny * nid * 2 + 1):(ny * nid * 3), ] <- so[(ny * nid * 2 + 1):(ny * nid * 3), ] + obs_err3
  
  ## Output ####
  output = list(
    true = x.t,
    observation = y.o,
    option = opt
  )
  
  return(output)
}
