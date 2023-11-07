SWE_Forecast <- function(s_l, opt){
  
  nx = opt$nx
  ny = opt$ny
  L = opt$L
  D = opt$D
  dt = opt$dt
  f = opt$f
  g = opt$g
  k.set = opt$k
  
  n = nx * ny
  dx <- L / nx
  dy <- D / (ny - 1)
  
  H0 <- matrix(1, nrow = ny + 2, ncol = nx + 2)
  U0 <- matrix(0, nrow = ny + 2, ncol = nx + 2)
  V0 <- matrix(0, nrow = ny + 2, ncol = nx + 2)
  Hx <- matrix(0, nrow = ny + 2, ncol = nx + 1)
  Ux <- matrix(0, nrow = ny + 2, ncol = nx + 1)
  Vx <- matrix(0, nrow = ny + 2, ncol = nx + 1)
  Hy <- matrix(0, nrow = ny + 1, ncol = nx + 2)
  Uy <- matrix(0, nrow = ny + 1, ncol = nx + 2)
  Vy <- matrix(0, nrow = ny + 1, ncol = nx + 2)
  
  # Reflective boundary conditions
  H0[2:(ny + 1), 2:(nx + 1)] <- array(s_l[(2 * n + 1):(3 * n)], dim = c(ny, nx))
  U0[2:(ny + 1), 2:(nx + 1)] <- array(s_l[1:n], dim = c(ny, nx)) * H0[2:(ny + 1), 2:(nx + 1)]
  V0[2:(ny + 1), 2:(nx + 1)] <- array(s_l[(n + 1):(2 * n)], dim = c(ny, nx)) * H0[2:(ny + 1), 2:(nx + 1)]
  
  H0[, 1] <- H0[, nx + 1]
  U0[, 1] <- U0[, nx + 1]
  V0[, 1] <- V0[, nx + 1]
  H0[, nx + 2] <- H0[, 2]
  U0[, nx + 2] <- U0[, 2]
  V0[, nx + 2] <- V0[, 2]
  H0[1, ] <- H0[2, ]
  U0[1, ] <- U0[2, ]
  V0[1, ] <- -V0[2, ]
  H0[ny + 2, ] <- H0[ny + 1, ]
  U0[ny + 2, ] <- U0[ny + 1, ]
  V0[ny + 2, ] <- -V0[ny + 1, ]
  
  # First half step
  # x direction
  i <- 1:(ny + 2)
  j <- 1:(nx + 1)
  # height
  Hx[i, j] <- (H0[i, j] + H0[i, j + 1]) / 2 - dt / (2 * dx) * (U0[i, j + 1] - U0[i, j])
  # x momentum
  Ux[i, j] <- (U0[i, j] + U0[i, j + 1]) / 2 - dt / (2 * dx) * ((U0[i, j + 1]^2 / H0[i, j + 1] + g / 2 * H0[i, j + 1]^2) - (U0[i, j]^2 / H0[i, j] + g / 2 * H0[i, j]^2))
  # y momentum
  Vx[i, j] <- (V0[i, j] + V0[i, j + 1]) / 2 - dt / (2 * dx) * ((U0[i, j + 1] * V0[i, j + 1] / H0[i, j + 1]) - (U0[i, j] * V0[i, j] / H0[i, j]))
  
  # y direction
  i <- 1:(ny + 1)
  j <- 1:(nx + 2)
  # height
  Hy[i, j] <- (H0[i, j] + H0[i + 1, j]) / 2 - dt / (2 * dy) * (V0[i + 1, j] - V0[i, j])
  # x momentum
  Uy[i, j] <- (U0[i, j] + U0[i + 1, j]) / 2 - dt / (2 * dy) * ((V0[i + 1, j] * U0[i + 1, j] / H0[i + 1, j]) - (V0[i, j] * U0[i, j] / H0[i, j]))
  # y momentum
  Vy[i, j] <- (V0[i, j] + V0[i + 1, j]) / 2 - dt / (2 * dy) * ((V0[i + 1, j]^2 / H0[i + 1, j] + g / 2 * H0[i + 1, j]^2) - (V0[i, j]^2 / H0[i, j] + g / 2 * H0[i, j]^2))
  
  # Second half step
  i <- 2:(ny + 1)
  j <- 2:(nx + 1)
  # height
  H0[i, j] <- H0[i, j] - (dt / dx) * (Ux[i, j] - Ux[i, j - 1]) - (dt / dy) * (Vy[i, j] - Vy[i - 1, j])
  # x momentum
  U0[i, j] <- U0[i, j] + dt * f * V0[i, j] - (dt / dx) * ((Ux[i, j]^2 / Hx[i, j] + g / 2 * Hx[i, j]^2) - (Ux[i, j - 1]^2 / Hx[i, j - 1] + g / 2 * Hx[i, j - 1]^2)) - (dt / dy) * ((Vy[i, j] * Uy[i, j] / Hy[i, j]) - (Vy[i - 1, j] * Uy[i - 1, j] / Hy[i - 1, j]))
  # y momentum
  V0[i, j] <- V0[i, j] - dt * f * U0[i, j] - (dt / dx) * ((Ux[i, j] * Vx[i, j] / Hx[i, j]) - (Ux[i, j - 1] * Vx[i, j - 1] / Hx[i, j - 1])) - (dt / dy) * ((Vy[i, j]^2 / Hy[i, j] + g / 2 * Hy[i, j]^2) - (Vy[i - 1, j]^2 / Hy[i - 1, j] + g / 2 * Hy[i - 1, j]^2))
  
  H0[, 1] <- H0[, nx + 1]
  U0[, 1] <- U0[, nx + 1]
  V0[, 1] <- V0[, nx + 1]
  H0[, nx + 2] <- H0[, 2]
  U0[, nx + 2] <- U0[, 2]
  V0[, nx + 2] <- V0[, 2]
  H0[1, ] <- H0[2, ]
  U0[1, ] <- U0[2, ]
  V0[1, ] <- -V0[2, ]
  H0[ny + 2, ] <- H0[ny + 1, ]
  U0[ny + 2, ] <- U0[ny + 1, ]
  V0[ny + 2, ] <- -V0[ny + 1, ]
  
  U0 <- U0 / H0
  V0 <- V0 / H0
  i <- 2:(ny + 1)
  j <- 2:(nx + 1)
  U0[i, j] <- U0[i, j] + dt / dx^2 * k.set * (U0[i, j + 1] + U0[i, j - 1] - 2 * U0[i, j]) + dt / dy^2 * k.set * (U0[i + 1, j] + U0[i - 1, j] - 2 * U0[i, j])
  V0[i, j] <- V0[i, j] + dt / dx^2 * k.set * (V0[i, j + 1] + V0[i, j - 1] - 2 * V0[i, j]) + dt / dy^2 * k.set * (V0[i + 1, j] + V0[i - 1, j] - 2 * V0[i, j])
  H0[i, j] <- H0[i, j] + dt / dx^2 * k.set * (H0[i, j + 1] + H0[i, j - 1] - 2 * H0[i, j]) + dt / dy^2 * k.set * (H0[i + 1, j] + H0[i - 1, j] - 2 * H0[i, j])
  
  output = c(as.vector(U0[2:(ny + 1), 2:(nx + 1)]), as.vector(V0[2:(ny + 1), 2:(nx + 1)]), as.vector(H0[2:(ny + 1), 2:(nx + 1)]))
  
  return(output)
}