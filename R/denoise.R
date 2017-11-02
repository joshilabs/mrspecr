gs1 <- function(u0, h = 1, beta = 5, lambda = 0.08, iter=20) {

    f <- u0
    u1 <- u0
    n <- length(u0) 
    u <- numeric(n+2)
    u[2:(n+1)] <- u0
    
    for (it in 1:iter) {
      u[1] <- u[2]
      u[n+2] <- u[n+1]
      for (jj in 2:(n+1)) {
        nu <- ( u[jj+1] - u[jj-1] )/(2*h)
        nu <- beta + nu*nu
        nu <- nu*sqrt(nu)
        nu <- 1/nu
        alpha <- beta*nu
        den <- (2*alpha)/(lambda*h*h)
        a <- 1/(1 + 2*den)
        b <- den/(1 + 2*den)
        u[jj] <- a*f[jj-1] + b*( u[jj-1] + u[jj+1])
      }
    }
    u1 <- u[2:(n+1)]
}

