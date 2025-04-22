# SIR Gillespie

# First we define our one-step function.
SIR.onestep <- function (x, params) { #function to calculate one step of stochastic SIR
  
  t  <- x[1] 
  S  <- x[2] 
  I  <- x[3]
  R  <- x[4]
  
  beta = params[1]
  gamma = params[2]
  
  N = S + I + R 
  
  rates <- c(beta*S*I/N, gamma*I)
  
  changes <- matrix(
    c(-1,1,0,  # beta*S*I/N
      0,-1,1   # gamma*I
    ),
    ncol=3, byrow=TRUE)
  
  tau <- rexp(n=1,rate=sum(rates)) # exponential waiting time
  U <- runif(1) #uniform random deviate 
  m <- min(which(cumsum(rates)>=U*sum(rates)))
  x <- x[2:4] + changes[m,]
  t <- t + tau
  return(out <- c(t, x))
}


# Next we write our function for simulating a whole epidemic.
SIR_Gillespie <- function (xstart, params, tfinal) { #function to simulate stochastic SIR
  output <- array(dim=c(1,4)) # set up array to store results
  colnames(output) <- c("time","S","I","R") #name variables
  output[1,] <- xstart # first record of output is initial condition
  k=1
  while (as.logical(output[k,1]<tfinal)){
    out.t <- SIR.onestep(output[k,],params)
    output <- rbind(output,out.t)
    checkEnd <- out.t[3] <= 0 
    if(checkEnd){break}
    k=k+1
    #cat("AFS Sim: ",checkEnd, "\n", sep=""); flush.console()
  }
  rval <- as.data.frame(output, row.names = F)
  return(rval) #return output
}

# # Example
# initial <- c(time=1,S=999,I=1,R=0) #initial conditions
# par <- c(beta=0.3, gamma=0.1)
# 
# # Run model
# out <- as.data.frame(
#   SIR_Gillespie(
#     xstart = initial,
#     params = par, 
#     tfinal = 100
#   )
# )
# 
# plot(out$time, out$I, type = 'l', col = 'red')
