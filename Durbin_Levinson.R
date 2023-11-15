## Durbin-Levinson:
DurbinLevinson = function(serie, ar = NULL, ma = NULL){
X = serie
## vector nu = c(nu1, nu2, ....)
nu = c()
## vector de predicciones a un paso
X.hat = c()
## largo del vector
N = length(X)
## Matriz de coeficientes phi[i,j]
Phi = matrix(0, ncol = N, nrow = N)
rho = ARMAacf(ar = ar, ma = ma, lag.max = N)[2:(N+1)]
Phi[1,1] = rho[1]
nu[1] = (1-Phi[1,1]^2)
X.hat[1] = 0
X.hat[2] = Phi[1,1]*X[1]
for(n in 2:(N-1)){
	Phi[n,n] = (rho[n] - sum(Phi[n-1,1:(n-1)]*rho[n-1:(n-1)]))/nu[n-1]
	Phi[n,1:(n-1)] = Phi[n-1,1:(n-1)]-Phi[n,n]*Phi[n-1,(n-1):1]
	nu[n] = nu[n-1]*(1-Phi[n,n]^2)
X.hat[n+1] = sum(X[n:1]*Phi[n,1:n])
}
nu = ((N-1)*var(serie)/N)*c(1,nu)
list(fitted = X.hat, nu = nu)
}
