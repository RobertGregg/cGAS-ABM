using DifferentialEquations, DiffEqOperators

N = 5
dx = 1.0
Dxx = CenteredDifference(2,2,dx,N)
Dyy = CenteredDifference{2}(2,2,dx,N)
Dzz = CenteredDifference{3}(2,2,dx,N)
L = Dxx + Dyy + Dzz

u = rand(N,N,N)
Q = Neumann0BC(dx,1)
sol = L*Q*u
