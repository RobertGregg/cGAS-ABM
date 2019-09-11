using DifferentialEquations, DiffEqOperators, PLots

# 5-point stencil
#Need to add the diffusion coefficient and the grid size
function ∇²(u)
    Δu = similar(u)
    n1, n2 = size(u)
    h = 1/(30.0)^2 #Grid spacing (diameter of cell in μm)

    # internal nodes
    for j = 2:n2-1
        for i = 2:n1-1
            @inbounds  Δu[i,j] = h*(u[i+1,j] + u[i-1,j] + u[i,j+1] + u[i,j-1] - 4*u[i,j])
        end
    end

    # left/right edges
    for i = 2:n1-1
        @inbounds Δu[i,1] = h*(u[i+1,1] + u[i-1,1] + 2*u[i,2] - 4*u[i,1])
        @inbounds Δu[i,n2] = h*(u[i+1,n2] + u[i-1,n2] + 2*u[i,n2-1] - 4*u[i,n2])
    end

    # top/bottom edges
    for j = 2:n2-1
        @inbounds Δu[1,j] = h*(u[1,j+1] + u[1,j-1] + 2*u[2,j] - 4*u[1,j])
        @inbounds Δu[n1,j] = h*(u[n1,j+1] + u[n1,j-1] + 2*u[n1-1,j] - 4*u[n1,j])
    end

    # corners
    @inbounds Δu[1,1] = h*(2*(u[2,1] + u[1,2]) - 4*u[1,1])
    @inbounds Δu[n1,1] = h*(2*(u[n1-1,1] + u[n1,2]) - 4*u[n1,1])
    @inbounds Δu[1,n2] = h*(2*(u[2,n2] + u[1,n2-1]) - 4*u[1,n2])
    @inbounds Δu[n1,n2] = h*(2*(u[n1-1,n2] + u[n1,n2-1]) - 4*u[n1,n2])

    return Δu
end


N = 100

u0 = zeros(N,N)
u0[1:50,1:50] .= 30.0
tspan = (0.0,200.0)

function f(du,u,p,t)
    du = ∇²(u)
end


prob = ODEProblem(f, u0, tspan)
sol = solve(prob,ROCK4())


#Loop through all time points to make an animation
Maxval = maximum(u0)
anim = @animate for i = tspan[1]:tspan[2]
    G = sol(i)
    plot(diag(G),ylims=(0,Maxval),title="Time = " * string(i) * " hrs",legend=false)
end

gif(anim,"CellGridAnimation.gif")
