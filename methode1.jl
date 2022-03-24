using JuMP
using GLPK
using SCIP
using Plots
using DelimitedFiles


function methode1(first::Float64, last::Float64, m::Int64, x::Float64)::Float64

    m -= 1
    breakpoints=collect(LinRange(first,last,m+1))
    println(breakpoints)
    fpoints=map(x->x^3 - x,breakpoints)

    model::Model = Model(GLPK.Optimizer)
    #model::Model = Model(SCIP.Optimizer)

    @variable(model, y[1:m], Bin)
    @variable(model, t[1:m+1] >= 0)

        @constraint(model, t[1] <= y[1])
        for k in 2:m
            @constraint(model, t[k] <= y[k-1] + y[k])
        end
        @constraint(model, t[m+1] <= y[m])
        @constraint(model, sum(y[k] for k in 1:m) == 1)
        @constraint(model, sum(t[k] for k in 1:(m+1)) == 1)
        @constraint(model, sum(breakpoints[k]*t[k] for k in 1:(m+1)) == x)
  
    #@constraint(model, sum(x[i]^2 for i in 1:n) <= 1)
    optimize!(model)

        return sum(fpoints[k] * value(t[k]) for k in 1:(m+1))      
end