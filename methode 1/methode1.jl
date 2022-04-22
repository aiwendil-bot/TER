using JuMP
using GLPK
using SCIP
using Plots
using DelimitedFiles
using LaTeXStrings


function methode1(first::Float64, last::Float64, m::Int64, x::Float64)::Float64

    breakpoints=collect(LinRange(first,last,m+1))
    fpoints=map(x-> x^2,breakpoints)

    model::Model = Model(GLPK.Optimizer)
    set_optimizer_attribute(model,"presolve",GLP_ON)
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

#points = [methode1(-2.0,2.0, 5, x) for x in -2.0:0.01:2.0]

#fp = [x^2 for x in -2.0:0.01:2.0]

#p2 = plot([x for x in -2.0:0.01:2.0], [points,fp],label=[L"PL_{sqr,5}" L"y=x^2"],legend=:bottomright)