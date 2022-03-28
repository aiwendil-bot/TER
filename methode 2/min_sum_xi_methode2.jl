using JuMP
using GLPK
using SCIP
using Plots
using DelimitedFiles


function methode_2(first::Float64, last::Float64, m::Int64, x::Float64)

    m -= 1
    breakpoints=collect(LinRange(first,last,m+1))
    fpoints=map(x->x^3 - x,breakpoints)

    s = [(fpoints[k+1] - fpoints[k])/(breakpoints[k+1] - breakpoints[k]) for k in 1:m]
    concave = Vector{Int64}(undef,0)
    convexe = Vector{Int64}(undef,0)

    for k in 2:m
        if s[k]>s[k-1]
            push!(convexe,k)
        end
        if s[k]< s[k-1]
            push!(concave,k)
        end
    end

    model::Model = Model(GLPK.Optimizer)
    #model::Model = Model(SCIP.Optimizer)              M = m + 1

    @variable(model, d[1:m] >= 0)
    @variable(model, z[1:m] >= 0)
    @variable(model, u[1:m], Bin)

            @constraint(model, x + sum(d[l] for l in 1:(m-1)) >= breakpoints[m])
            for l in 1:(m-1)
                @constraint(model, d[l] <= breakpoints[l+1] - breakpoints[l])
            end    
        for k in concave
            @constraint(model, x + last * (u[k] - 1 ) <= z[k])
        end    
    #@constraint(model, sum(x[i]^2 for i in 1:n) <= 1)
    optimize!(model)

    return fpoints[1] + s[1]*(x - breakpoints[1]) + 2* sum( (s[k] - s[k-1])*(x - breakpoints[k] + sum(value(d[l]) for l in 1:(k-1))) for k in convexe) + sum( ((s[k] - s[k-1])*(x - 2*value(z[k]) + 2*breakpoints[k]*value(u[k]) - breakpoints[k])) for k in concave)
end

#m nombre de breakpoints, n nombre de variables dans la somme
function min_sum_xi(first::Float64, last::Float64, n::Int64, m::Int64)::Vector{Float64}

    m -= 1
    breakpoints=collect(LinRange(first,last,m+1))
    fpoints=map(x->x^2,breakpoints)

    s = [(fpoints[k+1] - fpoints[k])/(breakpoints[k+1] - breakpoints[k]) for k in 1:m]

    concave = Vector{Int64}(undef,0)
    convexe = Vector{Int64}(undef,0)

    for k in 2:m
        if s[k]>s[k-1]
            push!(convexe,k)
        end
        if s[k]< s[k-1]
            push!(concave,k)
        end
    end            

    model::Model = Model(GLPK.Optimizer)
    #model::Model = Model(SCIP.Optimizer)              M = m + 1

    @variable(model, first <= x[1:n] <= last)

    @variable(model, d[1:n, 1:m] >= 0)
    @variable(model, z[1:n, 1:m] >= 0)
    @variable(model, u[1:n, 1:m], Bin)

    @objective(model, Min, sum(x[i] for i in 1:n))

    for i in 1:n
        @constraint(model, x[i] + sum(d[i,l] for l in 1:(m-1)) >= breakpoints[m])
        for l in 1:(m-1)
            @constraint(model, d[i,l] <= breakpoints[l+1] - breakpoints[l])
        end    
        for k in concave
            @constraint(model, x[i] + last * (u[i,k] - 1 ) <= z[i,k])
        end    
    end    

    @constraint(model, sum( (fpoints[1] + s[1]*(x[i] - breakpoints[1]) + sum( (s[k] - s[k-1])*(x[i] - breakpoints[k] + sum(d[i,l] for l in 1:(k-1))) for k in convexe) + 0.5 * sum( (s[k] - s[k-1])*(x[i] - 2*z[i,k] + 2*breakpoints[k]*u[i,k] - breakpoints[k] ) for k in concave)) for i in 1:n)  <= 1)
    #@constraint(model, sum(x[i]^2 for i in 1:n) <= 1)
    optimize!(model)

        return value.(x)       
end