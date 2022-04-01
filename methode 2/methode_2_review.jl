using JuMP
using GLPK
using SCIP
using Plots
using DelimitedFiles

#m = nb_breakpoints
function methode_2_review(x::Float64,m::Int64,bornes::Vector{Float64})

    breakpoints=collect(LinRange(bornes[1],bornes[2],m))
    fpoints= map(x->x^3-2*x^2-x+1,breakpoints)

    s::Vector{Float64} = [(fpoints[k+1] - fpoints[k])/(breakpoints[k+1] - breakpoints[k]) for k in 1:(m-1)]
    concave::Vector{Int64} = Vector{Int64}(undef,0)
    convexe::Vector{Int64} = Vector{Int64}(undef,0)

    for k in 2:(m -1)
        if s[k]>s[k-1]
            push!(convexe, k)
        end
        if s[k]< s[k-1]
            push!(concave, k)
        end
    end

    println(breakpoints)
    println(convexe)
    println(concave)
    println(s)

    model::Model = Model(GLPK.Optimizer)
    #model::Model = Model(SCIP.Optimizer)              M = m + 1
    @variable(model, d[1:(m-1)] >= 0)
    @variable(model, z[1:m-1])
    @variable(model, u[1:m-1], Bin)

    #contraintes partie convexe
    @constraint(model, x  + sum(d[l] for l in 1:(m -2)) >= breakpoints[m-1])
    for l in 1:(m -1)
        @constraint(model,d[l] <= breakpoints[l+1] - breakpoints[l])
    end

    #contraintes partie concave

    for k in concave 
        println(k)
        @constraint(model, (x  + bornes[2] * (u[k] - 1 )) <= z[k])
        @constraint(model, z[k] >= 0)
    end   
        
    
         
    #@constraint(model, sum(x ^2 for i in 1:n) <= 1)
    optimize!(model)

    L_f_x = fpoints[1] + s[1]*(x-breakpoints[1])
    
    if !isempty(convexe)
        L_f_x += sum(((s[k] - s[k-1])*(x - breakpoints[k] + sum(value(d[l]) for l in 1:(k-1)))) for k in convexe )
    end
    if !isempty(concave)    
        L_f_x += 0.5*sum(((s[k] - s[k-1])*(x -2*value(z[k]) +2*breakpoints[k]*value(u[k]) - breakpoints[k] )) for k in concave)
    end    
    return L_f_x    
end

#points_approches = Vector{Float64}(undef,0)

#for x in 1:0.01:2
#    push!(points_approches,methode_2_review(x,5,[1.0,2.0]))
#end
#f(x) = x->x^3-2*x^2-x+1

#plot([k for k in 1:0.01:2],[points_approches,f])