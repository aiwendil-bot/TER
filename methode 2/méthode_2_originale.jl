using JuMP
using GLPK
using SCIP
using Plots
using DelimitedFiles

#m = nb_breakpoints
function methode_2(functions,m::Vector{Int64},bornes::Vector{Vector{Float64}},upper_bound_x::Float64)

    n = length(functions)
    

    breakpoints=[collect(LinRange(bornes[i][1],bornes[i][2],m[i])) for i in 1:n]
    fpoints= [map(functions[i],breakpoints[i]) for i in 1:n]

    s::Vector{Vector{Float64}} = [[(fpoints[i][k+1] - fpoints[i][k])/(breakpoints[i][k+1] - breakpoints[i][k]) for k in 1:(m[i]-1)] for i in 1:n]
    concave::Vector{Vector{Int64}} = [Vector{Int64}(undef,0) for _ in 1:n]
    convexe::Vector{Vector{Int64}} = [Vector{Int64}(undef,0) for _ in 1:n]

    for i in 1:n
        for k in 2:(m[i]-1)
            if s[i][k]>s[i][k-1]
                push!(convexe[i],k)
            end
            if s[i][k]< s[i][k-1]
                push!(concave[i],k)
            end
        end
    end
    println(s)
    println(breakpoints[1])
    println(fpoints)
    println(convexe)
    println(concave)
    model::Model = Model(GLPK.Optimizer)
    #model::Model = Model(SCIP.Optimizer)              M = m + 1
    @variable(model,bornes[1][1] <= x[1:n] <= bornes[1][2])
    @variable(model, d[1:n,1:(maximum(m)-1)] >= 0)
    @variable(model, z[1:n,1:maximum(m)] >= 0)
    @variable(model, u[1:n,1:maximum(m)], Bin)

    @objective(model,Min,sum((fpoints[i][1] + s[i][1]*(x[i]-breakpoints[i][1]) + 2*sum(((s[i][k] - s[i][k-1])*(x[i] - breakpoints[i][k] 
    + sum(d[i,l] for l in 1:(k-1)))) for k in convexe[i]) + sum(((s[i][k] - s[i][k-1])*(x[i]-2*z[i,k] +2*breakpoints[i][k]*u[i,k] - breakpoints[i][k] ))
     for k in concave[i])) for i in 1:n))
    for i in 1:n
        #contraintes partie convexe
        @constraint(model, x[i] + sum(d[i,l] for l in 1:(m[i]-2)) >= breakpoints[i][m[i]-1])
        for l in 1:(m[i]-1)
            @constraint(model,d[i,l] <= breakpoints[i][l+1] - breakpoints[i][l])
        end

        #contraintes partie concave

        for k in concave[i]
            println(k)
            @constraint(model, (x[i] + upper_bound_x*(u[i,k] - 1 )) <= z[i,k])
            @constraint(model, z[i,k] >= 0)
        end   

    end        

    write_to_file(model, "methode2.lp", format = MOI.FileFormats.FORMAT_LP)    
    println(model)     
    #@constraint(model, sum(x[i]^2 for i in 1:n) <= 1)
    optimize!(model)

    return value.(x)
end

