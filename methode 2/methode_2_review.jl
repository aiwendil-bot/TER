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

    model::Model = Model(GLPK.Optimizer)
    #model::Model = Model(SCIP.Optimizer)              M = m + 1
    @variable(model, d[1:(m-1)] >= 0)
    @variable(model, z[1:m-1])
    @variable(model, u[1:m-1], Bin)

    #contraintes partie convexe

    if !isempty(convexe)
        @constraint(model, x  + sum(d[l] for l in 1:(m -2)) >= breakpoints[m-1])
        for l in 1:(m-1)
            @constraint(model,d[l] <= breakpoints[l+1] - breakpoints[l])
        end
    end    
    #contraintes partie concave

    for k in concave 
        @constraint(model, (x  + breakpoints[k+1] * (u[k] - 1 )) <= z[k])
        @constraint(model, z[k] >= 0)
    end   
        
    
         
    #@constraint(model, sum(x ^2 for i in 1:n) <= 1)
    optimize!(model)

    L_f_x = fpoints[1] + s[1]*(x-breakpoints[1])
    
    if !isempty(convexe)
        L_f_x += sum(((s[k] - s[k-1])*(x - breakpoints[k] + sum(value(d[l]) for l in 1:(k-1)))) for k in convexe )
    end
    if !isempty(concave)    
        L_f_x += sum(((s[k] - s[k-1])*(x -2*value(z[k]) +2*breakpoints[k]*value(u[k]) - breakpoints[k] )) for k in concave)
    end    
    return L_f_x    
end

#bornes
a = 0.0
b = 0.4


#nb_breakpoints
m = 5
bornes = [a,b]

breakpoints = collect(LinRange(bornes[1],bornes[2],m))
#fbreakpoints = map(x->-x^2 + 3*x + 3,breakpoints)
#fpoints = map(x->-x^2 + 3*x + 3,[k for k in a:0.01:b])
fpoints = map(x->x->x^3-2*x^2-x+1,[k for k in a:0.01:b])

s= [(fbreakpoints[k+1] - fbreakpoints[k])/(breakpoints[k+1] - breakpoints[k]) for k in 1:(m-1)]

points_approches = Vector{Float64}(undef,0)
points_approx = Vector{Float64}(undef,0)


for x = a:0.01:b
    #push!(points_approx,approx_avec_val_absolues(x,breakpoints,fpoints,s))
    push!(points_approches,methode_2_review(x,5,[a,b]))

end


plot([k for k in 0:0.01:0.4],[points_approches,fpoints])