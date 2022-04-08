using JuMP
using GLPK
using SCIP
using Plots
using DelimitedFiles

#m = nb_breakpoints
function methode_2_review_min(x::Float64,m::Int64,bornes::Vector{Float64})

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
    #M = m + 1
    @variable(model, d[1:(m-1)] >= 0)
    @variable(model, z[1:m-1]>=0)
    @variable(model, u[1:m-1], Bin)
    @variable(model, y)
    @objective(model,Min,y)

    #contraintes partie convexe

    if !isempty(convexe)
        @constraint(model, x  + sum(d[l] for l in 1:(m -2)) >= breakpoints[m-1])
        for l in 1:(m-1)
            @constraint(model,d[l] <= breakpoints[l+1] - breakpoints[l])
        end
    end    
    #contraintes partie concave

    for k in concave 
        @constraint(model, (x  + bornes[2]*(u[k] - 1 )) <= z[k])
    end   
    

    L_f_x = fpoints[1] + s[1]*(x-breakpoints[1])
    
    if !isempty(convexe)
        L_f_x += sum(((s[k] - s[k-1])*(x - breakpoints[k] + sum(d[l] for l in 1:(k-1)))) for k in convexe )
    end
    if !isempty(concave)    
        L_f_x += 0.5*sum(((s[k] - s[k-1])*(2*x -2*z[k] +2*breakpoints[k]*u[k] - 2*breakpoints[k] )) for k in concave)
    end
    @constraint(model, y == L_f_x)

    optimize!(model)
    if x == 0.0
        println(model)
        print(value.(u))
        print(value.(z))
        end
    return value(y)
end


function approx_avec_val_absolues(x::Float64, breakpoints::Vector{Float64}, fpoints::Vector{Float64}, s::Vector{Float64})::Float64

    return fpoints[1] + s[1]*(x-breakpoints[1]) + sum((s[k] - s[k-1])/2 * (abs(x-breakpoints[k]) + x - breakpoints[k]) for k = 2:length(s))

end


points_approches = Vector{Float64}(undef,0)

#bornes
a = 0.0
b = 2.0


#nb_breakpoints
m = 10
bornes = [a,b]

breakpoints = collect(LinRange(a,b,m))
fbreakpoints = map(x->x^3-2*x^2-x+1,breakpoints)
#fpoints = map(x->-x^2 + 3*x + 3,[k for k in a:0.01:b])
fpoints = map(x->x^3-2*x^2-x+1,[k for k in a:0.01:b])

println(breakpoints)
println(fbreakpoints)

s= [(fbreakpoints[k+1] - fbreakpoints[k])/(breakpoints[k+1] - breakpoints[k]) for k in 1:(m-1)]
println([s[k]-s[k-1] for k in 2:(length(s))])
points_approches = Vector{Float64}(undef,0)
points_approx = Vector{Float64}(undef,0)


for x = a:0.01:b
    push!(points_approx,approx_avec_val_absolues(x,breakpoints,fpoints,s))
    push!(points_approches,methode_2_review_min(x,m,[a,b]))

end

plot([k for k in a:0.01:b],[points_approches,fpoints],legend=:outerbottomright)
