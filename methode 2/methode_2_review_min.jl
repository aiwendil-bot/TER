using JuMP
using GLPK
using SCIP
using Plots
using DelimitedFiles

#m = nb_breakpoints
function methode_2_review_min(fonction,m::Int64,bornes::Vector{Float64})::Float64

    breakpoints=collect(LinRange(bornes[1],bornes[2],m))
    fpoints= map(fonction,breakpoints)

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
    set_optimizer_attribute(model,"presolve", GLP_ON)

    #M = m + 1
    @variable(model, d[1:(m-1)] >= 0)
    @variable(model, z[1:m-1]>=0)
    @variable(model, u[1:m-1], Bin)
    @variable(model, bornes[1] <= x <= bornes[2])
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

    return value(x)
end


function approx_avec_val_absolues(x::Float64, breakpoints::Vector{Float64}, fpoints::Vector{Float64}, s::Vector{Float64})::Float64

    return fpoints[1] + s[1]*(x-breakpoints[1]) + sum((s[k] - s[k-1])/2 * (abs(x-breakpoints[k]) + x - breakpoints[k]) for k = 2:length(s))

end


points_approches = Vector{Float64}(undef,0)

#bornes
a = 0
b = 3


#nb_breakpoints
m = 8
bornes = [a,b]

breakpoints = collect(LinRange(a,b,m))
fbreakpoints = map(x->-x^3 + 2*x^2+x-1,breakpoints)
#fpoints = map(x->-x^2 + 3*x + 3,[k for k in a:0.01:b])
fpoints = map(x->-x^3 + 2*x^2+x-1,[k for k in a:0.01:b])


s= [(fbreakpoints[k+1] - fbreakpoints[k])/(breakpoints[k+1] - breakpoints[k]) for k in 1:(m-1)]
points_approches = Vector{Float64}(undef,0)
points_approx = Vector{Float64}(undef,0)


#for x = a:0.01:b
 #   push!(points_approx,approx_avec_val_absolues(x,breakpoints,fpoints,s))
  #  push!(points_approches,methode_2_review_min(x,m,[a,b]))

#end


function graphiques_erreurs()

    bornes::Vector{Float64} = [0.0,3.0]

    f(x::Float64) = (x-2)^2
    g(x::Float64) = -(x-2)^2 + 3
    h(x::Float64) = (x-1)^3 -x^2 + 0.5x +2

    minimums::Vector{Float64} = [2.0,0.0,2.12]


    fonctions::Vector{Function} = [f,g,h]
    nb_breakpoints::Vector{Int64} = [k for k in 3:150]

    distance_minimum = Array{Float64, 2}(undef,length(nb_breakpoints),length(fonctions))
    temps = Array{Float64, 2}(undef,length(nb_breakpoints),length(fonctions))

    for i in 1:length(fonctions)

        for j in 1:length(nb_breakpoints)
            distance_minimum[j,i] = abs(minimums[i] - methode_2_review_min(fonctions[i],nb_breakpoints[j],bornes))
            
            temps[j,i] = @elapsed methode_2_review_min(fonctions[i],nb_breakpoints[j],bornes)
        end
    end

    writedlm("methode 2/distance_minimum.txt", distance_minimum)
    writedlm("methode 2/temps.txt", temps)

    labels = ["f" "g" "h"]

    pdistance = plot(nb_breakpoints,distance_minimum, label=labels)
    ptemps = plot(nb_breakpoints,temps, label=labels,legend=:topleft)
    savefig(pdistance,"methode 2/distance_minimum.png")
    savefig(ptemps,"methode 2/temps.png")

end

