using JuMP
using GLPK
using SCIP
using Plots
using DelimitedFiles


function methode1(fonction,first::Float64, last::Float64, m::Int64)::Float64

    breakpoints=collect(LinRange(first,last,m+1))
    fpoints=map(fonction,breakpoints)

    model::Model = Model(GLPK.Optimizer)
    #model::Model = Model(SCIP.Optimizer)

    @variable(model, y[1:m], Bin)
    @variable(model, t[1:m+1] >= 0)
    @variable(model,z)
    @variable(model, first <= x <= last)
    @objective(model,Min,z)

        @constraint(model, t[1] <= y[1])
        for k in 2:m
            @constraint(model, t[k] <= y[k-1] + y[k])
        end
        @constraint(model, t[m+1] <= y[m])
        @constraint(model, sum(y[k] for k in 1:m) == 1)
        @constraint(model, sum(t[k] for k in 1:(m+1)) == 1)
        @constraint(model, sum(breakpoints[k]*t[k] for k in 1:(m+1)) == x)
  
    @constraint(model, z==  sum(fpoints[k] * t[k] for k in 1:(m+1)) )
    optimize!(model)
    return  sum(breakpoints[k] * value(t[k]) for k in 1:(m+1))    
end

function graphiques_erreurs()

    bornes::Vector{Float64} = [0.0,3.0]

    f(x::Float64) = (x-2)^2
    g(x::Float64) = -(x-2)^2 + 3
    h(x::Float64) = (x-1)^3 -x^2 + 0.5x +2

    minimums::Vector{Float64} = [2.0,0.0,2.12]


    fonctions::Vector{Function} = [f,g,h]
    nb_breakpoints::Vector{Int64} = [k for k in 3:500]

    distance_minimum = Array{Float64, 2}(undef,length(nb_breakpoints),length(fonctions))
    temps = Array{Float64, 2}(undef,length(nb_breakpoints),length(fonctions))

    for i in 1:length(fonctions)

        for j in 1:length(nb_breakpoints)
            distance_minimum[j,i] = abs(minimums[i] - methode1(fonctions[i],bornes[1],bornes[2],nb_breakpoints[j]))
            
            temps[j,i] = @elapsed methode1(fonctions[i],bornes[1],bornes[2],nb_breakpoints[j])
        end
    end

    writedlm("methode 1/distance_minimum.txt", distance_minimum)
    writedlm("methode 1/temps.txt", temps)

    labels = ["f" "g" "h"]

    pdistance = plot(nb_breakpoints,distance_minimum, label=labels)
    ptemps = plot(nb_breakpoints,temps, label=labels,legend=:topleft)
    savefig(pdistance,"methode 1/distance_minimum.png")
    savefig(ptemps,"methode 1/temps.png")

end