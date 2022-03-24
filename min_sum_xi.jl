using JuMP
using GLPK
using SCIP
using Plots
using DelimitedFiles
using MathOptInterface

#m nombre de breakpoints, n nombre de variables dans la somme
function min_sum_xi(first::Float64, last::Float64, n::Int64, m::Int64)::Vector{Float64}

    m -= 1
    breakpoints::Vector{Float64} =collect(LinRange(first,last,m+1))
    fpoints::Vector{Float64} =map(x->x^2,breakpoints)

    model::Model = Model(GLPK.Optimizer)
    #model::Model = Model(SCIP.Optimizer)

    @variable(model, y[1:n, 1:m], Bin)
    @variable(model, t[1:n, 1:m+1] >= 0)

    @variable(model, first <= x[1:n] <= last)
    @objective(model, Min, sum(x[i] for i in 1:n))

    for i in 1:n
        @constraint(model, t[i,1] <= y[i,1])
        for k in 2:m
            @constraint(model, t[i,k] <= y[i,k-1] + y[i,k])
        end
        @constraint(model, t[i,m+1] <= y[i,m])
        @constraint(model, sum(y[i,k] for k in 1:m) == 1)
        @constraint(model, sum(t[i,k] for k in 1:(m+1)) == 1)
        @constraint(model, sum(breakpoints[k]*t[i,k] for k in 1:(m+1)) == x[i])
    end    

    @constraint(model, sum(sum(fpoints[k] * t[i,k] for k in 1:(m+1)) for i in 1:n) <= 1)
    #@constraint(model, sum(x[i]^2 for i in 1:n) <= 1)

    write_to_file(model, "focus/$n"*"_"*"$m.lp", format = MOI.FileFormats.FORMAT_LP)

    optimize!(model)

    return value.(x)       
end

function graphiques_erreurs()

    a = min_sum_xi(-2.0,2.0,5,5)

    nb_valeurs::Vector{Int64} = [2,3,4,5,10,15,20]
    nb_breakpoints::Vector{Int64} = [k for k in 5:5:500]

    #valeurs = Array{Float64, 2}(undef,length(nb_breakpoints),length(nb_valeurs))
    temps = Array{Float64, 2}(undef,length(nb_breakpoints),length(nb_valeurs))

    for i in 1:length(nb_valeurs)

        for j in 1:length(nb_breakpoints)
            #valeurs[j,i] = sqrt(sum((min_sum_xi(-2.0,2.0,nb_valeurs[i],nb_breakpoints[j]) .- [-1/sqrt(nb_valeurs[i]) for k in 1:nb_valeurs[i]]).^2))
            temps[j,i] = @elapsed min_sum_xi(-2.0,2.0,nb_valeurs[i],nb_breakpoints[j])
        end
    end

    #writedlm("resultats_calculs.txt", valeurs)
    writedlm("glpk_jump_temps.txt", temps)
    labels = ["n = 2" "n = 3" "n = 4" "n = 5" "n = 10" "n = 15" "n = 20" ]
    ptempsloglog = plot(nb_breakpoints, temps, title="Min sum x_i, sum x_i^2 <= 1", label=labels,xlabel="nb_breakpoints", ylabel="temps execution", yaxis=:log, xaxis=:log,legend=:outerbottomright)
    ptempslog = plot(nb_breakpoints, temps, title="Min sum x_i, sum x_i^2 <= 1", label=labels,xlabel="nb_breakpoints", ylabel="temps execution", yaxis=:log,legend=:outerbottomright)
    
    ptemps = plot(nb_breakpoints, temps, title="Min sum x_i, sum x_i^2 <= 1", label=labels,xlabel="nb_breakpoints", ylabel="temps execution",legend=:outerbottomright)
    savefig(ptempsloglog, "glpk_jump_tempsloglog.png")

    savefig(ptempslog, "glpk_jump_tempslog.png")
    savefig(ptemps, "glpk_jump_temps.png")

end

function focus()

    temps = Vector{Float64}(undef,11)
    breakpoints = [k for k = 320:330]
    for k = 1:11
        temps[k] = @elapsed min_sum_xi(-2.0,2.0,20,breakpoints[k])
    end
    writedlm("focus/temps.txt", temps)
    ptemps = plot(breakpoints, temps, title="Min sum x_i, sum x_i^2 <= 1", label=labels,xlabel="nb_breakpoints", ylabel="temps execution",legend=:outerbottomright)
    savefig(ptemps, "focus/temps.png")

end
labels = ["n = 2" "n = 3" "n = 4" "n = 5" "n = 10" "n = 15" "n = 20" ]
breakpoints = [k for k = 320:330]
temps = readdlm("focus/temps.txt")
ptemps = plot(breakpoints, temps, title="Min sum x_i, sum x_i^2 <= 1", label=labels,xlabel="nb_breakpoints", ylabel="temps execution",legend=:outerbottomright)
savefig(ptemps, "focus/temps.png")