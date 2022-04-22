using JuMP
using GLPK
using SCIP
using Plots
using DelimitedFiles
using MathOptInterface
using LaTeXStrings
#m nombre de breakpoints, n nombre de variables dans la somme
function min_sum_xi(first::Float64, last::Float64, n::Int64, m::Int64)::Vector{Float64}

    m -= 1
    breakpoints::Vector{Float64} =collect(LinRange(first,last,m+1))
    fpoints::Vector{Float64} =map(x->x^2,breakpoints)

    #model::Model = Model(GLPK.Optimizer)
    model::Model = Model(SCIP.Optimizer)
    #set_optimizer_attribute(model,"msg_lev", GLPK.GLP_MSG_ALL)
    #set_optimizer_attribute(model,"presolve", GLP_ON)

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

    #@constraint(model, sum(sum(fpoints[k] * t[i,k] for k in 1:(m+1)) for i in 1:n) <= 1)
    @constraint(model, sum(x[i]^2 for i in 1:n) <= 1)

    #write_to_file(model, "focus/$n"*"_"*"$m.lp", format = MOI.FileFormats.FORMAT_LP)

    optimize!(model)

    return value.(x)       
end

function graphiques_erreurs()

    a = min_sum_xi(-2.0,2.0,5,5)

    nb_valeurs::Vector{Int64} = [2,3,4,5,10,15,20]
    nb_breakpoints::Vector{Int64} = [k for k in 5:5:300]

    valeurs = Array{Float64, 2}(undef,length(nb_breakpoints),length(nb_valeurs))
    temps = Array{Float64, 2}(undef,length(nb_breakpoints),length(nb_valeurs))

    for i in 1:length(nb_valeurs)

        for j in 1:length(nb_breakpoints)
            valeurs[j,i] = sqrt(sum((min_sum_xi(-2.0,2.0,nb_valeurs[i],nb_breakpoints[j]) .- [-1/sqrt(nb_valeurs[i]) for k in 1:nb_valeurs[i]]).^2))
            temps[j,i] = @elapsed min_sum_xi(-2.0,2.0,nb_valeurs[i],nb_breakpoints[j])
        end
    end

    writedlm("methode 1/resultats_calculs.txt", valeurs)
    writedlm("methode 1/temps.txt", temps)
    labels = ["n = 2" "n = 3" "n = 4" "n = 5" "n = 10" "n = 15" "n = 20" ]
    ptempsloglog = plot(nb_breakpoints, temps, title=L"Min sum x_i, PL_{sqr} <= 1", label=labels,xlabel="nb_breakpoints", ylabel="temps execution", yaxis=:log, xaxis=:log,legend=:outerbottomright)
    ptempslog = plot(nb_breakpoints, temps, title=L"Min sum x_i, PL_{sqr} <= 1", label=labels,xlabel="nb_breakpoints", ylabel="temps execution", yaxis=:log,legend=:outerbottomright)
    
    ptemps = plot(nb_breakpoints, temps, title=L"Min sum x_i, PL_{sqr} <= 1", label=labels,xlabel="nb_breakpoints", ylabel="temps execution",legend=:outerbottomright)
    savefig(ptempsloglog, "methode 1/glpk_jump_tempsloglog_presolve.png")

    savefig(ptempslog, "methode 1/glpk_jump_tempslog_presolve.png")
    savefig(ptemps, "methode 1/glpk_jump_temps_presolve.png")

    pdistanceloglog = plot(nb_breakpoints, valeurs, title=L"Min sum x_i, PL_{sqr} <= 1", label=labels,xlabel="nb_breakpoints", ylabel="distance à la solution exacte", yaxis=:log, xaxis=:log,legend=:outerbottomright)
    pdistancelog = plot(nb_breakpoints, valeurs, title=L"Min sum x_i, PL_{sqr} <= 1", label=labels,xlabel="nb_breakpoints", ylabel="distance à la solution exacte", yaxis=:log,legend=:outerbottomright)
    
    pdistance = plot(nb_breakpoints, valeurs, title=L"Min sum x_i, PL_{sqr} <= 1", label=labels,xlabel="nb_breakpoints", ylabel="distance à la solution exacte",legend=:outerbottomright)
    savefig(pdistanceloglog, "methode 1/valeursloglog.png")

    savefig(pdistancelog, "methode 1/valeurslog.png")
    savefig(pdistance, "methode 1/valeurs.png")

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
labels = ["n = 2" "n = 3" "n = 5" "n = 10" "n = 15" "n = 20" ]
breakpoints = [k for k = 4:5:500]
distance = readdlm("res_jump/resultats_calculs.txt")
pdistance = plot(breakpoints, distance[:,[1,2,4,5,6,7]], title="Min sum x_i, sum x_i^2 <= 1", label=labels,xlabel="nb_breakpoints", ylabel="distance à la valeur exacte",yaxis=:log,xaxis=:log,legend=:outerbottomright)
savefig(pdistance, "res_jump/distancevaleurexacte_log_log.png")
