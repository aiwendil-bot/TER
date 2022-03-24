using Plots
using DelimitedFiles

function recup_scip()

    nb_valeurs::Vector{Int64} = [20]
    nb_breakpoints::Vector{Int64} = [k for k in 320:330]

    temps = Array{Float64, 2}(undef,length(nb_breakpoints),length(nb_valeurs))

    rx = r"^[0-9]{1,2}_[0-9]{1,3}"
    rx2 = r"[0-9]{1,2}\.[0-9]{1,2}$"

    open("focus/scip/output_temps.txt") do file

        for ln in eachline(file)
            m1 = match(rx,ln)
            m2 = match(rx2, ln)
            nm = split(m1.match,'_')
            n= parse(Int64,nm[1])
            m= parse(Int64,nm[2])
            time = parse(Float64,m2.match)
            println("n = $n")
            println("m = $m")
            temps[findfirst(x->x==m, nb_breakpoints),findfirst(x->x==n, nb_valeurs)] = time
            
        end
    end

    writedlm("focus_scip_temps.txt", temps)
    #labels = ["n = 2" "n = 3" "n = 4" "n = 5" "n = 10" "n = 15" "n = 20" ]
    labels = ["n = 20" ]
    #ptempsloglog = plot(nb_breakpoints, temps, title="Min sum x_i, sum x_i^2 <= 1", label=labels,xlabel="nb_breakpoints", ylabel="temps execution", yaxis=:log, xaxis=:log)
    #ptempslog = plot(nb_breakpoints, temps, title="Min sum x_i, sum x_i^2 <= 1", label=labels,xlabel="nb_breakpoints", ylabel="temps execution", yaxis=:log)
    
    ptemps = plot(nb_breakpoints, temps, title="Min sum x_i, sum x_i^2 <= 1", label=labels,xlabel="nb_breakpoints", ylabel="temps execution",legend=:topleft)
    #savefig(ptempsloglog, "scip_tempsloglog.png")

    #savefig(ptempslog, "scip_tempslog.png")
    savefig(ptemps, "focus/scip/scip_temps.png")


end

recup_scip()
