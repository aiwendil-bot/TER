using Plots
using DelimitedFiles

function recup_glpk()

    nb_valeurs::Vector{Int64} = [2,3,4,5,10,15,20]
    nb_breakpoints::Vector{Int64} = [k for k in 4:5:500]

    temps = Array{Float64, 2}(undef,length(nb_breakpoints),length(nb_valeurs))

    rx = r"^[0-9]{1,2}_[0-9]{1,3}"
    rx2 = r"[0-9]{1,2}\.[0-9]"

    open("res_glpk/output_temps.txt") do file

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

    writedlm("glpk_temps.txt", temps)
    #labels = ["n = 2" "n = 3" "n = 4" "n = 5" "n = 10" "n = 15" "n = 20" ]
    #ptempslog = plot(nb_breakpoints, temps, title="Min sum x_i, sum x_i^2 <= 1", label=labels,xlabel="nb_breakpoints", ylabel="temps execution", yaxis=:log, xaxis=:log)
    #ptemps = plot(nb_breakpoints, temps, title="Min sum x_i, sum x_i^2 <= 1", label=labels,xlabel="nb_breakpoints", ylabel="temps execution", yaxis=:log)
    #savefig(ptempslog, "glpk_tempsloglog.png")
    #savefig(ptemps, "glpk_temps.png")


end

recup_glpk()
nb_breakpoints = [k for k in 5:5:500]
labels = ["n = 2" "n = 3" "n = 4" "n = 5" "n = 10" "n = 15" "n = 20" ]

ptemps = plot(nb_breakpoints, readdlm("glpk_temps.txt",'\t',Float64,'\n'),label=labels,xlabel="nb_breakpoints", ylabel="temps d'execution", legend=:outerbottomright)
savefig(ptemps, "glpk_temps.png")

#m =match(rx,"20_389.log:Time used:   10.6 secs")
#m =match(rx2,"n_20_m_384.log:Time used:   9.9 secs")