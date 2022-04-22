using Plots
using DelimitedFiles

function pl_square(x::Float64, breakpoints::Vector{Float64}, fpoints::Vector{Float64}, s::Vector{Float64})::Float64

    return fpoints[1] + s[1]*(x-breakpoints[1]) + sum((s[k] - s[k-1])/2 * (abs(x-breakpoints[k]) + x - breakpoints[k]) for k = 2:length(s))

end    



function cercle_approche(first::Float64, last::Float64, m::Int64)::Vector{Vector{Float64}}

    breakpoints::Vector{Float64} =collect(LinRange(first,last,m+1))
    fpoints::Vector{Float64} =map(x->x^2,breakpoints)

    s = [(fpoints[k+1] - fpoints[k])/(breakpoints[k+1] - breakpoints[k]) for k in 1:m]

    approx = [pl_square(x,breakpoints,fpoints,s) for x in first:0.0001:last]

    res = Vector{Vector{Float64}}(undef,0)

    for x = first:0.001:last, y = first:0.001:last
        if abs(pl_square(x,breakpoints,fpoints,s) + pl_square(y,breakpoints,fpoints,s)-1) <= 10^(-3)
            push!(res,[x,y])
        end
    end
    
    return res

end     

function genere_cercles_approches()
    
    for m in 2:20
        res = cercle_approche(-2.0,2.0,m)
        writedlm("points_cercle_approche_$m.txt", res)

        x = Vector{Float64}(undef,length(res))
        y = Vector{Float64}(undef,length(res))

        for k in 1:length(res)
            x[k] = res[k][1]
            y[k] = res[k][2]
        end

        p = plot(x,y,seriestype = :scatter)
        plot!(size=(400,400))
        savefig(p, "autres problèmes/cercle_approche_$m.png")
    end    


end
