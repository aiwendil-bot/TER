using Plots
using DelimitedFiles
using ImplicitEquations
function pl_square(x::Float64, breakpoints::Vector{Float64}, fpoints::Vector{Float64}, s::Vector{Float64})::Float64

    return fpoints[1] + s[1]*(x-breakpoints[1]) + sum((s[k] - s[k-1])/2 * (abs(x-breakpoints[k]) + x - breakpoints[k]) for k = 2:length(s))

end    



function cercle_approche(first::Float64, last::Float64, m::Int64)::Vector{Vector{Float64}}

    m -= 1
    breakpoints::Vector{Float64} =collect(LinRange(first,last,m+1))
    fpoints::Vector{Float64} =map(x->x^2,breakpoints)

    s = [(fpoints[k+1] - fpoints[k])/(breakpoints[k+1] - breakpoints[k]) for k in 1:m]

    res = Vector{Vector{Float64}}(undef,0)

    display(p)
    for x = first:0.0001:last, y = first:0.0001:last
        if abs(pl_square(x,breakpoints,fpoints,s) + pl_square(y,breakpoints,fpoints,s)-1) <= 10^(-3)
            push!(res,[x,y])
        end
    end
    
    return res

end    

res = cercle_approche(-2.0,2.0,20)
writedlm("points_cercle_approche.txt", res)

x = Vector{Float64}(undef,length(res))
y = Vector{Float64}(undef,length(res))

for k in 1:length(res)
    x[k] = res[k][1]
    y[k] = res[k][2]
end

p = plot(x,y)

savefig(p, "cercle_approche.png")