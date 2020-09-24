module BosonSampling

# IMPORTS
using Immanants
import RandomMatrices: Haar
import StatsBase: sample, Weights
import IterTools: subsets
# EXPORTS
export bosonsampling, per, Haar
export generadistribucion

function per(a::Array{Float64,2})
    _,n::Int64 = size(a)
    immanant(Partition([n]),a)
end

function per(a::Array{Complex{Float64},2})
    _,n::Int64 = size(a)
    immanant(Partition([n]),a)
end

function pasoalgoritmo(k::Int64, mat::Array{Float64,2}, rp::Array{Int64,1})
    @assert length(rp) == k - 1
    _,n::Int64 = size(mat)

    probas::Array{Float64,1} = zeros(Float64,n)
    for i in 1:n
        total::Float64 = 0.0
        r = [rp; i]
        for c in subsets(1:n, k)
            total += abs( per( mat[r,c] ) )^2
        end
        probas[i] = total
    end
    probas/n
end

function pasoalgoritmo(k::Int64, mat::Array{Complex{Float64},2}, rp::Array{Int64,1})
    @assert length(rp) == k - 1
    _,n::Int64 = size(mat)

    probas::Array{Float64,1} = zeros(Float64,n)
    for i in 1:n
        total::Float64 = 0.0
        r = [rp; i]
        for c in subsets(1:n, k)
            total += Float64(abs( per( mat[r,c] ) )^2)
        end
        probas[i] = total
    end
    probas/n
end

function bosonsampling(mat::Array{Complex{Float64},2})
    _,n::Int64 = size(mat)
    proba::Array{Float64,1} = zeros(Float64, n)
    lista::Array{Int64,1} = Int64[]
    for i in 1:n
        proba = pasoalgoritmo(i, mat, lista)
        push!(lista, sample(1:n, proba |> Weights))
    end
    lista
end

function bosonsampling(mat::Array{Complex{Float64},2},m::Int64)
    _,n::Int64 = size(mat)
    proba::Array{Float64,1} = zeros(Float64, n)
    lista::Array{Int64,1} = Int64[]
    for i in 1:m
        proba = pasoalgoritmo(i, mat, lista)
        push!(lista, sample(1:n, proba |> Weights))
    end
    lista
end

#function estimadorP(mat::Array{Complex{Float64},2}, n::Int64, input::Array{Int64,1}, filas::Array{Int64,1})
#    @assert length(input) == length(filas)
#    tmp = mat[input, filas]
#    total::Float64 = 0.0
#    prod::Float64 = 1.0
#    for i in 1:n
#        total = 0.0
#        for j in 1:n
#            total += abs(tmp[i,j])^2
#        end
#        prod *= total
#    end
#    prod
#end
function estimadorP(mat::Array{Complex{Float64},2}, n::Int64, filas::Array{Int64,1})
    tmp = mat[filas, 1:n]
    total::Float64 = 0.0
    prod::Float64 = 1.0
    for i in 1:n
        total = 0.0
        for j in 1:n
            total += abs(tmp[i,j])^2
        end
        prod *= total
    end
    prod    
end

function generadistribucion(mat::Array{Complex{Float64},2})
    _, n::Int64 = size(mat)
    @assert n == 9
    lista::Array{Array{Int,1},1} = Array{Int,1}[]
    proba::Array{Float64,1} = Float64[]
    for i in 1:9, j in 1:9, k in 1:9
        pp = abs( per( mat[[i,j,k], 1:3] ) )^2
        push!(proba, pp/(*(factorial.([i,j,k])...)))
        push!(lista, [i,j,k])
    end
    lista, proba
end

end # module
