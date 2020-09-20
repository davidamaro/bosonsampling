module BosonSampling

# IMPORTS
using Immanants
import RandomMatrices: Haar
import StatsBase: sample, Weights
import IterTools: subsets
# EXPORTS
export bosonsampling, per, Haar

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

end # module
