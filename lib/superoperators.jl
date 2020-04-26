using LinearAlgebra

if !@isdefined σₙ
    const σₙ = [
        [1 0; 0 1],
        [0 1; 1 0],
        [0 -im; im 0],
        [1 0; 0 -1]
    ]/√2
    σ(n::Int) = σₙ[n+1]
    σ(n::AbstractArray{<:Real}) = sum(n[i]*σₙ[i+(length(n)==3)] for i ∈ 1:length(n))
    ⊗(A::AbstractArray, B::AbstractArray) = kron(A, B)
end

function superop(𝓔,u=1)
    if 𝓔 isa Tuple || size(𝓔) == (3,3)
        if 𝓔 isa Tuple
            g,G = 𝓔
        else
            g = zeros(3)
            G = 𝓔
        end
        return real.([[u;g] [0;G]])
    elseif !isempty(methods(𝓔)) && size(𝓔(σₙ[1])) == (2,2) # Check if 𝓔 is a linear map
        𝓖 = [tr(σₙ[i]'*𝓔(σₙ[j])) for i ∈ 1:4,j ∈ 1:4]
        𝓖[1] = u
        return real.(𝓖)
    elseif size(𝓔) == (2,2) # If 𝓔 is an operator
        𝓖 = [tr(σₙ[i]'*𝓔*σₙ[j]*𝓔') for i ∈ 1:4, j ∈ 1:4]
        𝓖[1] = u
        return real.(𝓖)
    elseif length(𝓔) > 1 && all(size(Eₖ) == (2,2) for Eₖ ∈ 𝓔)
        return sum(superop(Eₖ,u) for Eₖ ∈ 𝓔)
    else
        # Error
    end
end

function supervec(ρ)
    return [tr(σₙ[i]'*ρ) for i ∈ 1:4]
end

function rotate3D(m::AbstractArray{<:Real},ϑ::Real)
    m = normalize(vec(m))
    mₓ = zeros(3,3)
    for i ∈ 1:3
        mₓ[:,i] = m × (1:3 .== i)
    end
    return cos(ϑ*π)I + sin(ϑ*π)mₓ + (1-cos(ϑ*π))m*m'
end

function rotate(m::AbstractArray{<:Real}, ϑ::Real)
    m = normalize(vec(m))
    return cos.(ϑ*π/2)*I - sin.(ϑ*π/2)im*σ(m)*√2
end

function infidelity(G,G̃)

    if G isa UniformScaling
        G *= I(3)
    elseif G isa Tuple
        G = G[2]
    elseif size(G) == (4,4)
        G = G[2:end,2:end]
    elseif size(G) == (2,2)
        G = superop(G)[2:end,2:end]
    end

    if G̃ isa UniformScaling
        G̃ *= I(3)
    elseif G̃ isa Tuple
        G̃ = G̃[2]
    elseif size(G̃) == (4,4)
        G̃ = G̃[2:end,2:end]
    elseif size(G̃) == (2,2)
        G̃ = superop(G̃)[2:end,2:end]
    end

    GG̃ = G'G̃
    GG̃ isa UniformScaling && (GG̃ *= I(3))

    return 1/2-real(tr(GG̃))/6
end
