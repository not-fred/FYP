using LinearAlgebra
include("../lib/generateBetaNoise.jl")

if !@isdefined ħ
    const ħ = 0.658211951E-3 # ħ ÷ 1meV, units of 1/ω
    const ∅ = nothing
    # Euler pulse sequences
    const ∠HZH = [
         2    ∅    ∅ ; # I
         ∅  -1/2   ∅ ; # R(x̂,⋅)
         ∅   1/2   ∅ ;
         ∅    1    ∅ ;
       -1/2 -1/2  1/2; # R(ŷ,⋅)
        1/2 -1/2 -1/2;
        1/2   1  -1/2;
       -1/2   ∅    ∅ ; # R(ẑ,⋅)
        1/2   ∅    ∅ ;
         1    ∅    ∅ ;
        1/2  1/2  1/2; # R(x̂+ẑ,π) = H
       -1/2  1/2 -1/2; # R(-x̂+ẑ,π) = H₋
        3/4   1  -3/4; # R(x̂+ŷ,⋅)
       -3/4   1   3/4;
         1  -1/2   ∅ ; # R(ŷ+ẑ,⋅)
         ∅  -1/2   1 ;
         ∅   1/2  1/2; # R(x̂+ŷ+ẑ,⋅)
        1/2  1/2   1 ;
        1/2 -1/2   1 ; # R(-x̂-ŷ+ẑ,⋅)
         ∅  -1/2  1/2;
        1/2  1/2   ∅ ; # R(x̂-ŷ+ẑ,⋅)
         1   1/2  1/2;
        1/2 -1/2   ∅ ; # R(-x̂+ŷ+ẑ,⋅)
         1  -1/2  1/2
    ]
    # Improved pulse sequences
    const ∠ZHH₋ = [
         2   ∅   ∅ ; # I
         0  3/2  0 ; # R(x̂,⋅)
         0  1/2  0 ;
         0   1   0 ;
         0   1   ∅ ; # R(ŷ,⋅)
         1   0   ∅ ;
        -1   0   ∅ ;
        3/2  ∅   ∅ ; # R(ẑ,⋅)
        1/2  ∅   ∅ ;
         1   ∅   ∅ ;
         0   ∅   ∅ ; # R(x̂+ẑ,π) = H
        -1   ∅   ∅ ; # R(-x̂+ẑ,π) = H₋
        1/2  0  -1 ; # R(x̂+ŷ,⋅)
        -1   0  1/2;
        -1  1/2  0 ; # R(ŷ+ẑ,⋅)
         0  1/2 -1 ;
        3/2  0   ∅ ; # R(x̂+ŷ+ẑ,⋅)
         0  1/2  ∅ ;
        -1  1/2  ∅ ; # R(-x̂-ŷ+ẑ,⋅)
        3/2 -1   ∅ ;
         0  3/2  ∅ ; # R(x̂-ŷ+ẑ,⋅)
        1/2  0   ∅ ;
        -1  3/2  ∅ ; # R(-x̂+ŷ+ẑ,⋅)
        1/2 -1   ∅ ;
    ]
    const cτ = denominator(rationalize(sqrt(2/3),tol=1E-4))
    const Nₛ = ceil(Int,10^4/cτ) * cτ
end

# Choice for discretising time for fourier transform
discreteτ(τ;p...)   = Int.(round.(τ / ħ / π / abs(p[:U]-p[:V]) * p[:t]^2 * 32 * cτ))
discreteτ⁻¹(τ;p...) = τ/cτ/32/p[:t]^2*abs(p[:U]-p[:V])*π*ħ

# Takes a series of Clifford indices, returns the pulse
# sequence for the entire gate sequence
function tqdpulseseq(S::AbstractArray{<:Int},τₖ=0;p...)
    tₗ = zeros(1)    # left tunneling parameter
    tᵣ = zeros(1)    # right tunneling parameter
    τ  = Float64[τₖ] # time
    τI = [1]         # index of gate. τ[τI[1]] is where
                     #     the first gate starts,
                     #     τ[τI[2]] is where the
                     #     second gate starts, etc.
    for i ∈ S
        # Essentially just appending the
        # pulse sequence for each individual gate
        # and collating them
        tₗₖ,tᵣₖ,τₖ = tqdpulseC(i,τ[end];p...)
        append!(tₗ,tₗₖ)
        append!(tᵣ,tᵣₖ)
        append!(τ,τₖ)
        push!(τI,length(τ))
    end
    return tₗ,tᵣ,τ,τI
end

# Introduce noise into the pulse sequences
function tqdnoiseypulseseq(S::AbstractArray{<:Int},τ₀=0;p...)
    # Obtain the “clean” pulse sequence
    tₗₖ,tᵣₖ,τₖ,τI = tqdpulseseq(S,τ₀;p...)
    keysP = keys(p)

    # Discretise the time into equal steps
    dτ = diff(τₖ)
    ndτ = discreteτ(dτ;p...)
    nτ = cumsum(ndτ)
    pushfirst!(nτ,0)
    N = nτ[end]

    # If target dB not met with chosen discretisation,
    # increase the number of time steps by some factor
    # to meet it
    dB = :dB ∈ keysP ? p[:dB] : 30
    Nfactor = N < βN(dB) ? ceil(Int,βN(dB)/N) : 1

    # Readjust discretisation based on the factor
    N   *= Nfactor
    ndτ *= Nfactor
    nτ  *= Nfactor
    dτ  = discreteτ⁻¹(ndτ;p...)/Nfactor
    τₖ  = τₖ[1] .+ discreteτ⁻¹(nτ;p...)/Nfactor
    Δτ = (τₖ[end]-τₖ[1])/N

    β = :β ∈ keys(p) ? p[:β] : 1

    # Initialise arrays in which to keep
    # the parameters
    tₗ = Float64[tₗₖ[1]]
    tᵣ = Float64[tᵣₖ[1]]
    ε  = Float64[p[:ε]]
    εₘ = Float64[p[:εₘ]]
    τ  = τₖ[1] .+ discreteτ⁻¹(0:N;p...)/Nfactor

    δt  = :δt  ∈ keysP ? p[:δt]  : p[:t]
    δtₗ = :δtₗ ∈ keysP ? p[:δtₗ] : δt
    δtᵣ = :δtᵣ ∈ keysP ? p[:δtᵣ] : δt
    if :δε ∈ keysP
        δεεₘ = p[:δε]
    elseif :δεₘ ∈ keysP
        δεεₘ = p[:δεₘ]
    elseif p[:ε] ≈ 0
        δεεₘ = p[:εₘ]
    elseif p[:εₘ] ≈ 0
        δεεₘ = p[:ε]
    end
    δε  = :δε  ∈ keysP ? p[:δε]  : δεεₘ
    δεₘ = :δεₘ ∈ keysP ? p[:δεₘ] : δεεₘ

    δt  *= βnoise(τₖ[end]-τₖ[1],dB,N=N,β=β,Nₛ=Nₛ)
    δtₗ *= βnoise(τₖ[end]-τₖ[1],dB,N=N,β=β,Nₛ=Nₛ)
    δtᵣ *= βnoise(τₖ[end]-τₖ[1],dB,N=N,β=β,Nₛ=Nₛ)
    δε  *= βnoise(τₖ[end]-τₖ[1],dB,N=N,β=β,Nₛ=Nₛ)
    δεₘ *= βnoise(τₖ[end]-τₖ[1],dB,N=N,β=β,Nₛ=Nₛ)

    for k ∈ 2:length(nτ)
        arrI = 1+nτ[k-1]:nτ[k]
        append!(tₗ, tₗₖ[k] .+ p[:δ]*δtₗ[arrI]*tₗₖ[k]/p[:t])
        append!(tᵣ, tᵣₖ[k] .+ p[:δ]*δtᵣ[arrI]*tᵣₖ[k]/p[:t])
        append!(ε , p[:ε]  .+ p[:δ]*δε[arrI])
        append!(εₘ, p[:εₘ] .+ p[:δ]*δεₘ[arrI])
    end
    return tₗ,tᵣ,ε,εₘ,τ,(nτ[τI].+1)
end

# Return pulse sequence of a specific
# Clifford gate
function tqdpulseC(i::Int,τₖ::Real;p...)
    tₗ = zeros(0)
    tᵣ = zeros(0)
    τ  = zeros(0)
    if :improved ∈ keys(p)
        ∠  = ∠ZHH₋[i,:]
        for k ∈ 1:3
             ∠[k] == ∅ && continue
             if ∠[k] == -1
                 tₗₖ,tᵣₖ,τₖ = tqdpulseH₋(τₖ;p...)
             elseif ∠[k] == 0
                 tₗₖ,tᵣₖ,τₖ = tqdpulseH(τₖ;p...)
             else
                 tₗₖ,tᵣₖ,τₖ = tqdpulseZ(∠[k],τₖ;p...)
             end
             push!(tₗ,tₗₖ)
             push!(tᵣ,tᵣₖ)
             push!(τ,τₖ)
        end
    else
        ∠ = ∠HZH[i,:]
        if i == 11
            tₗₖ,tᵣₖ,τₖ = tqdpulseH(τₖ;p...)
            push!(tₗ,tₗₖ)
            push!(tᵣ,tᵣₖ)
            push!(τ,τₖ)
        elseif i == 12
            tₗₖ,tᵣₖ,τₖ = tqdpulseH₋(τₖ;p...)
            push!(tₗ,tₗₖ)
            push!(tᵣ,tᵣₖ)
            push!(τ,τₖ)
        else
            for k ∈ 1:3
                ∠[k] == nothing && continue
                if k == 2
                    tₗₖ,tᵣₖ,τₖ = tqdpulseH(τₖ;p...)
                    push!(tₗ,tₗₖ)
                    push!(tᵣ,tᵣₖ)
                    push!(τ,τₖ)
                end
                tₗₖ,tᵣₖ,τₖ = tqdpulseZ(∠[k],τₖ;p...)
                push!(tₗ,tₗₖ)
                push!(tᵣ,tᵣₖ)
                push!(τ,τₖ)
                if k == 2
                    tₗₖ,tᵣₖ,τₖ = tqdpulseH(τₖ;p...)
                    push!(tₗ,tₗₖ)
                    push!(tᵣ,tᵣₖ)
                    push!(τ,τₖ)
                end
            end
        end
    end
    return tₗ,tᵣ,τ
end
tqdpulseH(τ::Real;p...)  = tqdpulseabstract([1,0,1],1,τ;p...)
tqdpulseH₋(τ::Real;p...) = tqdpulseabstract([-1,0,1],1,τ;p...)
tqdpulseZ(θ::Real,τ::Real;p...) = tqdpulseabstract([0,0,1],θ,τ;p...)

# Return the tunneling parameters and time
# given the rotation desired (Equation 11.4)
function tqdpulseabstract(n::AbstractArray{<:Real},θ::Real,τ::Real;p...)
    n = normalize(vec(n))
    θ = (2ceil(θ/2)-θ)%2
    θ ≈ 0 && (θ = 2)
    g₊ = 1/(p[:U]-2p[:V]-p[:εₘ]+p[:ε]) + 1/(p[:U]+p[:εₘ]-p[:ε])
    g₋ = 1/(p[:U]-2p[:V]-p[:εₘ]-p[:ε]) + 1/(p[:U]+p[:εₘ]+p[:ε])

    N = p[:t]^2*(g₋+g₊)
    τ += ħ/2N * θ*π
    tₗ = √( N/2g₊ * (n[3] + n[1]/√3) )
    tᵣ = √( N/2g₋ * (n[3] - n[1]/√3) )
    return tₗ,tᵣ,τ
end

# Returns the AEON TQD Hamiltonian given the
# control parameters. (Equation 11.1)
function tqdH(tₗ::Real,tᵣ::Real,ε::Real,εₘ::Real,τ::Real;p...)
    g₊ = 1 ./ (p[:U]-2p[:V].-εₘ+ε) + 1 ./ (p[:U].+εₘ-ε)
    g₋ = 1 ./ (p[:U]-2p[:V].-εₘ-ε) + 1 ./ (p[:U].+εₘ+ε)
    n₁ = √3( tₗ.^2 .* g₊ - tᵣ.^2 .* g₋ )
    n₃ =   ( tₗ.^2 .* g₊ + tᵣ.^2 .* g₋ )
    N  = √(n₁^2 + n₃^2)
    # Using exp[–ixn̂⋅σ⃗] = cos(x) I – i sin(x) n̂⋅σ⃗
    if N ≈ 0
        return I
    else
        return cos(τ*N/ħ)I + √2im*sin(τ*N/ħ)*σ([n₁/N,0,n₃/N])
    end
end
