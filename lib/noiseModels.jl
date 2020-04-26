using LinearAlgebra
!(@isdefined σₙ) && include("superoperators.jl")

function boundinfidelity(r, rᵤ=2/3; type="", rᵤtext="≤ ⅔")
    if r < 0 || r > rᵤ
        type != "" && (type = " for $type")
        throw(DomainError(r, "The gate infidelity$type must be bounded by 0 ≤ r $rᵤtext"))
    end
end

function randomΛ(r,Nₖ=5;maxLoop=1000)
    boundinfidelity(
        r, 1/2-eps(),
        type="random errors", rᵤtext="< ½"
    )
    # δ is just some parameterisation that signifies
    # the magnitude of the random matrices generated
    # found numerically, and meant to overshoot
    # the target infidelity (and scaled back to the
    # target later on)
    δ = 10atanh(2*r)/7
    r̄ = 0
    λ,Λ = zeros(3),I

    # Function to generate random real 2×2 matrices
    A() = (2rand(2,2) .- 1)

    # Keep generating Λ until its infidelity exceeds
    # the target infidelity
    while r̄ < r && maxLoop > 0
        # Generate random process matrices Lₖ,
        # then calculate B = Σₖ Lₖ' Lₖ to first enforce
        # B ≲ 1 (See 4.3.5 Quantum Processes), slightly
        # less than 1 so a one more Lₖ can be
        # included to make them sum to identity.
        Lₖ = [I + δ*(A()+im*A()) for k ∈ 1:Nₖ-1]
        B = sum(Lₖ[k]'Lₖ[k] for k ∈ 1:Nₖ-1)
        Lₖ = Lₖ./√eigmax(B)*(1-1e-9)

        # Now find the last C needed to
        # turn Σₖ [ Lₖ' Lₖ ] + C = I.
        C = I - sum(Lₖ[k]'Lₖ[k] for k ∈ 1:Nₖ-1)
        # Then, use Cholesky decomposition
        # to find the upper triangular matrix
        # Lₗ that satisfies C = Lₗ' Lₗ,
        # and add a random rotation
        # so it wouldn't be triangular
        # (just looks nicer and more random).
        # It's alright to add a rotation because
        #    Lₗ → U Lₗ ⟹ (U Lₗ)' (U Lₗ)
        #                = Lₗ' U' U Lₗ
        #                = Lₗ' Lₗ
        #                = C
        U = rotate(2rand(3).-1,2δ)
        push!(Lₖ,U*cholesky(C).U)
        Λ = superop(Lₖ)

        # Convert into the e⃗,E form (Equation 5.18,5.19)
        λ,Λ = Λ[2:end,1],Λ[2:end,2:end]
        r̄ = infidelity(I,(λ,Λ))
        maxLoop -= 1
    end
    if r̄ > r
        # If the actual infidelity is above the target infidelity,
        # proportionally include the identity channel
        # to make the actual infidelity exactly
        # the target infidelity
        λ,Λ = r/r̄*λ, r/r̄*Λ + (1-r/r̄)*I
    end
    return λ,Λ
end


# ================================================================= #
#   See Table 6.1 for bitflip, depolarise, amplitudedamp, random𝓤   #
# ================================================================= #

function bitflip(r)
    boundinfidelity(r)
    p = 1-3r/2
    return zeros(3), p*I + (1-p)*rotate3D([1,0,0],1)
end

function depolarise(r)
    boundinfidelity(r)
    return zeros(3), (1-2r)*I
end

function amplitudedamp(r)
    boundinfidelity(r, 1/2, type="amplitude damping", rᵤtext="≤ ½")
    p = (2*√(1-3r/2)-1)^2
    return [0,0,1-p], diagm([√p,√p,p])
end

function random𝓤(r)
    boundinfidelity(r)
    δ = acos(1-3r)/π
    return zeros(3),rotate3D(2rand(3).-1,δ)
end


function addnoise(𝔾,r₀,dep=0,r=r₀;genΛ=randomΛ,genΛₖ=genΛ)
    # Writing the dirty gate as
    # 𝒢ᵢ = (𝓛 + ℒᵢ)𝓖ᵢ(𝓡 + ℛᵢ),
    # where 𝓛,𝓡 is the independent noise and
    # ℒᵢ,ℛᵢ is the dependent noise,
    # r₀ is the total infidelity of 𝓛 and 𝓡,
    # r  is the total infidelity of ℒᵢ and ℛᵢ.
    #
    # Ranges can be given, in which case
    # a random infidelity between the range
    # will be used.
    r₀ = extrema(r₀)
    r₀ = (r₀[2]-r₀[1])*rand() + r₀[1]
    r = extrema(r)

    𝓵,𝓛 = zeros(3),I
    if genΛ isa Tuple
        pₗ = rand()
        𝓵,𝓛 = genΛ[1](pₗ*r₀)
        𝓻,𝓡 = genΛ[2]((1-pₗ)*r₀)
    else
        𝓻,𝓡 = genΛ(r₀)
    end

    𝕘 = [(zeros(3),g) for g ∈ 𝔾]
    for k ∈ 1:length(𝕘)
        if dep == 0
            𝓡ₖ = 𝓡
            𝓛ₖ = 𝓛
            𝓻ₖ = 𝓻
            𝓵ₖ = 𝓵
        else
            rₖ = (r[2]-r[1])*rand() + r[1]
            𝓵ₖ,𝓛ₖ = zeros(3),I
            if genΛₖ isa Tuple
                δₗ = rand()
                𝓵ₖ,𝓛ₖ = genΛₖ[1](δₗ*rₖ)
                𝓻ₖ,𝓡ₖ = genΛₖ[2]((1-δₗ)*rₖ)
            else
                𝓻ₖ,𝓡ₖ = genΛₖ(rₖ)
            end
            # (Equation 9.2)
            𝓡ₖ = (1-dep)*𝓡 + dep*𝓡ₖ
            𝓛ₖ = (1-dep)*𝓛 + dep*𝓛ₖ
            𝓵ₖ = (1-dep)*𝓵 + dep*𝓵ₖ
            𝓻ₖ = (1-dep)*𝓻 + dep*𝓻ₖ
        end
        𝕘[k] = 𝓵ₖ+𝓻ₖ, 𝓛ₖ*𝕘[k][2]*𝓡ₖ
    end
    return 𝕘
end
