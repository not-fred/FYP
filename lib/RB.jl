using Statistics, LsqFit, GLM
import DataFrames.DataFrame

function standardRB(ℂ,𝕔,M,K=500;n=[0,0,1],q=copy(n),Δ=false)
    F  = zeros(length(M))
    σF = zeros(length(M))
    r̄ₖ = zeros(sum((M.+1)*K))
    𝔾 = Any[]
    𝕘 = Any[]
    Δ𝕘 = Δ != false
    if (Δ𝕘)
        ϵ = zeros(length(M))
        σϵ = zeros(length(M))
    else
        ϵ = false
        σϵ = false
    end
    for m ∈ 1:length(M)
        Fₘ = zeros(K)
        if (Δ𝕘 != false)
            ϵₘ = zeros(K)
        end
        for k ∈ 1:K
            S  = rand(1:24,M[m])
            Sᵢ = reduce(getCproduct, reverse(S))
            push!(S, getCinverse(Sᵢ))
            nₖ = copy(n)
            Δ𝕘 && (nₑ = [1;copy(n)])
            for i ∈ 1:M[m]+1
                nₖ = 𝕔[S[i]][1] + 𝕔[S[i]][2]nₖ
                Δ𝕘 && (nₑ = Δ[S[i]]nₑ)
                r̄ₖ[
                    K*sum(M[1:m-1].+1) +
                    (k-1)*(M[m]+1) +
                    i
                ] = infidelity(ℂ[S[i]],𝕔[S[i]])
                push!(𝔾,ℂ[S[i]])
                push!(𝕘,𝕔[S[i]])
            end
            Fₘ[k] = (1 + q⋅nₖ)/2
            Δ𝕘 && (ϵₘ[k] = abs(([1;q]⋅nₑ)/2))
        end
        F[m]  = mean(Fₘ)
        σF[m] = std(Fₘ)/√(K+1)
        if (Δ𝕘)
            ϵ[m]  = mean(ϵₘ)
            σϵ[m] = std(ϵₘ)/√(K+1)
        end
    end
    return F,σF,mean(r̄ₖ),std(r̄ₖ),𝔾,𝕘,ϵ,σϵ
end

modelF(m,p) = p[2]*p[1].^m .+ p[3]
jacobF(m,p) = [ p[2] .* m .* p[1].^(m.-1) p[1].^m  ones(length(m)) ]
function fitRB(M,F,σF)
    p₀ = [1.,.5,.5]
    pₗ = [0.,0.,0.]
    pᵤ = [1.,1.,1.]
    p = fill(-1,3)
    σp = fill(-1,3)
    try
        fitF = curve_fit(modelF,jacobF,M,F,σF.^(-2),p₀,lower=pₗ,upper=pᵤ)
        !fitF.converged && fail()
        p = coef(fitF)
        σp = margin_error(fitF)
    catch
        try
            fitF = curve_fit(modelF,jacobF,M,F,p₀,lower=pₗ,upper=pᵤ)
            p = coef(fitF)
            σp = margin_error(fitF)
        catch
            F₁ = copy(F)
            F₁[F .<= .5] .= (F₁[F .<= .5] + σF[F .<= .5] .+ .5)/2
            wts = copy(σF)
            wts[F .<= .5] .=  F₁[F .<= .5] .- .5
            wts ./= (F₁.-.5)
            fitF = lm(@formula(y~x),DataFrame(
                        y=log.(F₁[F₁ .> .5].-.5),
                        x=M[F₁ .> .5]),
                    wts=wts[F₁ .> .5].^(-2))
            p = exp.(coef(fitF))
            σp = p .* 2stderror(fitF)
            p = [p[end:-1:1]...,.5]
            σp = [σp[end:-1:1]...,0.]
        end
    end
    return p,σp
end

function eigenRB(𝔾,𝕘;ρ=[1,0,0,1]/√2,Q=copy(ρ))

    # Averaging 𝓖ᵤ ⊗ 𝒢 and 𝒢† ⊗ 𝓖ᵤ†
    𝓖𝒢 = sum(superop(𝔾[k],0) ⊗ superop(𝕘[k])    for k ∈ 1:length(𝔾))/length(𝔾)
    𝒢𝓖 = sum(superop(𝕘[k])'  ⊗ superop(𝔾[k],0)' for k ∈ 1:length(𝔾))/length(𝔾)

    # Averaging 𝒢
    𝒢 = sum(superop(𝕘[k]) for k ∈ 1:length(𝔾))/length(𝔾)

    # Finding L₀ and R₀
    _, L₀ₖ = eigen(𝒢)
    L₀ = L₀ₖ[:,end]/L₀ₖ[1,end]
    _,R₀ₖ = eigen(𝒢')
    R₀ = R₀ₖ[:,end]'/R₀ₖ[1,end]

    # Finding 𝓛ᵤ, 𝓡ᵤ
    p, v𝓛ᵤ = eigen(𝓖𝒢)
    _, v𝓡ᵤ = eigen(𝒢𝓖)
    p = p[end] # Largest value of p
    𝓛ᵤ = √3reshape(v𝓛ᵤ[:,end],4,4)
    𝓡ᵤ = √3reshape(v𝓡ᵤ[:,end],4,4)

    # Making sure left-hand-rule is preserved
    if det(real(𝓛ᵤ[2:end,2:end])) < 0
        𝓛ᵤ *= -1
    end
    if det(real(𝓡ᵤ[2:end,2:end])) < 0
        𝓡ᵤ *= -1
    end

    A = Q⋅(𝓛ᵤ*𝓡ᵤ*ρ)
    B = (Q⋅L₀)*(R₀⋅ρ)
    𝓡 = 𝓡ᵤ + [1,0,0,0]*R₀
    𝓛 = 𝓛ᵤ + L₀*[1 0 0 0]

    # Δ = 𝒢 - 𝓛𝓖𝓡
    Δ = [ real.(superop(𝕘[k]) - 𝓛*superop(𝔾[k])*𝓡) for k ∈ 1:length(𝔾) ]
    return real.((p,A,B))..., Δ
end
