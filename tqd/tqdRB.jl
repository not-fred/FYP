import Statistics.mean, Statistics.std

include("tqd.jl")
include("../lib/superoperators.jl")
include("../lib/cliffordGates.jl")
include("../lib/RB.jl")

function tqdRB(M,K=500;ρ=[1,0]*[1,0]',Q=copy(ρ),params...)
    #M = M₀:ΔM:M
    F = zeros(length(M), K)
    τ = zeros(length(M), K)
    r̄ = zeros(sum((M.+1)*K))

    𝔾 = []
    𝕘 = []
    for m ∈ 1:length(M)
        for r ∈ 1:K
            S  = rand(1:24,M[m])
            Sᵢ = reduce(getCproduct, reverse(S))
            push!(S, getCinverse(Sᵢ))

            tₗ,tᵣ,ε,εₘ,τₖ,τI = tqdnoiseypulseseq(S;params...)

            G = I
            for k ∈ 1:M[m]+1

                C̃ = I

                for j ∈ τI[k]+1:τI[k+1]
                    R = tqdH(tₗ[j],tᵣ[j],ε[j],εₘ[j],τₖ[j]-τₖ[j-1];params...)
                    C̃ = R*C̃
                end

                C = rotate(CliffordAngles[S[k]]...)

                push!(𝔾,C)
                push!(𝕘,C̃)

                r̄[
                    K*sum(M[1:m-1].+1) +
                    (r-1)*(M[m]+1) +
                    k
                ] = infidelity(C,C̃)

                G = C̃*G
            end
            F[m,r] = tr(Q*G*ρ*G')
            τ[m,r] = τₖ[end]
        end
    end
    p̄₀,A₀,B₀ = eigenRB(𝔾,𝕘,ρ=supervec(ρ),Q=supervec(Q))
    p̄ = [p̄₀,A₀,B₀]
    return F,τ,p̄,mean(r̄),std(r̄)
end
