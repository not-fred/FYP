using FFTW
FFTW.set_num_threads(1)

function βN(dB;β=1,Nₛ=10^4)
    return ceil(Int, min(2*10^(dB/10β), Nₛ))
    if N == Nₛ
        @warn "Target noise of −$(dB)dB probably not met, capped at $N"
    end
end

function βnoise(
        dτ=1,dB=30;
        β=1,f₀=0,
        N=nothing, Nₛ=10^4, # maximum N in case β is small
        returnFreq=false
    )

    # If not provided, get the sample size
    # according to the target magnitude
    if N == nothing
        N = βN(dB,β=β,Nₛ=Nₛ)
    end

    # Get the fourier frequencies,
    # cut off at f₀
    # NOTE: Used `collect` here because
    # indexing doesn't work with frequency
    f = collect(rfftfreq(N,(N-1)/dτ))
    samples = length(f)

    f₀ = max(f[2], f₀)
    f[f .< f₀] .= f₀

    # Calculate S(f) ∝ f⁻ᵝ
    S = f.^-β

    # Generate normal distributions
    Sᵣ = randn(samples) .* .√S
    if N%2 == 0
        # For even data points,
        # frequency is always real
        Sᵢ = 0
    else
        # The “f=0” component should be real
        Sᵢ = randn(samples) .* .√S
        Sᵢ[1] = 0
    end

    # Inverse Fourier transform
    # to get back the amplitudes
    y = irfft(Sᵣ .+ im*Sᵢ, N)

    # Normalise so that √2×RMS(y) = 1
    y /= sqrt(2sum(y.^2)/N)

    return returnFreq ? (f,y) : y
end
