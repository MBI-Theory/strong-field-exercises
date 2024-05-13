plot_dipole_moment(::Nothing) = nothing

function plot_dipole_moment(t, Fv, Av, z, vz=nothing, az=nothing; kwargs...)
    au2fs = auconvert(u"fs", 1)
    tplot = au2fs*t

    pE = plot(tplot, Fv, ylabel=L"$E(t)$ [au]")
    pA = plot(tplot, Av, ylabel=L"$A(t)$ [au]")

    pz = plot(tplot, z, label=L"$z$")
    isnothing(vz) || plot!(pz, tplot, vz, label=L"$v_z$")
    isnothing(az) || plot!(pz, tplot, az, label=L"$a_z$")

    plot(pE, pA, pz;
         xlabel=L"$t$ [fs]", layout=@layout([a;b;c]),
         size=(700,800), kwargs...)
end

function plot_dipole_moment(results; kwargs...)
    t = results.t

    z = real(results.z)
    vz = real(results.vz)
    az = real(results.az)

    plot_dipole_moment(t, results.Ez, results.vp_mag, z, vz, az; kwargs...)
end

function ωaxis(ωkind, ωunit, Iₚ, ω₀, Uₚ)
    offset,englabel = if ωkind == :kinetic
        Iₚ,L"$W_k$"
    else
        0,L"\omega"
    end
    if ωunit == "q"
        ω -> ω/ω₀, L"Harmonic order of $\omega_0$", :identity
    elseif ωunit == "Up"
        ω -> (ω .- offset)/Uₚ, englabel*L" ($U_p$)", :identity
    elseif 1ωunit isa Unitful.Energy
        ω -> auconvert(ωunit, 1)*(ω .- offset), englabel, :identity
    elseif 1ωunit isa Unitful.Length
        ω -> (auconvert(ωunit, 1)*2π*austrip(1u"c"))./ω, L"\lambda", :log10
    else
        @error "ωunit must either be an energy (e.g. u\"hartree\", u\"eV\", &c.) or \"q\" (indicating harmonic orders)"
        identity, L"\omega", :identity
    end
end

round_quantity(v::Quantity) = ElectricFields.si_round(v)
round_quantity(v) = format("{1:.4f}", v)

function plot_dipole_spectrum(t,
                              δt, ω₀, Iₚ, Uₚ, cutoff,
                              z, vz=nothing, az=nothing;
                              ωkind=:total, ωunit=u"eV", window=hanning,
                              kwargs...)

    ω = 2π*rfftfreq(length(t), 1/δt)

    w = window(length(t))

    Z = -ω.^2 .* rfft(w .* real(z))
    Vz = isnothing(vz) ? nothing : im*ω .* rfft(w .* real(vz))
    Az = isnothing(az) ? nothing : rfft(w .* real(az))

    # Hack to please Plots.jl that refuses to deal with non-positive
    # values for logarithmic axes.
    nan_map(v) = v ≤ 0 ? NaN : v

    ωf,ωlabel,ωax = ωaxis(ωkind, ωunit, Iₚ, ω₀, Uₚ)

    pz = plot(ωf(ω), nan_map.(abs.(Z));
              xlabel=ωlabel,
              ylabel="Dipole spectrum [arb.u.]",
              label=L"$|-\omega^2 Z|$",
              xaxis=ωax, yaxis=:log10,
              size=(700,600),
              kwargs...)

    isnothing(vz) || plot!(pz, ωf(ω), nan_map.(abs.(Vz)),
                           label=L"$|-\mathrm{i}\omega V_z|$")
    isnothing(az) || plot!(pz, ωf(ω), nan_map.(abs.(Az)),
                           label=L"$|A_z|$")

    Iₚ = ωf(Iₚ)
    vline!(pz, ([Iₚ]), label="Ionization threshold $(round_quantity(Iₚ))")
    cutoff = ωf(cutoff)
    vline!(pz, ([cutoff]), label="HHG cut-off $(round_quantity(cutoff))")
    pz
end

plot_dipole_spectrum(::Nothing; kwargs...) = nothing
function plot_dipole_spectrum(results; kwargs...)
    ip = results.input_params
    δt = ip.dt*ip.detail_frequency
    ω₀ = ip.vp_param[1][2][1]

    plot_dipole_spectrum(results.t, δt, ω₀,
                         results.Iₚ, results.Uₚ, results.cutoff,
                         results.z, results.vz, results.az;
                         kwargs...)
end

function wind_transf(t, x, w)
    f = rfftfreq(length(t),1.0/(t[2]-t[1]))
    ecat(a,b) = cat(a,b,dims=ndims(x)+1)
    nt = length(t)
    s = max(1, floor(Int, nt/1000))
    tsel = 1:s:nt
    Cₓ = zeros(length(t), length(tsel))
    @withprogress name="Windowed FFT" begin
        for (ii,i) in enumerate(tsel)
            Cₓ[:,ii] .= x .* circshift(w,i)
            @logprogress ii/length(tsel)
        end
    end
    Wₓ = rfft(Cₓ, 1)
    t[tsel], f,Wₓ
end
gabor(t, x, σ) = wind_transf(t, x, fftshift(exp.(-(t .- mean(t)).^2/2σ^2)))

plot_time_frequency_analysis(::Nothing, args...; kwargs...) = nothing
function plot_time_frequency_analysis(t, v, F, Iₚ, Uₚ, cutoff;
                                      ωkind=:total, ωunit=u"eV", window_length=0.05, dynamic_range=2)
    T = austrip(period(F))
    ω₀ = photon_energy(F)

    w = hanning(length(t))

    v = real(v)
    tt,f,W = gabor(t, w.*v, window_length*T)

    au2fs = auconvert(u"fs", 1)
    tplot = au2fs*tt

    ω = 2π*f
    ωf,ωlabel,ωax = ωaxis(ωkind, ωunit, Iₚ, ω₀, Uₚ)

    y = ωf(ω)
    Z = log10.(abs.(W))
    Zmin,Zmax = extrema(Z)

    p = sortperm(y)
    p = p[findall(isfinite, y[p])]

    ph = heatmap(ustrip.(tplot), y[p], Z[p,:],
                 yaxis=ωax,
                 xlabel=L"$t$ [fs]",
                 ylabel=ωlabel,
                 clim=(Zmax-dynamic_range,Zmax),
                 title=L"\log_{10} |G(t,W_k)|",
                 size=(700,600))

    Iₚ = ωf(Iₚ)
    hline!(ph, ([Iₚ]), label="Ionization threshold $(round_quantity(Iₚ))")
    cutoff = ωf(cutoff)
    hline!(ph, ([cutoff]), label="HHG cut-off $(round_quantity(cutoff))")
end

plot_time_frequency_analysis(results, F; kwargs...) =
    plot_time_frequency_analysis(results.t, results.az, F,
                                 results.Iₚ, results.Uₚ, results.cutoff; kwargs...)
