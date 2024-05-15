plot_dipole_moment!(fig, ::Nothing; kwargs...) = fig

function plot_dipole_moment!(fig, t, Fv, Av, z, vz=nothing, az=nothing;
                             title=nothing, yaxisposition=:left, kwargs...)
    au2fs = auconvert(u"fs", 1)
    tplot = ustrip.(au2fs*t)

    ax1 = Axis(fig[1,1], xticklabelsvisible=false, ylabel=L"$E(t)$ [au]"; title=title,
               yaxisposition=yaxisposition)
    lines!(ax1, tplot, Fv)
    ax2 = Axis(fig[2,1], xticklabelsvisible=false, ylabel=L"$A(t)$ [au]",
               yaxisposition=yaxisposition)
    lines!(ax2, tplot, Av)

    ax3 = Axis(fig[3,1], xlabel=L"$t$ [fs]",
               yaxisposition=yaxisposition)
    lines!(ax3, tplot, z, label=L"$z$")
    isnothing(vz) || lines!(ax3, tplot, vz, label=L"$v_z$")
    isnothing(az) || lines!(ax3, tplot, az, label=L"$a_z$")
    axislegend(ax3)

    fig
end

function plot_dipole_moment!(fig, results; kwargs...)
    t = results.t

    z = real(results.z)
    vz = real(results.vz)
    az = real(results.az)

    plot_dipole_moment!(fig, t, results.Ez, results.vp_mag, z, vz, az; kwargs...)
end

plot_dipole_moment(args...; kwargs...) =
    plot_dipole_moment!(Figure(size=(900,800)), args...; kwargs...)

function ωaxis(ωkind, ωunit, Iₚ, ω₀, Uₚ)
    offset,englabel = if ωkind == :kinetic
        Iₚ,L"$W_k$"
    else
        0,L"\omega"
    end
    if ωunit == "q"
        ω -> ω/ω₀, L"Harmonic order of $\omega_0$", identity
    elseif ωunit == "Up"
        ω -> (ω .- offset)/Uₚ, L"%$(englabel) ($U_p$)", identity
    elseif 1ωunit isa Unitful.Energy
        ω -> auconvert(ωunit, 1)*(ω .- offset), L"%$(englabel) [%$(ωunit)]", identity
    elseif 1ωunit isa Unitful.Length
        ω -> (auconvert(ωunit, 1)*2π*austrip(1u"c"))./ω, L"$\lambda$ [%$(ωunit)]", log10
    else
        @error "ωunit must either be an energy (e.g. u\"hartree\", u\"eV\", &c.) or \"q\" (indicating harmonic orders)"
        identity, L"\omega", identity
    end
end

round_quantity(v::Quantity) = ElectricFields.si_round(v)
round_quantity(v) = format("{1:.4f}", v)

function plot_dipole_spectrum!(fig, t,
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

    x = ustrip.(ωf(ω))

    ax = Axis(fig[1,1]; xlabel=ωlabel,
              ylabel="Dipole spectrum [arb.u.]",
              xscale=ωax, yscale=log10, kwargs...)

    lines!(ax, x, nan_map.(abs.(Z));
           label=L"$|-\omega^2 Z|$")

    isnothing(vz) || lines!(ax, x, nan_map.(abs.(Vz)),
                            label=L"$|-\mathrm{i}\omega V_z|$")
    isnothing(az) || lines!(ax, x, nan_map.(abs.(Az)),
                            label=L"$|A_z|$")

    Iₚ = ωf(Iₚ)
    vlines!(ax, [ustrip(Iₚ)], label="Ionization threshold $(round_quantity(Iₚ))")
    cutoff = ωf(cutoff)
    vlines!(ax, [ustrip(cutoff)], label="HHG cut-off $(round_quantity(cutoff))")

    axislegend(ax)

    fig
end

plot_dipole_spectrum!(fig, ::Nothing; kwargs...) = fig
function plot_dipole_spectrum!(fig, results; kwargs...)
    ip = results.input_params
    δt = ip.dt*ip.detail_frequency
    ω₀ = ip.vp_param[1][2][1]

    plot_dipole_spectrum!(fig, results.t, δt, ω₀,
                          results.Iₚ, results.Uₚ, results.cutoff,
                          results.z, results.vz, results.az;
                          kwargs...)
end

plot_dipole_spectrum(args...; kwargs...) =
    plot_dipole_spectrum!(Figure(size=(700,600)), args...; kwargs...)

function wind_transf(t, x, w; maxsteps=500, verbosity=0, kwargs...)
    f = rfftfreq(length(t),1.0/(t[2]-t[1]))
    ecat(a,b) = cat(a,b,dims=ndims(x)+1)
    nt = length(t)
    s = max(1, floor(Int, nt/maxsteps))
    tsel = 1:s:nt
    verbosity > 0 && @info "Time" nt s tsel
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
gabor(t, x, σ; kwargs...) = wind_transf(t, x, fftshift(exp.(-(t .- mean(t)).^2/2σ^2)); kwargs...)

plot_time_frequency_analysis!(fig, ::Nothing, args...; kwargs...) = fig
function plot_time_frequency_analysis!(fig, t, v, F, Iₚ, Uₚ, cutoff;
                                       ωkind=:total, ωunit=u"eV", window_length=0.05,
                                       dynamic_range=2, maxsteps=500,
                                       kwargs...)
    T = austrip(period(F))
    ω₀ = photon_energy(F)

    w = hanning(length(t))

    v = real(v)
    tt,f,W = gabor(t, w.*v, window_length*T; maxsteps=maxsteps)

    au2fs = auconvert(u"fs", 1)
    tplot = au2fs*tt

    ω = 2π*f
    ωf,ωlabel,ωax = ωaxis(ωkind, ωunit, Iₚ, ω₀, Uₚ)

    y = ustrip.(ωf(ω))
    Z = log10.(abs.(W))
    Zmin,Zmax = extrema(Z)

    p = sortperm(y)
    p = p[findall(isfinite, y[p])]

    ax = Axis(fig[1,1]; yscale=ωax,
              xlabel=L"$t$ [fs]", ylabel=ωlabel, kwargs...)

    ph = heatmap!(ax, ustrip.(tplot), y[p], Z[p,:]',
                  colorrange=(Zmax-dynamic_range,Zmax),
                  colormap=:inferno)

    Iₚ = ωf(Iₚ)
    hlines!(ax, [ustrip(Iₚ)], label="Ionization threshold $(round_quantity(Iₚ))")
    cutoff = ωf(cutoff)
    hlines!(ax, [ustrip(cutoff)], label="HHG cut-off $(round_quantity(cutoff))")
    axislegend(ax)

    fig
end

plot_time_frequency_analysis!(fig, results, F; kwargs...) =
    plot_time_frequency_analysis!(fig, results.t, results.az, F,
                                  results.Iₚ, results.Uₚ, results.cutoff; kwargs...)

plot_time_frequency_analysis(args...; kwargs...) =
    plot_time_frequency_analysis!(Figure(size=(700,600)), args...; kwargs...)
