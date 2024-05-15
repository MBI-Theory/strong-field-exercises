### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ ce907c44-0391-11ef-123a-1d651bb74fcc
begin
    using Pkg
    Pkg.activate(".")

    using PlutoUI
    import PlutoUI: combine
    using HypertextLiteral
    using ProgressLogging

    using CairoMakie
    using LaTeXStrings
    using DomainColoring

    using Statistics
    using FFTW
    using AbstractFFTs
    using DSP

    using Unitful
    using UnitfulAtomic

    using ElectricFields

    using ComplexRegions
    using RationalFunctionApproximation

    using Revise

    using NotebookHelper    

    using StrongFieldApproximation
end

# ╔═╡ 8e113bbb-6d30-481d-884b-11ecc81ae6e8
@bind inputs let
    is = InputSection("Calculation inputs", [
        InputSection("Electric field", [
            Input("I₀", NumberIntervalField(0.0..1000, default=300.0), "Intensity", unit="TW/cm²"),
            Input("λorω", Select([:λ => "λ", :ω => "ħω"], default=:λ), "Carrier"),
            Input("λ", NumberIntervalField(1.0..1000, default=800.0), "Wavelength", unit="nm"),
            Input("ω", NumberIntervalField(0.01..10, default=1.0), "``\\hbar\\omega``", unit="Ha"),
            Input("cycles", NumberIntervalField(1.0..30, default=2.0), "Duration", unit="cycles"),
            Input("envelope", Select([:gaussian => "Gaussian", :tophat => "Top-hat pulse"],
                                     default=:tophat), "Pulse envelope")
        ]),
    ])

    confirm(format_input(is))
end

# ╔═╡ 99655ba4-b323-41d3-ad36-484db298ebe4
F = get_electric_field(inputs; gauss_kind=:gauss, ramp=0.0)

# ╔═╡ bab736ac-5aa5-478b-88b7-c20ee328f122
ndt = 100 # Steps per cycle

# ╔═╡ 95624be6-9f58-4a0e-953b-545d1a9caaa0
Iₚ = 0.5 # Hydrogen

# ╔═╡ 29554dea-2b64-4976-b850-48b5aaeca0b8
d = Base.Fix2(d_hyd, hyd_α(Iₚ));

# ╔═╡ 4c28f702-f278-4337-bd5e-6f8a16e723d2
Uₚ = ponderomotive_potential(F) |> u"eV"

# ╔═╡ 180fd636-4acd-4d74-af7a-c520e83a407c
t = timeaxis(F, ndt);

# ╔═╡ 43789066-1806-400e-a94d-cb66c11eff78
tim = 1austrip(period(F))*range(-1,stop=1,length=401);

# ╔═╡ f5551674-12ab-4386-925e-de85043450f6
dlim = (t[1],t[end],tim[1],tim[end]);

# ╔═╡ f8c36bf3-6ef9-4039-8acb-09fb3d6ba93e
begin
    pixels = let
        np = 600
        ceil(Int, np*(dlim[2]-dlim[1])/100), ceil(Int, np*(dlim[4]-dlim[3])/100)
    end
end;

# ╔═╡ d9d59e03-19ff-4669-8e83-6ddf73a5982b
tc = t' .+ im*tim;

# ╔═╡ 75169dc6-99a0-44e6-a8cb-520aa6c633b2
au2fs = auconvert(u"fs", 1);

# ╔═╡ fbee840d-14ed-4fd3-86d4-66c2c53e5ebd
tplot = au2fs*t;

# ╔═╡ 59c9b684-d5a9-4d98-9706-539bc2d7ea90
timplot = tim*au2fs;

# ╔═╡ ebc718dc-d506-4168-8b04-ee696aa1e1ab
kx = 3range(-1, stop=1, length=301);

# ╔═╡ 8578f4a0-e50d-4789-a5e4-5a621e3f066c
lines(kx, imag(d.(kx)), axis=(xlabel=L"$k_x$ [Bohr$^{-1}$]", title=L"\Im[d(k_x)]"))

# ╔═╡ aa86939e-7add-4379-ae4f-29c51dbf7cfa
Fv = field_amplitude(F, t);

# ╔═╡ bab8dc3b-7755-4487-875a-f189ed40335a
Av = vector_potential(F, t);

# ╔═╡ 1da023e0-b0d6-48e2-88f0-97809af5ef03
Fvc = field_amplitude.(F, tc);

# ╔═╡ 4c3d86a2-08d3-4c4d-aea3-b3e2204279c8
Avc = vector_potential.(F, tc);

# ╔═╡ 978c2eeb-d9af-47bc-8340-e6f411e392f3
vp = StrongFieldApproximation.VolkovPhases(F, t);

# ╔═╡ 916feb89-509e-466e-ba7e-8a27b90a986d
let
    fig = Figure(size=(700,900))
    ax1 = Axis(fig[1,1], ylabel=L"F(t)", xticklabelsvisible=false)
    lines!(ax1, ustrip.(tplot), Fv)
    ax2 = Axis(fig[2,1], ylabel=L"A(t)", xticklabelsvisible=false)
    lines!(ax2, ustrip.(tplot), Av
    )

    ax3 = Axis(fig[3,1], xlabel=L"$t$ [fs]")
    lines!(ax3, ustrip.(tplot), vp.∫A.(t), label=L"\int_t^T\mathrm{d}\tau A(\tau)")
    lines!(ax3, ustrip.(tplot), vp.∫A².(t), label=L"\int_t^T\mathrm{d}\tau A^2(\tau)")
    axislegend(ax3, position=:rb)

    # save("sfa-ati-fields.pdf", fig)
    fig
end

# ╔═╡ ef12d365-2b25-42a5-8c29-0ce5ec1a115c
vpc = StrongFieldApproximation.VolkovPhases(F, t, tim);

# ╔═╡ b06bc47a-f5eb-4b0f-961c-fdb50618a0f9
let
    fig = Figure(size=(700,1000))
    ax1 = Axis(fig[1,1], title=L"A(t)",
               xticklabelsvisible=false, ylabel=L"$ℑ(t)$ [au]")
    domaincolor!(ax1, Base.Fix1(vector_potential, F), dlim, pixels=pixels, rasterize=true)

    ax2 = Axis(fig[2,1], title=L"\int_t^T\mathrm{d}\tau A(\tau)",
               xticklabelsvisible=false, ylabel=L"$ℑ(t)$ [au]")
    domaincolor!(ax2, z -> vpc.∫A(z), dlim, pixels=pixels, rasterize=true)

    ax3 = Axis(fig[3,1], title=L"\int_t^T\mathrm{d}\tau A^2(\tau)",
               xlabel=L"$\Re(t)$ [au]", ylabel=L"$ℑ(t)$ [au]")
    domaincolor!(ax3, z -> vpc.∫A²(z), dlim, pixels=pixels, rasterize=true)

    # save("sfa-ati-fields-complex.pdf", fig)
    fig
end

# ╔═╡ be1211d9-4fed-4e48-8e7a-13b27204cff0
let
    fig = Figure(size=(700,1400))

    minmax(Z) = (Zmax = maximum(abs, Z); (-Zmax,Zmax))
    pcan = kx .+ Av'

    ax1 = Axis(fig[1,1],
               title=L"k_x + A(t)",
               xticklabelsvisible=false, ylabel=L"$k$ [Bohr$^{-1}$]")
    Colorbar(fig[1,2],
             heatmap!(ax1, ustrip.(tplot), kx, pcan',
                      colorrange=minmax(pcan), colormap=:RdBu, rasterize=true))

    dcan = imag.(d.(pcan))
    ax2 = Axis(fig[2,1],
               title=L"\Im(d[k_x + A(t)])",
               xticklabelsvisible=false, ylabel=L"$k$ [Bohr$^{-1}$]")
    Colorbar(fig[2,2],
             heatmap!(ax2, ustrip.(tplot), kx, dcan',
                      colorrange=minmax(dcan), colormap=:RdBu, rasterize=true))

    Fdcan = Fv' .* d.(pcan)
    imFdcan = imag.(Fdcan)
    ax3 = Axis(fig[3,1],
               title=L"\Im(F_x(t) d[k_x + A(t)])",
               xticklabelsvisible=false, ylabel=L"$k$ [Bohr$^{-1}$]")
    Colorbar(fig[3,2],
             heatmap!(ax3, ustrip.(tplot), kx, imFdcan',
                      colorrange=minmax(imFdcan), colormap=:RdBu, rasterize=true))

    Λ = zeros(length(kx), length(t))
    for (i,k) in enumerate(kx)
        for (j,t) in enumerate(t)
            Λ[i,j] = StrongFieldApproximation.volkov_phase(k, vp, t)
        end
    end
    S = Λ .+ Iₚ*(t[end] .- t)'

    ax4 = Axis(fig[4,1],
               title=L"\cos[I_p(T-t) + \Lambda(k_x,T)-\Lambda(k_x,t)]",
               xticklabelsvisible=false, ylabel=L"$k$ [Bohr$^{-1}$]")
    heatmap!(ax4, ustrip.(tplot), kx, cos.(S)',
             colorrange=(-1,1), colormap=:RdBu, rasterize=true)
    ax5 = Axis(fig[5,1],
               title=L"\sin[I_p(T-t) + \Lambda(k_x,T)-\Lambda(k_x,t)]",
               xlabel=L"$t$ [fs]", ylabel=L"$k$ [Bohr$^{-1}$]")
    heatmap!(ax5, ustrip.(tplot), kx, sin.(S)',
             colorrange=(-1,1), colormap=:RdBu, rasterize=true)

    # save("sfa-ati-action.pdf", fig)
    fig
end

# ╔═╡ 8ed2ebf7-43ff-42c8-8299-99ca0bb94377
kkx = 1.0

# ╔═╡ 26465e4e-a6e6-42ed-b5a6-3da23421bb01
let
    fig = Figure(size=(800,700))

    minmax(Z) = (Zmax = maximum(abs, Z); (-Zmax,Zmax))
    pcan = kkx .+ Av
    ax1 = Axis(fig[1,1], xticklabelsvisible=false, ylabel=L"k_x + A(t_1)")
    pcanonical = lines!(ax1, ustrip.(tplot), pcan)

    dcan = imag.(d.(pcan))
    ax2 = Axis(fig[2,1], xticklabelsvisible=false,
               ylabel=L"\Im(d[k_x + A(t_1)])")
    lines!(ax2, ustrip.(tplot), dcan)

    Fdcan = imag.(Fv .* d.(pcan))
    ax3 = Axis(fig[3,1], xlabel=L"$t_1$ [fs]",
               ylabel=L"\Im(F_x(t_1) d[k_x + A(t_1)])")
    lines!(ax3, ustrip.(tplot), Fdcan)

    Λ = zeros(length(t))
    for (j,t) in enumerate(t)
        Λ[j] = StrongFieldApproximation.volkov_phase(kkx, vp, t)
    end
    S = Λ .+ Iₚ*(t[end] .- t)

    ax4 = Axis(fig[1,2], xticklabelsvisible=false,
               ylabel=L"S", yaxisposition=:right)
    lines!(ax4, ustrip.(tplot), S)

    c = cos.(S)
    s = sin.(S)
    ax5 = Axis(fig[2,2], xticklabelsvisible=false, ylabel=L"\exp(-\mathrm{i}S)", yaxisposition=:right)
    lines!(ax5, ustrip.(tplot), c,
           label=L"\cos(S)")
    lines!(ax5, ustrip.(tplot), s,
           label=L"\sin(S)")
    axislegend(ax5, position=:rb)

    ax6 = Axis(fig[3,2], xlabel=L"$t_1$ [fs]", ylabel=L"\exp(-\mathrm{i}S)(F_x(t_1) d[k_x + A(t_1)])", yaxisposition=:right)
    peSFd = lines!(ax6, ustrip.(tplot), c .* Fdcan,
                   label=L"\cos(S)\Im(F_x(t_1) d[k_x + A(t_1)])")
    lines!(ax6, ustrip.(tplot), s .* Fdcan,
           label=L"\sin(S)\Im(F_x(t_1) d[k_x + A(t_1)])")
    axislegend(ax6, position=:rb)

    # save("sfa-ati-action-one-momentum.pdf", fig)
    fig
end

# ╔═╡ 5cdb2f14-58a4-4b33-bedb-7a991c51d059
∂Sf(z,k) = -Iₚ - (k + vector_potential(F, z))^2 / 2;

# ╔═╡ ede51caa-e20b-4420-9193-158aec3e830d
get_r(k) = approximate(Base.Fix2(∂Sf, k),
                       rectangle(dlim[1]+im*max(dlim[3],0), dlim[2]+im*max(dlim[4],0)));

# ╔═╡ 053c9566-ca55-440d-9544-905944212774
begin
    saddle_points_kkx = let
        rf = get_r(kkx)
        # @info "Rational approximation" rf roots(rf) poles(rf)
        filter(z -> imag(z) > 0, roots(rf))
    end
end;

# ╔═╡ 92adcd20-eecd-4737-bf25-6dd76ce878c7
function plot_saddle_points!(args...)
    for ax in args
        li = [ax.xaxis.attributes.limits[]..., ax.yaxis.attributes.limits[]...]
        scatter!(ax, real(saddle_points_kkx), imag(saddle_points_kkx),
                 color=("white", 0), strokewidth=1, strokecolor="white",
                 markersize=7, marker=:diamond)
        limits!(ax, li...)
    end
end;

# ╔═╡ d6d87b41-55a2-4b19-a6b1-715beb770e32
let
    Sf = z -> StrongFieldApproximation.volkov_phase(kkx, vpc, z) + Iₚ*(t[end] - z)
    pcanf = z -> kkx + vector_potential(F, z)
    dcanf = z -> d(pcanf(z))

    fig = Figure(size=(800,1000))

    ax1 = Axis(fig[1,1], title=L"S",
               xticklabelsvisible=false, ylabel=L"$\Im[t_1]$ [jiffies]")
    domaincolor!(ax1, Sf, dlim, pixels=pixels, rasterize=true)
    ax2 = Axis(fig[2,1], title=L"\exp(-\mathrm{i}S)",
               xticklabelsvisible=false, ylabel=L"$\Im[t_1]$ [jiffies]")
    domaincolor!(ax2, z -> exp(-im*Sf(z)), dlim, pixels=pixels, abs=Inf, rasterize=true)
    ax3 = Axis(fig[3,1], title=L"\exp(-\mathrm{i}S)F(t)d[k+A(t_1)]",
               xlabel=L"$\Re[t_1]$ [jiffies]", ylabel=L"$\Im[t_1]$ [jiffies]",
               xgridvisible=false, ygridvisible=false)
    domaincolor!(ax3, z -> exp(-im*Sf(z))*field_amplitude(F, z)*dcanf(z), dlim, pixels=pixels, abs=Inf, rasterize=true)

    plot_saddle_points!(ax1, ax2, ax3)

    # save("sfa-ati-action-one-momentum-complex.pdf", fig)
    fig
end

# ╔═╡ b47d87c4-3e42-41d9-bc2e-6b32007e5f1e
let
    ∂Sfk = Base.Fix2(∂Sf, kkx)
    rf = get_r(kkx)

    fig = Figure(size=(800,1000))

    ax1 = Axis(fig[1,1][1,1], title=L"\partial S",
               xticklabelsvisible=false, ylabel=L"$\Im[t_1]$ [jiffies]")
    domaincolor!(ax1, ∂Sfk, dlim, pixels=pixels, rasterize=true)

    ax2 = Axis(fig[1,1][2,1], title="Rational approximation",
               xticklabelsvisible=false, ylabel=L"$\Im[t_1]$ [jiffies]")
    domaincolor!(ax2, rf, dlim, pixels=pixels, rasterize=true)

    ax3 = Axis(fig[1,1][3,1], title=L"\log_{10}\Delta", #
               # xlims=(dlim[1],dlim[2]),ylims=(dlim[3],dlim[4])
               xlabel=L"$\Re[t_1]$ [jiffies]", ylabel=L"$\Im[t_1]$ [jiffies]")
    hm = heatmap!(ax3, t, tim, log10.(abs2.(∂Sfk.(tc) .- rf.(tc)))',
                  colormap=:plasma, rasterize=true)
    Colorbar(fig[1,2][3,1], hm)

    plot_saddle_points!(ax1, ax2, ax3)

    # save("sfa-ati-saddle-points-one-momentum-complex.pdf", fig)
    fig
end

# ╔═╡ 562f4b5a-0319-48ec-98cb-6cf65820748f
let
    Sf = z -> StrongFieldApproximation.volkov_phase(kkx, vpc, z) + Iₚ*(t[end] - z)
    ∂Sfk = Base.Fix2(∂Sf, kkx)

    fig = Figure(size=(800,650))

    ax1 = Axis(fig[1,1], title=L"S",
               xticklabelsvisible=false, ylabel=L"$\Im[t_1]$ [jiffies]")
    domaincolor!(ax1, Sf, dlim, pixels=pixels, rasterize=true)

    ax2 = Axis(fig[2,1], title=L"\partial S",
               xlabel=L"$\Re[t_1]$ [jiffies]", ylabel=L"$\Im[t_1]$ [jiffies]")
    domaincolor!(ax2, ∂Sfk, dlim, pixels=pixels, rasterize=true)

    plot_saddle_points!(ax1, ax2)

    # save("sfa-ati-action-and-saddle-points-one-momentum-complex.pdf", fig)
    fig
end

# ╔═╡ 2714f784-a671-4661-9ce6-a8594e2ae85f
let
    pcan = z -> kkx + vector_potential(F, z)
    dcan = z -> d(pcan(z))
    Fdcan = z -> field_amplitude(F,z) * dcan(z)

    fig = Figure(size=(800,800))

    ax1 = Axis(fig[1,1],
               title=L"k_x + A(t_1)",
               xticklabelsvisible=false, ylabel=L"$\Im[t_1]$ [jiffies]")
    domaincolor!(ax1, pcan, dlim, abs=Inf, pixels=pixels, rasterize=true)

    ax2 = Axis(fig[2,1],
               title=L"d[k_x + A(t_1)]",
               xlabel=L"$\Re[t_1]$ [jiffies]", ylabel=L"$\Im[t_1]$ [jiffies]")
    domaincolor!(ax2, dcan, dlim, abs=Inf, pixels=pixels, rasterize=true)

    plot_saddle_points!(ax1, ax2)

    # save("sfa-dipoles-at-saddle-points-one-momentum-complex.pdf", fig)

    fig
end

# ╔═╡ 48e09719-6e76-42ec-9c59-2de6c517ed7b
# ╠═╡ disabled = true
#=╠═╡
let
    fig = Figure(size=(600,400))
    ax = Axis(fig[1,1])

    f = d
    # f = z -> z/(z^2 + 1)^3
    rf = approximate(f, unit_circle)
    @info "Rational approximation" rf roots(rf) poles(rf)

    domaincolor!(ax, f, 2, abs=Inf, pixels=pixels, rasterize=true)
    fig
end
  ╠═╡ =#

# ╔═╡ 2822bf7d-429a-492c-83b6-57f5d42eadef
begin
    saddle_points = Vector{Vector{ComplexF64}}()
    for kk in kx
        rf = get_r(kk)
        push!(saddle_points, filter(z -> imag(z) > 0, roots(rf)))
    end
    saddle_points
end

# ╔═╡ 3001c781-54a1-4a07-be8d-68309009a107
let
    function plot_saddle_points!(ax, kx, saddle_points)
        kmin,kmax = extrema(kx)
        s = nothing
        for (j,kk) in enumerate(kx)
            sj = saddle_points[j]
            s = scatter!(ax, real(sj), imag(sj),
                     color=kk*ones(length(sj)),
                     colorrange=(kmin,kmax),
                     strokewidth=0, markersize=6,
                     colormap=:RdBu)
        end

        limits!(ax, dlim[1], dlim[2], 0, dlim[4])
        s
    end

    sp = findall(>(0), kx)
    sm = findall(<(0), kx)

    fig = Figure(size=(700,500))

    # ax1 = Axis(fig[1,1],
    #            xticklabelsvisible=false, ylabel=L"$\Im[t_1]$ [jiffies]")
    # s1 = plot_saddle_points!(ax1, kx[sp], saddle_points[sp])
    # Colorbar(fig[1,2], s1, label=L"$k$ [Bohr$^{-1}$]")

    # ax2 = Axis(fig[2,1],
    #            xlabel=L"$\Re[t_1]$ [jiffies]", ylabel=L"$\Im[t_1]$ [jiffies]")
    # s2 = plot_saddle_points!(ax2, kx[sm], saddle_points[sm])
    # Colorbar(fig[2,2], s2, label=L"$k$ [Bohr$^{-1}$]")

    ax1 = Axis(fig[1,1],
               xticklabelsvisible=false, ylabel=L"$\Im[t_1]$ [jiffies]")
    s1 = plot_saddle_points!(ax1, kx, saddle_points)
    Colorbar(fig[1,2], s1, label=L"$k$ [Bohr$^{-1}$]")

    ax2 = Axis(fig[2,1],
               xlabel=L"$\Re[t_1]$ [jiffies]", ylabel=L"$F(t_1)$ [au]",
               limits=((dlim[1],dlim[2]), nothing))
    lines!(ax2, t, Fv)

    # save("sfa-ati-saddle-points-complex.pdf", fig)
    fig
end


# ╔═╡ ee9cecf7-640e-479c-9610-612715feed4f
md"# HHG saddle points"

# ╔═╡ 00616098-d6fb-40c8-80c7-06446d12b707
let
    Ω = 0.5
    nt = length(t)
    Z = zeros(nt,nt)

    for (i,tf) in enumerate(t)
        Af = vector_potential(F, tf)
        for (j,ti) in enumerate(t)
            ti < tf || continue
            Ai = vector_potential(F, tf)
            kst = StrongFieldApproximation.stationary_momentum(vp, ti, tf)
            Z[i,j] = abs(kst + Af) + abs(kst + Ai) - 2Ω
        end
    end

    Zmax = maximum(abs, filter(!isnan, Z))

    fig = Figure(size=(700,700))
    ax = Axis(fig[1,1], xlabel=L"$t_i$ [fs]", ylabel=L"$t$ [fs]")
    Colorbar(fig[1,2],
             heatmap!(ax, ustrip.(tplot), ustrip.(tplot), Z,
                      colorrange=(-Zmax,Zmax), colormap=:RdBu))
    fig
end

# ╔═╡ e2ad5d65-62fc-49ba-a774-c337307b7619
md"# Code loading"

# ╔═╡ Cell order:
# ╟─8e113bbb-6d30-481d-884b-11ecc81ae6e8
# ╟─99655ba4-b323-41d3-ad36-484db298ebe4
# ╟─bab736ac-5aa5-478b-88b7-c20ee328f122
# ╟─95624be6-9f58-4a0e-953b-545d1a9caaa0
# ╟─29554dea-2b64-4976-b850-48b5aaeca0b8
# ╟─4c28f702-f278-4337-bd5e-6f8a16e723d2
# ╟─180fd636-4acd-4d74-af7a-c520e83a407c
# ╟─43789066-1806-400e-a94d-cb66c11eff78
# ╟─f5551674-12ab-4386-925e-de85043450f6
# ╟─f8c36bf3-6ef9-4039-8acb-09fb3d6ba93e
# ╟─d9d59e03-19ff-4669-8e83-6ddf73a5982b
# ╟─75169dc6-99a0-44e6-a8cb-520aa6c633b2
# ╟─fbee840d-14ed-4fd3-86d4-66c2c53e5ebd
# ╟─59c9b684-d5a9-4d98-9706-539bc2d7ea90
# ╟─ebc718dc-d506-4168-8b04-ee696aa1e1ab
# ╟─8578f4a0-e50d-4789-a5e4-5a621e3f066c
# ╟─aa86939e-7add-4379-ae4f-29c51dbf7cfa
# ╟─bab8dc3b-7755-4487-875a-f189ed40335a
# ╟─1da023e0-b0d6-48e2-88f0-97809af5ef03
# ╟─4c3d86a2-08d3-4c4d-aea3-b3e2204279c8
# ╟─978c2eeb-d9af-47bc-8340-e6f411e392f3
# ╟─916feb89-509e-466e-ba7e-8a27b90a986d
# ╟─ef12d365-2b25-42a5-8c29-0ce5ec1a115c
# ╟─b06bc47a-f5eb-4b0f-961c-fdb50618a0f9
# ╟─be1211d9-4fed-4e48-8e7a-13b27204cff0
# ╠═8ed2ebf7-43ff-42c8-8299-99ca0bb94377
# ╟─26465e4e-a6e6-42ed-b5a6-3da23421bb01
# ╟─5cdb2f14-58a4-4b33-bedb-7a991c51d059
# ╟─ede51caa-e20b-4420-9193-158aec3e830d
# ╟─053c9566-ca55-440d-9544-905944212774
# ╟─92adcd20-eecd-4737-bf25-6dd76ce878c7
# ╟─d6d87b41-55a2-4b19-a6b1-715beb770e32
# ╟─b47d87c4-3e42-41d9-bc2e-6b32007e5f1e
# ╟─562f4b5a-0319-48ec-98cb-6cf65820748f
# ╟─2714f784-a671-4661-9ce6-a8594e2ae85f
# ╟─48e09719-6e76-42ec-9c59-2de6c517ed7b
# ╟─2822bf7d-429a-492c-83b6-57f5d42eadef
# ╟─3001c781-54a1-4a07-be8d-68309009a107
# ╟─ee9cecf7-640e-479c-9610-612715feed4f
# ╟─00616098-d6fb-40c8-80c7-06446d12b707
# ╟─e2ad5d65-62fc-49ba-a774-c337307b7619
# ╟─ce907c44-0391-11ef-123a-1d651bb74fcc
