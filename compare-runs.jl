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

# ╔═╡ ee488d3d-b9d3-458a-a8ee-37f43f5d7b2a
begin
    using Pkg
    Pkg.activate(".")

    using PlutoUI

    using CairoMakie
    using LaTeXStrings
    using DelimitedFiles

    using IntervalSets
    using Statistics
    using FFTW
    using AbstractFFTs
    using DSP

    using Unitful
    using UnitfulAtomic

    using ElectricFields

    using Revise

    using NotebookHelper

    using SCIDWrapper
    using StrongFieldApproximation
end

# ╔═╡ 2c5978fe-9f86-4fab-98b6-4c77b02ab4ed
@bind runs select_previous_runs(extra_files=["output.detail", "run.inp", "run.out", "volkov.pes",
                                             "sfa_hhg_direct.txt", "sfa_c_direct.txt"])

# ╔═╡ 858cee8f-78fb-44a2-aeb5-49c02f4374aa


# ╔═╡ b00dce58-620d-4234-b386-8cd8f673a6e7
all_results = let
    au2fs = auconvert(u"fs", 1)
    map(runs) do r
        results = SCIDWrapper.load_calculation(r)
        inputs = load_inputs(r)

        F = get_electric_field(inputs)
        ω₀ = austrip(photon_energy(F))

        sfa_t = timeaxis(F, inputs["sfa_ndt"])
        sfa_d = readdlm(joinpath(r, "sfa_hhg_direct.txt"),',',Float64,'\n')

        k,kmag,θ = momentum_grid(0, inputs["surff_kmax"], inputs["surff_nk"], inputs["surff_nθ"], spacing=:momentum)

        sfa_c_direct = readdlm(joinpath(r, "sfa_c_direct.txt"),',',ComplexF64,'\n')
        sfa_c_rescattered = let f = joinpath(r, "sfa_c_rescattered.txt")
            isfile(f) ? readdlm(f,',',ComplexF64,'\n') : nothing
        end

        sfa_Iₚ = (isnothing(results) ? inputs["charge"]/2 : -real(results.E₀))*u"hartree"
        sfa_Uₚ = ponderomotive_potential(F) |> u"eV"
        sfa_hhg_cutoff = 3.17sfa_Uₚ + sfa_Iₚ |> u"eV"

        merge(results, (;inputs=inputs,
                        F=F, ω₀=ω₀,
                        sfa_t=sfa_t,sfa_d=sfa_d,
                        sfa_tplot=au2fs*sfa_t,
                        sfa_Iₚ=sfa_Iₚ,
                        sfa_Uₚ=sfa_Uₚ,
                        sfa_hhg_cutoff=sfa_hhg_cutoff,
                        sfa_c_direct=sfa_c_direct,
                        sfa_c_rescattered=sfa_c_rescattered,
                        k=k, kmag=kmag, θ=θ))
    end
end;

# ╔═╡ fd11ba16-8899-4259-82e7-5e5047e0fb95
function plot_many_runs(fun::Function, v=runs; column_width=600, column_height=500,
                        layout=(1,length(v)))
    nr,nc = layout
    @assert nr*nc ≥ length(v)
    fig = Figure(size=(column_width*nc,column_height*nr))
    ci = CartesianIndices(layout)
    for i in eachindex(v)
        I = ci[i]
        fun(i, GridLayout(fig[I[1],I[2]]))
    end
    fig
end;

# ╔═╡ 29482791-5d88-474b-84fe-c7bbbbf08c09
md"# HHG"

# ╔═╡ 0ce086f7-f380-4a49-b1f4-bfca7c0f801f
let
    plot_many_runs(column_height=1000) do i, fig
        r = all_results[i]
        gl = GridLayout(fig[1:3,1])
        SCIDWrapper.plot_dipole_moment!(gl, r, title=runs[i])
        ax = last(contents(gl))
        ax.xaxis.attributes.labelvisible[] = false
        ax.xaxis.attributes.ticklabelsvisible[] = false

        ax = Axis(fig[4,1], title="SFA", xlabel=L"$t$ [fs]")
        lines!(ax, ustrip.(r.sfa_tplot), -vec(r.sfa_d))
    end
end

# ╔═╡ 38b3e22a-f80a-428e-a6ec-2e83338a67f4
begin
    E = @bind energy_kind Select([:total => "Total energy", :kinetic => "Kinetic energy"])
    u = @bind freq_unit Select([u"hartree" => "Ha ≈ 27.211 eV", u"eV" => "eV", "Up" => "Ponderomotive potential",
                                "q" => "Harmonic order", u"nm" => "Wavelength (nm)"])
    md"Preferred unit for frequency axes: $(E) $(u)"
end

# ╔═╡ da8b8473-b13c-493f-9616-9651eaa5472f
begin
    aw = @bind apodizing_window Select([hanning => "Hann", hamming => "Hamming",
                                        blackman => "Blackman",
                                        rect => "Rectangular"])
    md"[Window function](https://en.wikipedia.org/wiki/Window_function): $(aw)"
end

# ╔═╡ cd796434-38ff-48fb-a136-2c0089aabc43
let
    kw = (;ωkind=energy_kind, ωunit=freq_unit, window=apodizing_window)

    plot_many_runs(column_height=1000) do i, fig
        r = all_results[i]
        fig[1,1] = gl1 = GridLayout()
        SCIDWrapper.plot_dipole_spectrum!(gl1, r; kw..., title="TDSE",
                                          xticklabelsvisible=false, xlabelvisible=false)
                                              ax1 = first(contents(gl1[1,1]))
        xl = SCIDWrapper.get_limits!(ax1)[1]

        δt = step(r.sfa_t)

        fig[2,1] = gl2 = GridLayout()
        SCIDWrapper.plot_dipole_spectrum!(gl2, r.sfa_t, δt, r.ω₀,
                                          austrip(r.sfa_Iₚ), austrip(r.sfa_Uₚ), austrip(r.sfa_hhg_cutoff),
                                          -vec(r.sfa_d);
                                          kw..., limits=((xl...,), nothing),
                                          title="SFA")
        ax2 = first(contents(gl2[1,1]))
        linkxaxes!(ax1, ax2)
    end
end

# ╔═╡ 952d668e-9d2b-4ddc-b762-4a1a601bde0b
begin
    wl = @bind gabor_window_length NumberIntervalField(0.0..3, default=0.05)
    dr = @bind gabor_dynamic_range NumberIntervalField(0.0..10, default=3)
    md"Time–frequency analysis, window length: $(wl) cycles, dynamic range: $(dr) orders of magnitude"
end

# ╔═╡ de0abe8c-e1d0-427a-a446-9ed181e6890d
let
    kw = (ωkind=energy_kind, ωunit=freq_unit,
          window_length=gabor_window_length,
          maxsteps=500)

    plot_many_runs(column_height=1000) do i,fig
        r = all_results[i]
        gl1 = GridLayout(fig[1,1])
        SCIDWrapper.plot_time_frequency_analysis!(gl1, r, r.F; kw..., dynamic_range=gabor_dynamic_range, title="TDSE")

        ax1 = first(contents(gl1[1,1]))
        yl = SCIDWrapper.get_limits!(ax1)[2]
        ax1.xaxis.attributes.labelvisible[] = false
        ax1.xaxis.attributes.ticklabelsvisible[] = false

        gl2 = GridLayout(fig[2,1])
        SCIDWrapper.plot_time_frequency_analysis!(gl2, r.sfa_t, r.sfa_d, r.F,
                                                  austrip(r.sfa_Iₚ), austrip(r.sfa_Uₚ), austrip(r.sfa_hhg_cutoff);
                                                  kw..., dynamic_range=gabor_dynamic_range+1.5,
                                                  limits=(nothing, (yl...,)),

                                                  title="SFA")
    end
end

# ╔═╡ 5c0a2b02-1577-4145-9be0-cbaf925d2288
md"## Population decomposition"

# ╔═╡ 94297747-7dad-4631-9515-1ed3850a417c
md"### Angular momentum decomposition"

# ╔═╡ 417c71b5-7512-46cf-b30f-59a2d8901de2
plot_many_runs() do i,fig
    SCIDWrapper.plot_angular_decomposition!(fig, all_results[i])
end

# ╔═╡ 35813159-49b6-4072-aef0-dfa07901e89b
md"### Eigendecomposition"

# ╔═╡ cf180888-276e-4e27-ab21-8fde2d5d6372
plot_many_runs() do i,fig
    SCIDWrapper.plot_eigen_decomposition!(fig, all_results[i])
end

# ╔═╡ a81e0544-bd74-472e-a660-5e471211a7ab
begin
    dp = @bind draw_photons CheckBox(default=true)
    md"Draw photons: $(dp)"
end

# ╔═╡ 896f9649-3da3-416d-be29-aa618d6ec1fc
let
    plot_many_runs(column_height=500) do i, fig
        SCIDWrapper.plot_photon_diagram!(fig, all_results[i], draw_photons=draw_photons)
    end
end

# ╔═╡ 7a8f6c45-a03e-4406-8dc2-3aac325409e5
md"## Photoelectron spectra"

# ╔═╡ 9d2f0c3a-cc2a-4c39-8e20-882ef40415b2
let
    dr = @bind pes_dynamic_range NumberIntervalField(0.0..10, default=6)
    projection = @bind pes_projection Select([:box => "Box", :polar => "Polar"])
    md"Photoelectron spectrum, dynamic range: $(dr) orders of magnitude, projection: $(projection)"
end

# ╔═╡ 200a827a-5b6d-441d-ade8-fbc65aee6c88
let
    plot_many_runs(column_height=1000) do i, fig
        r = all_results[i]
        SCIDWrapper.plot_pes!(fig, r.volkov_pes,
                              dynamic_range=pes_dynamic_range, projection=pes_projection,
                              Uₚ=r.Uₚ, title="TDSE")
    end
end

# ╔═╡ eb8c3fb1-1ea3-4c6a-af0d-5ff167c34dac
let
    plot_many_runs(column_height=1000) do i, fig
        r = all_results[i]
        sfa_pes = (A=r.sfa_c_direct,θ=rad2deg.(r.θ),k=r.kmag)

        p = SCIDWrapper.plot_pes!(fig, sfa_pes,
                                  dynamic_range=pes_dynamic_range, projection=pes_projection,
                                  Uₚ=r.sfa_Uₚ, title="SFA direct only", nolegend=true)
        ax = first(contents(p[2,1]))
        if !isnothing(r.sfa_c_rescattered)
            lines!(ax, r.kmag, vec(abs2.(r.sfa_c_rescattered)), label=L"Rescattered, $\theta=0^\circ$")
            ax.title="SFA, direct and rescattered"
        end
        axislegend(ax, position=:lb)
    end
end

# ╔═╡ 4227dd14-1134-11ef-27e8-91c16eb7a3a3
md"# Helper code"

# ╔═╡ 77c362e5-04bc-4ba8-8aea-a5ba452aef37
notebook_styling()

# ╔═╡ Cell order:
# ╟─2c5978fe-9f86-4fab-98b6-4c77b02ab4ed
# ╟─858cee8f-78fb-44a2-aeb5-49c02f4374aa
# ╟─b00dce58-620d-4234-b386-8cd8f673a6e7
# ╟─fd11ba16-8899-4259-82e7-5e5047e0fb95
# ╟─29482791-5d88-474b-84fe-c7bbbbf08c09
# ╟─0ce086f7-f380-4a49-b1f4-bfca7c0f801f
# ╟─38b3e22a-f80a-428e-a6ec-2e83338a67f4
# ╟─da8b8473-b13c-493f-9616-9651eaa5472f
# ╟─cd796434-38ff-48fb-a136-2c0089aabc43
# ╟─952d668e-9d2b-4ddc-b762-4a1a601bde0b
# ╟─de0abe8c-e1d0-427a-a446-9ed181e6890d
# ╟─5c0a2b02-1577-4145-9be0-cbaf925d2288
# ╟─94297747-7dad-4631-9515-1ed3850a417c
# ╟─417c71b5-7512-46cf-b30f-59a2d8901de2
# ╟─35813159-49b6-4072-aef0-dfa07901e89b
# ╟─cf180888-276e-4e27-ab21-8fde2d5d6372
# ╟─a81e0544-bd74-472e-a660-5e471211a7ab
# ╟─896f9649-3da3-416d-be29-aa618d6ec1fc
# ╟─7a8f6c45-a03e-4406-8dc2-3aac325409e5
# ╟─9d2f0c3a-cc2a-4c39-8e20-882ef40415b2
# ╟─200a827a-5b6d-441d-ade8-fbc65aee6c88
# ╟─eb8c3fb1-1ea3-4c6a-af0d-5ff167c34dac
# ╟─4227dd14-1134-11ef-27e8-91c16eb7a3a3
# ╟─ee488d3d-b9d3-458a-a8ee-37f43f5d7b2a
# ╟─77c362e5-04bc-4ba8-8aea-a5ba452aef37
