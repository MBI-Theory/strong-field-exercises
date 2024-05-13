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

# ╔═╡ c509765c-0079-11ef-10e0-29f04d6cc915
begin
    using Pkg
    Pkg.activate(".")

    using PlutoUI

    using Plots
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

# ╔═╡ e95a3092-5273-4d97-9dc6-221f79f45bbc
md"# HHG and ATI using TDSE and SFA"

# ╔═╡ a535b157-3a44-42ae-99d1-c11409d7cb28
md"Previous calculations:"

# ╔═╡ 5a92b770-33f8-45fd-8164-09bde6eb16f0
@bind run_preset presets_input()

# ╔═╡ cb8d72d1-bc2e-4d0c-b659-e00fed62239c
@bind inputs saved_inputs(run_preset) do sv
    is = InputSection("Calculation inputs", [
        Input("run_dir", TextField(default=run_preset), "Run directory"),
        InputSection("Electric field", [
            Input("I₀", NumberIntervalField(0.0..1000, default=sv("I₀", 300.0)), "Intensity", unit="TW/cm²"),
            Input("λorω", Select([:λ => "λ", :ω => "ħω"], default=Symbol(sv("λorω", :λ))), "Carrier"),
            Input("λ", NumberIntervalField(1.0..1000, default=sv("λ", 800.0)), "Wavelength", unit="nm"),
            Input("ω", NumberIntervalField(0.01..10, default=sv("ω", 1.0)), "``\\hbar\\omega``", unit="Ha"),
            Input("cycles", NumberIntervalField(1.0..30, default=sv("cycles", 3.0)), "Duration", unit="cycles")
        ]),
        InputSection("Atom", [
            Input("potential", Select([:coulomb => "Coulomb, V(r)=-Z/r", :yukawa => "Yukawa, V(r)=-Z exp(-Ar)"], default=Symbol(sv("potential", :coulomb))), "Potential"),
            Input("charge", NumberIntervalField(0.0..137, default=sv("charge", 1)), "Charge ``Z``"),
            Input("yukawa_A", NumberIntervalField(0.0..137, default=sv("yukawa_A", 1)), "Yukawa ``A``")
        ]),
        InputSection("Time steps", [
            InputSection("TDSE", [
                Input("δt", NumberIntervalField(0.001..0.5, default=sv("δt", 0.2), step=0.001), "Max time step", unit="jiffies"),
                Input("δtobs", NumberIntervalField(0.001..50, default=sv("δtobs", 0.5), step=0.001), "Observables time step", unit="jiffies"),
                Input("default_sampling_factor", NumberIntervalField(1.0..1000, default=sv("default_sampling_factor", 10.0)), "Sampling factor (10 is good)"),
                Input("auto_decrease_timestep", CheckBox(default=sv("auto_decrease_timestep", true)), "Auto decrease time step",
                      description="The code might decide to decrease the time step, if deemed necessary."),
            ]),
            InputSection("SFA", [
                Input("sfa_ndt", NumberField(1:700, default=sv("sfa_ndt", 350)), "Time steps per cycle"),
                Input("sfa_memory", NumberIntervalField(0.0..4, default=sv("sfa_memory", 0.65)), "SFA integral memory", unit="cycles"),
                Input("sfa_rescattered", CheckBox(default=sv("sfa_rescattered", false)), "Compute SFA rescattered ATI (on-axis only)")
            ])
        ]),
        InputSection("Photoelectron spectra", [
            Input("volkov_tsurff", CheckBox(default=sv("volkov_tsurff", true)), "Volkov spectrum",
                  description="Project wavefunction onto Volkov waves during time propagation."),
            Input("volkov_isurff", CheckBox(default=sv("volkov_isurff", true)), "Infinite-time correction",
                  description="Compute correction due to flux that did not leave the computational box during time propagation."),
            Input("coulomb_isurff", CheckBox(default=sv("coulomb_isurff", false)), "Coulomb spectrum",
                  description="Project wavefunction onto Coulomb waves after time propagation; requires the whole wavefunction to remain inside the computational box."),
            Input("surff_kmax", NumberIntervalField(0.01..10, default=sv("surff_kmax", 4.0)), "``k_{\\textrm{max}}``"),
            Input("surff_nk", NumberField(2:1000, default=sv("surff_nk", 300)), "``n_k``"),
            Input("surff_nθ", NumberField(1:150, default=sv("surff_nθ", 50)), "``n_\\theta``"),
            Input("surff_nϕ", NumberField(1:300, default=sv("surff_nϕ", 1)), "``n_\\phi``")
        ], detail=true),
        InputSection("Computational box", [
            Input("rmax", NumberIntervalField(33.0..1000, default=sv("rmax", 60.0), step=0.5), "Radial box size", unit="Bohr"),
            Input("dr", NumberIntervalField(0.05..0.5, default=sv("dr", 0.2), step=0.05), "Radial step", unit="Bohr"),
            Input("ℓmax", NumberField(1:100, default=sv("ℓmax", 6)), "Maximum ``\\ell``"),
            InputSection(" Complex absorbing potential (CAP)", [
                Input("cap_kind", Select(["Manolopoulos", "none"], default=sv("cap_kind", "Manolopoulos")), "CAP kind")
                Input("cap_kmin", NumberIntervalField(0.0..1, default=sv("cap_kmin", 0.2)), "``k_{\\textrm{min}}``", unit="Bohr⁻¹",
                      description=md"""The lowest momentum to be absorbed with the guaranteed efficiency. All higher momenta will be absorbed more efficiently. Smaller values of ``k_{\textrm{min}}`` result in thicker absorbing boundary.""")
                Input("cap_delta", NumberIntervalField(0.0..1, default=sv("cap_delta", 0.2)), "``\\delta``",
                      description=md"The JWKB scaling parameter. The default is 0.2, corresponding to 1% reflection probability for ``k=k_{\textrm{min}}``. Smaller values of ``\delta`` decrease the reflection probability, at the cost of a thicker absorbing boundary. For further discussion of the ``k_{\textrm{min}}`` and ``\delta`` parameters, see the [original publication](http://dx.doi.org/10.1063/1.1517042).")
            ], detail=true),
        ]),
        InputSection("Debugging", [
            Input("overwrite", CheckBox(default=sv("overwrite", false)), "Overwrite previous results"),
            Input("debug_plots", CheckBox(default=sv("debug_plots", false)), "Debug plots"),
            Input("composition_threshold", NumberIntervalField(-1.0..1, default=sv("composition_threshold", 1e-5)), "Composition threshold"),
            Input("omp_num_threads", NumberField(1:16, default=sv("omp_num_threads", 8)), "`OMP_NUM_THREADS`"),
            Input("omp_stacksize", TextField(default=sv("omp_stacksize", "65520k")), "`OMP_STACKSIZE`"),
            Input("run_tdse", CheckBox(default=sv("run_tdse", true)), "Run TDSE calculation", description="(If unchecked: only generate SCID input)")
        ], detail=true),
    ])

    confirm(format_input(is), label="Save inputs and run calculation")
end

# ╔═╡ 44649b9f-6629-41ad-a4d1-3ab93703c0c9
inputs

# ╔═╡ 4f5b0481-a592-469e-9ab0-1664703aca03
begin
    I₀ = inputs.I₀*u"TW/cm^2"
    λ = inputs.λ*u"nm"
    ω = inputs.ω
    F = if inputs.λorω == :λ
        @field(F) do
            I₀ = I₀
            λ = λ
            τ = inputs.cycles*u"fs"(inputs.λ*u"nm"/u"c")
            σoff = 3.0
            σmax = 4.0
            env = :trunc_gauss
        end
    else
        @field(F) do
            I₀ = I₀
            ω = ω
            τ = inputs.cycles*(2π/inputs.ω)
            σoff = 3.0
            σmax = 4.0
            env = :trunc_gauss
        end
    end
end

# ╔═╡ b8478347-69eb-47b7-8531-23716024eb4f
potential = let Z = inputs.charge
    if isinteger(Z)
        Z = Int(Z)
    end
    if inputs.potential == :coulomb
        PointCharge(Z)
    elseif inputs.potential == :yukawa
        A = inputs.yukawa_A
        T = promote_type(typeof(Z), typeof(A))
        Yukawa(T(Z), T(A), one(T))
    else
        @error "Unknown potential kind $(inputs.potential)"
    end
end

# ╔═╡ 846e23f7-c430-48b2-92cd-7df3dc12aa44
md"# TDSE"

# ╔═╡ 523e76f6-dc84-4b2f-a270-1e824e75f1c8
run_dir = SCIDWrapper.run_calc(SCIDWrapper.System(F, inputs.δt, potential),
                               ENV["SCID_EXE"];
                               generate_input_only=!inputs.run_tdse,
                               inputs...)

# ╔═╡ dd127149-bf06-4910-933c-eb22a012243e
save_inputs(inputs, run_dir);

# ╔═╡ 464a422f-0083-41e6-a4fd-14a4088ab689
results = isnothing(run_dir) ? nothing : SCIDWrapper.load_calculation(run_dir);

# ╔═╡ 27d02d96-d85b-465b-a5a2-96ced9efae69
md"# SFA"

# ╔═╡ bd960418-b383-4cd9-a268-be4e1acad804
# If we have TDSE results, take the ionization potential as the
# binding energy of the initial state, otherwise we use the Coulombic
# value.
sfa_Iₚ = (isnothing(results) ? inputs.charge/2 : -real(results.E₀))*u"hartree"

# ╔═╡ 58c7a714-8f93-407c-a7f4-3c47bfec2424
sfa_Uₚ = ponderomotive_potential(F) |> u"eV"

# ╔═╡ e25b6ead-e745-48ba-9e7f-329c101ceec9
sfa_hhg_cutoff = 3.17sfa_Uₚ + sfa_Iₚ |> u"eV"

# ╔═╡ 8ce0c0c2-d1a3-4bc8-8774-16bf451a33af
sfa_t = timeaxis(F, inputs.sfa_ndt)

# ╔═╡ 683182fe-cbe4-45fb-b36d-a9b0ff6a79f3
Fv = field_amplitude(F, sfa_t)

# ╔═╡ 08282bff-8f46-4185-a2a9-a9a43d39e1b7
k,kmag,θ = momentum_grid(0, inputs.surff_kmax, inputs.surff_nk, inputs.surff_nθ, spacing=:momentum)

# ╔═╡ 5ebac048-1163-4399-bb88-cdec5a267163
sfa_system = let
    ic = StrongFieldApproximation.IonizationChannel(sfa_Iₚ, F, inputs.sfa_ndt)
    cc = StrongFieldApproximation.CoulombCoupling((𝐤,𝐩) -> yukawa_fourier(𝐩-𝐤, 1, 0, 1))
    StrongFieldApproximation.System([ic], [[cc;;]], F, inputs.sfa_ndt)
end

# ╔═╡ c6b2c550-7715-4ee5-802a-f4bc532f90d8
hhg_direct_diagram = Diagram([(1,0),(1,0)], sfa_system)

# ╔═╡ 7b7b1ef0-d814-4c1d-af6f-ccb0d38a9cd3
sfa_d = vec(cached_calculation(Float64, joinpath(run_dir, "sfa_hhg_direct.txt"), overwrite=inputs.overwrite) do
    induced_dipole(sfa_system, hhg_direct_diagram,
                   memory=floor(Int, inputs.sfa_memory*inputs.sfa_ndt))
end)

# ╔═╡ baac7bb3-1e36-4007-b450-84561b6803fd
ati_direct_diagram = Diagram(sfa_system)

# ╔═╡ 7ff2c1ce-7a93-46d4-84d7-1dd3340f27df
sfa_c_direct = cached_calculation(ComplexF64, joinpath(run_dir, "sfa_c_direct.txt"), overwrite=inputs.overwrite) do
    photoelectron_spectrum(k, sfa_system, ati_direct_diagram)
end

# ╔═╡ 6ecabf08-446d-49fe-96ca-d32f99eb207d
ati_rescattering_diagram = Diagram([(1,1),(1,0)], sfa_system)

# ╔═╡ 7f3629db-943d-4fe7-a0bc-f7ab9a77408d
sfa_c_rescattered = if inputs.sfa_rescattered
    cached_calculation(ComplexF64, joinpath(run_dir, "sfa_c_rescattered.txt"), overwrite=inputs.overwrite) do
        photoelectron_spectrum(kmag, sfa_system, ati_rescattering_diagram, memory=1*inputs.sfa_ndt)
    end
else
    @warn "Calculation of SFA rescattered ATI disabled"
end

# ╔═╡ a176e7cd-a3f4-4d2b-94ea-b8b1cb36601b
md"# Results"

# ╔═╡ 22577c81-2684-436e-9fa6-1c70bfe26021
md"""## (High-order) Harmonic Generation

The radiation emitted is proportional to the acceleration of the
induced dipole moment. SCID computes the time-dependent dipole moment
``\langle\vec{d}(t)\rangle``, its velocity
``\langle\vec{v}(t)\rangle=\partial_t\langle\vec{d}(t)\rangle``, and
its acceleration
``\langle\vec{a}(t)\rangle=\partial_t^2\langle\vec{d}(t)\rangle``. By
Fourier-transforming these time series and using the Fourier identity
``\partial_t \leftrightarrow \mathrm{i}\omega``, we can study the
spectrum of the emitted radiation.

Since we solved the TDSE for a linearly polarized field, only the
``z`` component is non-zero.
"""

# ╔═╡ 0947d3b4-2c76-4de4-a6cf-f91ae840c4ed
let
    ptdse = SCIDWrapper.plot_dipole_moment(results, title="TDSE")

    Av = vector_potential(F, sfa_t)
    psfa = SCIDWrapper.plot_dipole_moment(sfa_t, Fv, Av, -sfa_d, title="SFA")

    plot(ptdse, psfa, size=(900,800))
end

# ╔═╡ 1c3b3451-7348-43ee-9f8c-5161c36f2279
begin
    E = @bind energy_kind Select([:total => "Total energy", :kinetic => "Kinetic energy"])
    u = @bind freq_unit Select([u"hartree" => "Ha ≈ 27.211 eV", u"eV" => "eV", "Up" => "Ponderomotive potential",
                                "q" => "Harmonic order", u"nm" => "Wavelength (nm)"])
    md"Preferred unit for frequency axes: $(E) $(u)"
end

# ╔═╡ 57ab7dae-a0ee-494f-9bfa-9edb264de363
begin
    aw = @bind apodizing_window Select([hanning => "Hann", hamming => "Hamming",
                                        blackman => "Blackman",
                                        rect => "Rectangular"])
    md"[Window function](https://en.wikipedia.org/wiki/Window_function): $(aw)"
end

# ╔═╡ 98dce21a-7df0-4059-857c-2777c6e145e0
let
    kw = (;ωkind=energy_kind, ωunit=freq_unit, window=apodizing_window)
    ptdse = SCIDWrapper.plot_dipole_spectrum(results; kw..., title="TDSE")

    xl = xlims(ptdse)

    ω₀ = austrip(photon_energy(F))
    δt = step(sfa_t)

    psfa = SCIDWrapper.plot_dipole_spectrum(sfa_t, δt, ω₀,
                                            austrip(sfa_Iₚ), austrip(sfa_Uₚ), austrip(sfa_hhg_cutoff),
                                            -sfa_d;
                                            kw..., xlims=xl, title="SFA")

    plot(ptdse, psfa, size=(900,900),
         layout=@layout([a;b]))
end

# ╔═╡ 2ca2bf08-51ab-40b1-a5d3-eebb7c11b1f4
begin
    wl = @bind gabor_window_length NumberIntervalField(0.0..3, default=0.05)
    dr = @bind gabor_dynamic_range NumberIntervalField(0.0..10, default=3)
    md"Time–frequency analysis, window length: $(wl) cycles, dynamic range: $(dr) orders of magnitude"
end

# ╔═╡ 6992276e-f690-4c4b-98e7-5e08bc08f745
let

    kw = (ωkind=energy_kind, ωunit=freq_unit,
          window_length=gabor_window_length)
    ptdse = plot(SCIDWrapper.plot_time_frequency_analysis(results, F; kw..., dynamic_range=gabor_dynamic_range),
                 title="TDSE")

    yl = ylims(ptdse)
    psfa = plot(SCIDWrapper.plot_time_frequency_analysis(sfa_t, sfa_d, F,
                                                         austrip(sfa_Iₚ), austrip(sfa_Uₚ), austrip(sfa_hhg_cutoff);
                                                         kw..., dynamic_range=gabor_dynamic_range+1.5),
                ylims=yl,
                title="SFA")

    plot(ptdse, psfa, size=(900,1000),
         layout=@layout([a;b]))
end

# ╔═╡ 021d2604-5ce5-4858-827a-06f7f72d204b
md"""
## Population decomposition

At the end of the calculation, we may project the wavefunction onto
angular symmetries (``\ell`` and ``m_\ell``) to analyze its structure,
as well as the eigenstates of the field-free Hamiltonian.
"""

# ╔═╡ 81fe22ab-dfe7-4605-9a14-e44db976218f
md"### Angular momentum decomposition"

# ╔═╡ dff4b697-ec4d-42f2-8fbd-3309eb108879
if !isnothing(results) && !isnothing(results.final_population_m_resolved)
    real.(results.final_population_m_resolved)
end

# ╔═╡ ef8559d3-57ef-49f7-8c75-7dd93de41c31
SCIDWrapper.plot_angular_decomposition(results)

# ╔═╡ cfc8fab1-2695-40ef-8d8c-54851fd5edb1
md"### Eigendecomposition"

# ╔═╡ 4a00dcb9-9885-432b-993e-e495ae2c7d04
if !isnothing(results) && !isnothing(results.large_amplitudes)
    real.(results.large_amplitudes)
end

# ╔═╡ 78b0682a-d1e9-41cd-89a6-34ffa5dd6ac3
SCIDWrapper.plot_eigen_decomposition(results)

# ╔═╡ 7055184b-cc62-4f21-9381-3974ceaff888
begin
    dp = @bind draw_photons CheckBox(default=true)
    md"Draw photons: $(dp)"
end

# ╔═╡ f3328b5b-9540-4a65-8463-69e7b63bd63d
SCIDWrapper.plot_photon_diagram(results, draw_photons=draw_photons)

# ╔═╡ e4a48216-1551-4214-a127-69d46cf97864
md"## Photoelectron spectra"

# ╔═╡ 8757101f-5361-4072-b054-23fcd3af9bc9
let
    dr = @bind pes_dynamic_range NumberIntervalField(0.0..10, default=6)
    projection = @bind pes_projection Select([:box => "Box", :polar => "Polar"])
    md"Photoelectron spectrum, dynamic range: $(dr) orders of magnitude, projection: $(projection)"
end

# ╔═╡ 1b8115bb-1b44-488a-b1db-a2e14e4f95af
let
    plot(SCIDWrapper.plot_pes(results.volkov_pes,
                              dynamic_range=pes_dynamic_range, projection=pes_projection,
                              Uₚ=results.Uₚ), title="TDSE")
end

# ╔═╡ 74cc9cf1-152c-4be3-afc5-90785a1c1fb0
let
    sfa_pes = (A=sfa_c_direct,θ=rad2deg.(θ),k=kmag)

    p = plot(SCIDWrapper.plot_pes(sfa_pes,
                                  dynamic_range=pes_dynamic_range, projection=pes_projection,
                                  Uₚ=sfa_Uₚ), title="SFA direct only")
    if !isnothing(sfa_c_rescattered)
        plot!(p[2], kmag, abs2.(sfa_c_rescattered), label=L"Rescattered, $\theta=0^\circ$",
            title="SFA, direct and rescattered")
    end
    p
end

# ╔═╡ c589c676-b22f-455e-a63f-39131123d5b1
md"# Helper code"

# ╔═╡ 6bf8557c-f7fd-4059-8e29-a430400577c2
notebook_styling()

# ╔═╡ 5e3d475d-54e9-4e40-af42-764a1f8759aa


# ╔═╡ Cell order:
# ╟─e95a3092-5273-4d97-9dc6-221f79f45bbc
# ╟─a535b157-3a44-42ae-99d1-c11409d7cb28
# ╟─5a92b770-33f8-45fd-8164-09bde6eb16f0
# ╟─cb8d72d1-bc2e-4d0c-b659-e00fed62239c
# ╟─44649b9f-6629-41ad-a4d1-3ab93703c0c9
# ╟─4f5b0481-a592-469e-9ab0-1664703aca03
# ╟─b8478347-69eb-47b7-8531-23716024eb4f
# ╟─846e23f7-c430-48b2-92cd-7df3dc12aa44
# ╟─523e76f6-dc84-4b2f-a270-1e824e75f1c8
# ╟─dd127149-bf06-4910-933c-eb22a012243e
# ╟─464a422f-0083-41e6-a4fd-14a4088ab689
# ╟─27d02d96-d85b-465b-a5a2-96ced9efae69
# ╟─bd960418-b383-4cd9-a268-be4e1acad804
# ╟─58c7a714-8f93-407c-a7f4-3c47bfec2424
# ╟─e25b6ead-e745-48ba-9e7f-329c101ceec9
# ╟─8ce0c0c2-d1a3-4bc8-8774-16bf451a33af
# ╟─683182fe-cbe4-45fb-b36d-a9b0ff6a79f3
# ╟─08282bff-8f46-4185-a2a9-a9a43d39e1b7
# ╟─5ebac048-1163-4399-bb88-cdec5a267163
# ╟─c6b2c550-7715-4ee5-802a-f4bc532f90d8
# ╟─7b7b1ef0-d814-4c1d-af6f-ccb0d38a9cd3
# ╟─baac7bb3-1e36-4007-b450-84561b6803fd
# ╟─7ff2c1ce-7a93-46d4-84d7-1dd3340f27df
# ╟─6ecabf08-446d-49fe-96ca-d32f99eb207d
# ╟─7f3629db-943d-4fe7-a0bc-f7ab9a77408d
# ╟─a176e7cd-a3f4-4d2b-94ea-b8b1cb36601b
# ╟─22577c81-2684-436e-9fa6-1c70bfe26021
# ╟─0947d3b4-2c76-4de4-a6cf-f91ae840c4ed
# ╟─1c3b3451-7348-43ee-9f8c-5161c36f2279
# ╟─57ab7dae-a0ee-494f-9bfa-9edb264de363
# ╟─98dce21a-7df0-4059-857c-2777c6e145e0
# ╟─2ca2bf08-51ab-40b1-a5d3-eebb7c11b1f4
# ╟─6992276e-f690-4c4b-98e7-5e08bc08f745
# ╟─021d2604-5ce5-4858-827a-06f7f72d204b
# ╟─81fe22ab-dfe7-4605-9a14-e44db976218f
# ╟─dff4b697-ec4d-42f2-8fbd-3309eb108879
# ╟─ef8559d3-57ef-49f7-8c75-7dd93de41c31
# ╟─cfc8fab1-2695-40ef-8d8c-54851fd5edb1
# ╟─4a00dcb9-9885-432b-993e-e495ae2c7d04
# ╟─78b0682a-d1e9-41cd-89a6-34ffa5dd6ac3
# ╟─7055184b-cc62-4f21-9381-3974ceaff888
# ╟─f3328b5b-9540-4a65-8463-69e7b63bd63d
# ╟─e4a48216-1551-4214-a127-69d46cf97864
# ╟─8757101f-5361-4072-b054-23fcd3af9bc9
# ╟─1b8115bb-1b44-488a-b1db-a2e14e4f95af
# ╟─74cc9cf1-152c-4be3-afc5-90785a1c1fb0
# ╟─c589c676-b22f-455e-a63f-39131123d5b1
# ╟─c509765c-0079-11ef-10e0-29f04d6cc915
# ╟─6bf8557c-f7fd-4059-8e29-a430400577c2
# ╟─5e3d475d-54e9-4e40-af42-764a1f8759aa
