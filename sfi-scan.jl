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

# ╔═╡ 0c7bb53c-05fb-11ef-0fd7-b157332ae1ef
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

# ╔═╡ bdafba4b-32c4-4665-9511-dce60b30947e
Fs = range(0.01, stop=0.1, length=11)

# ╔═╡ 153f2360-4777-4a0b-8d33-4051b905b514
let
    fsa = @bind field_strength_axis Select([:invF => "1/F", :Fsquared => "I = F²"])
    md"Field strength axis: $(fsa)"
end

# ╔═╡ 6f1985e8-0976-4136-8061-58c3748580bd
md"# Helper code"

# ╔═╡ 2b170dab-8dd9-4e4a-b15a-ab51d056d954
function scan_inputs()
    input_tree = InputSection("Scan inputs", [
        InputSection("Electric field", [
            Input("λ", NumberIntervalField(1.0..1000, default=400.0), "Wavelength", unit="nm"),
            Input("cycles", NumberIntervalField(1.0..30, default=3.0), "Duration", unit="cycles")
        ]),
        InputSection("Atom", [
            Input("charge", NumberIntervalField(0.0..137, default=1), "Charge ``Z``"),
        ]),
        InputSection("Time steps", [
            InputSection("TDSE", [
                Input("δt", NumberIntervalField(0.001..0.5, default=0.2, step=0.001), "Max time step", unit="jiffies"),
                Input("δtobs", NumberIntervalField(0.001..50, default=0.5, step=0.001), "Observables time step", unit="jiffies"),
                Input("default_sampling_factor", NumberIntervalField(1.0..1000, default=10.0), "Sampling factor (10 is good)"),
                Input("auto_decrease_timestep", CheckBox(default=true), "Auto decrease time step",
                      description="The code might decide to decrease the time step, if deemed necessary."),
            ])
        ]),
        InputSection("Computational box", [
            Input("rmax", NumberIntervalField(33.0..1000, default=60.0, step=0.5), "Radial box size", unit="Bohr"),
            Input("dr", NumberIntervalField(0.05..0.5, default=0.2, step=0.05), "Radial step", unit="Bohr"),
            Input("ℓmax", NumberField(1:100, default=6), "Maximum ``\\ell``"),
            InputSection(" Complex absorbing potential (CAP)", [
                Input("cap_kind", Select(["Manolopoulos", "none"]), "CAP kind")
                Input("cap_kmin", NumberIntervalField(0.0..1, default=0.2), "``k_{\\textrm{min}}``", unit="Bohr⁻¹",
                      description=md"""The lowest momentum to be absorbed with the guaranteed efficiency. All higher momenta will be absorbed more efficiently. Smaller values of ``k_{\textrm{min}}`` result in thicker absorbing boundary.""")
                Input("cap_delta", NumberIntervalField(0.0..1, default=0.2), "``\\delta``",
                      description=md"The JWKB scaling parameter. The default is 0.2, corresponding to 1% reflection probability for ``k=k_{\textrm{min}}``. Smaller values of ``\delta`` decrease the reflection probability, at the cost of a thicker absorbing boundary. For further discussion of the ``k_{\textrm{min}}`` and ``\delta`` parameters, see the [original publication](http://dx.doi.org/10.1063/1.1517042).")
            ], detail=true),
        ]),
        InputSection("Debugging", [
            Input("overwrite", CheckBox(default=false), "Overwrite previous results"),
            Input("debug_plots", CheckBox(default=false), "Debug plots"),
            Input("composition_threshold", NumberIntervalField(-1.0..1, default=1e-5), "Composition threshold"),
            Input("omp_num_threads", NumberField(1:16, default=8), "`OMP_NUM_THREADS`"),
            Input("omp_stacksize", TextField(default="65520k"), "`OMP_STACKSIZE`"),
            Input("run_tdse", CheckBox(default=true), "Run TDSE scan", description="(If unchecked: only generate SCID input)")
        ], detail=true),
    ])

    confirm(format_input(input_tree))
end;

# ╔═╡ 5fcfb1a9-6112-4060-8683-19d4484540bd
@bind inputs scan_inputs()

# ╔═╡ d8e7aab4-1046-44cb-af83-422848ab7979
potential = let Z = inputs.charge
    if isinteger(Z)
        Z = Int(Z)
    end
    PointCharge(Z)
end

# ╔═╡ 373c0414-38f2-4690-935a-0c8bd58330e1
function get_field(F₀, λ, cycles)
    @field(F) do
        E₀ = F₀
        λ = λ
        τ = inputs.cycles*u"fs"(inputs.λ*u"nm"/u"c")
        σoff = 4.0
        σmax = 6.0
        env = :trunc_gauss
    end
end

# ╔═╡ 8f71c6c7-3fec-4a7d-aecb-c919eab0a88b
scan_dir = SCIDWrapper.run_scan(Fs, ENV["SCID_EXE"], "scans/sfi-intensity-scan-$(inputs.λ) nm-$(inputs.cycles) cycles"; inputs...) do F₀
    F = get_field(F₀, inputs.λ*u"nm", inputs.cycles)
    SCIDWrapper.System(F, inputs.δt, potential)
end

# ╔═╡ 9b928f14-b785-436e-8841-8b7bda8092d6
results = SCIDWrapper.load_scan(Fs, scan_dir) do F₀,data
    ppt(args...) = StrongFieldApproximation.IonizationRates.ionization_yield(args...)
    sum_norms(v) = sum(real(e) for e in v if real(e) > 0)

    bn = sum_norms(data.final_population[!,"Bound population"])
    (F₀=F₀,
     bound_norm=bn,
     ionization_yield = 1 - bn,
     gst_norm = max(0, real(data.large_amplitudes[1,"Wgt"])),
     # ppt_ionization_yield = ppt(get_field(F₀, inputs.λ*u"nm", inputs.cycles), -real(data.E₀), 0, 0, inputs.charge)
    )
end

# ╔═╡ 169cacd0-d9c6-46ac-be84-2ff0ac1f48dc
let
    x,xlabel,xaxis = if field_strength_axis == :invF
        1 ./ Fs, L"$1/F$ [au]", :identity
    elseif field_strength_axis == :Fsquared
        Fs .^2, L"$F^2$ [au]", :log10
    else
        throw(ArgumentError("Unknown axis $(field_strength_axis)"))
    end
    p = plot(x, results.gst_norm, markershape=:auto, label="Ground state")
    plot!(p, x, results.bound_norm .- results.gst_norm, markershape=:auto, label="Excited states")
    plot!(p, x, results.ionization_yield, markershape=:auto, label="Ionization yield")
    # plot!(p, x, results.ppt_ionization_yield, markershape=:auto, label="PPT")
    plot(p, xlabel=xlabel, ylabel="Yield",
         xaxis=xaxis,
         yaxis=(:log10, (1e-14,:auto)),
         legend=:outerbottom, size=(900,700))
end

# ╔═╡ f2defed4-09d4-49d2-81a8-4c4e6b0d75f2
notebook_styling()

# ╔═╡ Cell order:
# ╟─bdafba4b-32c4-4665-9511-dce60b30947e
# ╟─5fcfb1a9-6112-4060-8683-19d4484540bd
# ╟─d8e7aab4-1046-44cb-af83-422848ab7979
# ╟─373c0414-38f2-4690-935a-0c8bd58330e1
# ╟─8f71c6c7-3fec-4a7d-aecb-c919eab0a88b
# ╟─9b928f14-b785-436e-8841-8b7bda8092d6
# ╟─153f2360-4777-4a0b-8d33-4051b905b514
# ╟─169cacd0-d9c6-46ac-be84-2ff0ac1f48dc
# ╟─6f1985e8-0976-4136-8061-58c3748580bd
# ╟─0c7bb53c-05fb-11ef-0fd7-b157332ae1ef
# ╟─2b170dab-8dd9-4e4a-b15a-ab51d056d954
# ╟─f2defed4-09d4-49d2-81a8-4c4e6b0d75f2
