### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ 84fe36e4-1416-11ef-0c00-bf1237e41321
begin
    using Pkg
    Pkg.activate(".")

    using CairoMakie

    using Unitful
    using UnitfulAtomic

    using ElectricFields

    using Revise

    using StrongFieldApproximation

    function get_limits!(ax)
        reset_limits!(ax)
        box = ax.finallimits[]
        x₀,y₀ = box.origin
        wx,wy = box.widths
        x₀ .+ [0, wx], y₀ .+ [0, wy]
    end
end;

# ╔═╡ 189cc73a-4db0-4a56-9721-a510f37f5ce8
λ = 800u"nm"

# ╔═╡ 8960fa73-89bd-4017-b996-01256b90f0cd
ω = 2π*u"ħ*c" ./ λ .|> u"eV"

# ╔═╡ a2b3a163-92c1-47ad-a234-0be791c86ec7
I = 10.0 .^ range(0, stop=3, length=101)*u"TW/cm^2"

# ╔═╡ 01646e64-6b65-4d35-b804-24e03d6a1b2c
Iau = ElectricFields.Iaustrip.(I)

# ╔═╡ 0896346c-ced5-41c8-970e-e4c8a7130bbd
Uₚ = Iau / (4austrip(ω)^2)

# ╔═╡ 27559cec-6d58-4e5b-890b-077241947366
F = .√(Iau)

# ╔═╡ a9dc6dc6-cf5d-43cb-bd70-5ddedb30d4eb
Iₚ = 0.5u"hartree"

# ╔═╡ 398e2f27-f99e-4325-bdfc-0d070eb6ab53
γ = .√(austrip(Iₚ) ./ (2Uₚ))

# ╔═╡ bb2a0db9-ccf2-4bc8-9ff3-e2a0e9d4edcb
w_DC = exp.(-2*(2austrip(Iₚ)^(3/2)) ./ (3F))

# ╔═╡ 7134ece8-f19a-4f57-9d9c-15c08cb9a3cc
w_MPI = (sqrt(exp(1)) ./ (2γ)) .^ (2Iₚ/ω)

# ╔═╡ 2a9cb070-d091-4e31-8742-aeb54cfce6fe
w_Keldysh = exp.(-2*austrip(Iₚ)/austrip(ω)*((1 .+ 1 ./ (2γ.^2)) .* asinh.(γ) .- .√(γ .^2 .+ 1)./ (2γ)))

# ╔═╡ 7c148e18-b317-424d-80f1-2a6d603065f4
w_ppt = StrongFieldApproximation.IonizationRates.PPT.(austrip(Iₚ), Iau, austrip(ω), 0, 0)

# ╔═╡ 90d53ff4-bab2-4738-9f45-4b87819e13d4
let
    fig = Figure(size=(900,600))

    a = ceil(Int, minimum(γ))
    b = ceil(Int, maximum(γ))
    s = ceil(Int,(b-a)/4)
    γt = a:s:b
    Uₚt = (austrip(Iₚ) ./ (γt .^ 2))/2
    Iaut = Uₚt * 4austrip(ω)^2
    It = (Iaut*ElectricFields.Iau/u"TW/cm^2") .|> NoUnits
    Ft = .√(Iaut)

    @info "Ticks" It Ft

    for (j,(x, xscale, xlabel, xt)) in enumerate(((ustrip.(I), log10, L"$I$ [TW/cm$^2$]", It),
                                                  (1 ./ F, identity, L"$1 / F$ [au]", 1 ./ Ft)))
        ax = Axis(fig[1,j], xscale=xscale, yscale=log10, xlabel=xlabel,
                  ylabel = L"\Gamma", ylabelvisible=j == 1,
                  yticksvisible=j == 1, yticklabelsvisible=j == 1)

        lines!(x, w_Keldysh * w_ppt[end]/w_Keldysh[end], label="Keldysh", linewidth=3.0)
        lines!(x, w_DC * w_ppt[end]/w_DC[end], label=L"DC: $|\exp[-2(2I_p)^{3/2}/3F]|^2$", linestyle=:dash)
        lines!(x, w_MPI * w_ppt[1]/w_MPI[1], label=L"MPI: $|\gamma^{-2I_p/\omega}|^2$", linestyle=:dash)
        lines!(x, w_ppt, label="PPT", color="black")

        j == 2 && axislegend(ax, position=:lb)

        xl = first(get_limits!(ax))

        ax2 = Axis(fig[1,j], xaxisposition=:top, xscale=xscale, xlabel=L"\gamma")
        hidespines!(ax2)
        hidexdecorations!(ax2, label=false, ticklabels=false, ticks=false)
        hideydecorations!(ax2)
        linkxaxes!(ax, ax2)
        ax2.xticks = (xt, string.(γt))
    end

    save("../../../figures/sfi.pdf", fig)
    fig
end

# ╔═╡ 405ca664-288d-4efd-a57b-856d90314f63
let
        x = range(-1,stop=1,length=101)

        y = 2asinh.(x)

        fig = Figure()
        ax = Axis(fig[1,1], yscale=log10)
        lines!(x, abs.(sinh.(x) .- asin.(x)) .+ 1e-20)
        fig
end

# ╔═╡ Cell order:
# ╠═84fe36e4-1416-11ef-0c00-bf1237e41321
# ╠═189cc73a-4db0-4a56-9721-a510f37f5ce8
# ╠═8960fa73-89bd-4017-b996-01256b90f0cd
# ╠═a2b3a163-92c1-47ad-a234-0be791c86ec7
# ╠═01646e64-6b65-4d35-b804-24e03d6a1b2c
# ╠═0896346c-ced5-41c8-970e-e4c8a7130bbd
# ╠═27559cec-6d58-4e5b-890b-077241947366
# ╠═a9dc6dc6-cf5d-43cb-bd70-5ddedb30d4eb
# ╠═398e2f27-f99e-4325-bdfc-0d070eb6ab53
# ╠═bb2a0db9-ccf2-4bc8-9ff3-e2a0e9d4edcb
# ╠═7134ece8-f19a-4f57-9d9c-15c08cb9a3cc
# ╠═2a9cb070-d091-4e31-8742-aeb54cfce6fe
# ╠═7c148e18-b317-424d-80f1-2a6d603065f4
# ╠═90d53ff4-bab2-4738-9f45-4b87819e13d4
# ╠═405ca664-288d-4efd-a57b-856d90314f63
