plot_pes(::Nothing; kwargs...) = nothing
function plot_pes(pes; Uₚ=nothing, dynamic_range=3, projection=:box)
    pp = plot()

    nd = ndims(pes.A)
    A = if nd == 2
        pes.A
    else
        ∑A = sum(pes.A, dims=3:nd)
        for i = nd:-1:3
            ∑A = dropdims(∑A, dims=i)
        end
        ∑A
    end

    Z = log10.(abs2.(A))
    Zmin,Zmax = extrema(Z)

    if projection == :box
        heatmap!(pp, pes.θ, pes.k, Z, clims=(Zmax-dynamic_range,Zmax),
            xlabel=L"$\theta$ [$^\circ$]", ylabel=L"$|k|$ [Bohr$^{-1}$]")
    elseif projection == :polar
        θ = vcat(pes.θ, reverse(360 .- pes.θ))
        Z = hcat(Z,Z[:,end:-1:1])
        heatmap!(pp, deg2rad.(θ), pes.k, Z, clims=(Zmax-dynamic_range,Zmax), projection=:polar)
    else
        @error "Unknown projection $(projection)"
    end

    pl = plot()
    for θ in [0,90,180,rad2deg(atan(√2))]
        iθ = argmin(i -> abs(pes.θ[i]-θ), eachindex(pes.θ))
        plot!(pl, pes.k, abs2.(A[:,iθ]), yaxis=:log10, label=L"\theta = %$(round_quantity(θ))^\circ",
              xlabel=L"$|k|$ [Bohr$^{-1}$]", ylabel=L"$P$ [Ha$^{-1}$]",
              legend=:bottomleft)
    end

    if !isnothing(Uₚ)
        yl = ylims(pp)
        ns = [2,10]
        marks = filter(<(yl[2]), .√(2austrip(Uₚ)*ns))
        projection == :box && hline!(pp, marks, color="white", label=nothing)
        for (m,ls,n) in zip(marks,[:solid,:dash],ns)
            vline!(pl, [m], color="black", linestyle=ls, label=L"%$(n)U_p")
        end
    end

    plot(pp, pl, size=(900, 1200), layout=@layout([a; b]))
end
