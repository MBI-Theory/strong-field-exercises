plot_pes!(fig, ::Nothing; kwargs...) = fig
function plot_pes!(fig, pes; Uₚ=nothing, dynamic_range=3, projection=:box, title=nothing,
                   nolegend=false)

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

    θ = pes.θ

    At,akw = if projection == :box
        Axis, (;xlabel=L"$\theta$ [$^\circ$]", ylabel=L"$|k|$ [Bohr$^{-1}$]")
    elseif projection == :polar
        θ = deg2rad.(vcat(pes.θ, reverse(360 .- pes.θ)))
        Z = hcat(Z,Z[:,end:-1:1])
        if θ[1] ≠ 0
            θ = vcat(θ, θ[1])
            Z = hcat(Z,Z[:,1])
        end
        PolarAxis, (;)
    else
        @error "Unknown projection $(projection)"
    end
    ax1 = At(fig[1,1]; title=title, akw...)

    pkw = (;colormap=:inferno, colorrange=(Zmax-dynamic_range,Zmax))
    if projection == :box
        heatmap!(ax1, θ, pes.k, Z'; pkw...)
    elseif projection == :polar
        surface!(ax1, θ, pes.k, zeros(size(Z')); color=Z',
            shading=NoShading, pkw...)
        hiderdecorations!(ax1)
    end

    ax2 = Axis(fig[2,1], xlabel=L"$|k|$ [Bohr$^{-1}$]", ylabel=L"$P$ [Ha$^{-1}$]",
               yscale=log10)
    for θ in [0,90,180,rad2deg(atan(√2))]
        iθ = argmin(i -> abs(pes.θ[i]-θ), eachindex(pes.θ))
        lines!(ax2, pes.k, abs2.(A[:,iθ]), label=L"\theta = %$(round_quantity(θ))^\circ")
    end

    if !isnothing(Uₚ)
        xl = get_limits!(ax2)[1]
        ns = [2,10]
        marks = filter(<(xl[2]), .√(2austrip(Uₚ)*ns))
        projection == :box && hlines!(ax1, marks, color="white", label=nothing)
        for (m,ls,n) in zip(marks,[:solid,:dash],ns)
            vlines!(ax2, m, color="black", linestyle=ls, label=L"%$(n)U_p")
        end
    end
    nolegend || axislegend(ax2, position=:lb)

    fig
end
plot_pes(args...; kwargs...) =
    plot_pes!(Figure(size=(900, 1200)), args...; kwargs...)
