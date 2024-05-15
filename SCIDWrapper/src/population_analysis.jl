plot_angular_decomposition!(fig, ::Nothing; kwargs...) = fig
function plot_angular_decomposition!(fig, results; ℓmax=results.input_params.sd_lmax)
    fp = results.final_population_m_resolved
    isnothing(fp) && return

    ℓs = fp.L
    ms = fp.M
    mmin,mmax = extrema(ms)

    pop_kinds = ["Total population", "Bound population", "Continuum population"]
    marker = [:none, :circle, :cross]

    ax = Axis(fig[1,1],
              title="Angular momentum-resolved population",
              xlabel=L"\ell", ylabel="Population",
              yscale=log10)

    for m ∈ mmin:mmax
        pm = findall(==(m), ms)
        isempty(pm) && continue

        for pop_kind in pop_kinds
            w = real(fp[!, pop_kind])
            p = pm[findall(>(0), w[pm]) ∩ findall(<=(ℓmax), ℓs[pm])]
            scatterlines!(ax, ℓs[p], w[p], label=L"%$(pop_kind), $m_\ell = %$(m)$")
        end
    end
    axislegend(ax)

    fig
end
plot_angular_decomposition(args...; kwargs...) =
    plot_angular_decomposition!(Figure(size=(700,600)), args...; kwargs...)

function fill_missing(is, vs, f)
    nis = Vector{eltype(is)}()
    nvs = Vector{eltype(vs)}()

    iexp = 1
    for (ii,i) in enumerate(is)
        if i == iexp
            push!(nis, i)
            push!(nvs, vs[ii])
        else
            miss = iexp:i-1
            append!(nis, miss)
            append!(nvs, fill(f, length(miss)))
        end
        iexp = i + 1
    end

    nis,nvs
end

function plot_eigen_decomposition!(fig, results; predicate::Union{Function,Nothing}=nothing,
                                   Epredicate::Function=(E -> true), kwargs...)
    la = results.large_amplitudes
    isnothing(la) && return

    ℓmax = results.input_params.sd_lmax
    Es = real(la.E)
    ℓs = la.L
    ms = la.M
    is = la.i
    ws = real(la.Wgt)

    pE = findall(Epredicate, Es)

    ax = Axis(fig[1,1],
              title="Eigenstate-resolved population",
              xlabel=L"# + $\ell$", ylabel="Population",
              yscale=log10)

    for ℓ ∈ 0:ℓmax
        for m ∈ -ℓ:ℓ
            isnothing(predicate) || predicate(ℓ, m) || continue
            pℓ = findall(==(ℓ), ℓs)
            pm = findall(==(m), ms)
            p = pℓ ∩ pm ∩ pE
            isempty(p) && continue
            p = p[findall(>(0), ws[p])]
            x,y = fill_missing(is[p], ws[p], NaN)
            scatterlines!(ax, x .+ ℓ, y;
                          label=L"$\ell = %$(ℓ)$, $m_\ell = %$(m)$",
                          marker=:circle,
                          kwargs...)
        end
    end
    axislegend(ax)

    fig
end

plot_eigen_decomposition!(fig, ::Nothing; kwargs...) = fig
plot_eigen_decomposition(args...; kwargs...) =
    plot_eigen_decomposition!(Figure(size=(900,600)), args...; kwargs...)

plot_photon_diagram!(fig, ::Nothing; kwargs...) = fig
function plot_photon_diagram!(fig, results; draw_photons=true)
    la = results.large_amplitudes
    isnothing(la) && return

    Es = real(la.E)
    ℓs = la.L
    ℓmax = maximum(ℓs)
    ms = la.M

    ℓ₀,m₀,i₀ = results.calc_params.initial_wfn_index
    E₀ = real(results.E₀)
    ω₀ = results.input_params.vp_param[1][2][1]

    ax = Axis(fig[1,1], xlabel=L"\ell", ylabel=L"$E$ [Ha]")

    for (E,ℓ) in zip(Es, ℓs)
        E > 0 && continue
        lines!(ax, ℓ .+ [0, 0.5], E*[1,1], linewidth=2.0, label=nothing)
    end

    xs = Float64[]
    ys = Float64[]
    us = Float64[]
    vs = Float64[]

    if draw_photons
        photon_stack = Tuple{Int,Int}[]
        push!(photon_stack, (0,ℓ₀))
        max_photon = ceil(Int, (0 - E₀)/ω₀) + 1

        η = 0.95
        while !isempty(photon_stack)
            n,ℓ = pop!(photon_stack)
            n > max_photon && continue
            for Δℓ = (-1,1)
                ℓn = ℓ+Δℓ
                ℓn < 0 && continue
                ℓn == ℓmax && continue
                a = [ℓ+0.25, E₀+ω₀*n]
                b = [ℓn+0.25, E₀+ω₀*(n+1)]
                v = η*(b-a)
                push!(xs, a[1])
                push!(ys, a[2])
                push!(us, v[1])
                push!(vs, v[2])
                push!(photon_stack, (n+1,ℓn))
            end
        end
    end

    arrows!(ax,
            xs, ys, us, vs,
            color=:black, linewidth=1.0)

    hlines!(ax, [0.0], color=:black, linewidth=3.0, label=nothing)


    fig
end
plot_photon_diagram(args...; kwargs...) =
    plot_photon_diagram!(Figure(size=(700,600)), args...; kwargs...)
