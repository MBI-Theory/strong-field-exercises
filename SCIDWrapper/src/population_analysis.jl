plot_angular_decomposition(::Nothing; kwargs...) = nothing
function plot_angular_decomposition(results; ℓmax=results.input_params.sd_lmax)
    fp = results.final_population_m_resolved
    isnothing(fp) && return

    ℓs = fp.L
    ms = fp.M
    mmin,mmax = extrema(ms)

    pop_kinds = ["Total population", "Bound population", "Continuum population"]
    markershapes = [:none, :circle, :cross]

    pa = plot()

    for m ∈ mmin:mmax
        pm = findall(==(m), ms)
        isempty(pm) && continue

        for (pop_kind, markershape) in zip(pop_kinds, markershapes)
            w = real(fp[!, pop_kind])
            p = pm[findall(>(0), w[pm]) ∩ findall(<=(ℓmax), ℓs[pm])]
            plot!(pa, ℓs[p], w[p], label=L"%$(pop_kind), $m_\ell = %$(m)$",
                  markershape=markershape, linestyle=:auto)
        end
    end

    plot!(pa,
          title="Angular momentum-resolved population",
          xlabel=L"\ell", ylabel="Population",
          yaxis=:log10, size=(700,600))
end

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

function plot_eigen_decomposition!(pa, results; predicate::Union{Function,Nothing}=nothing,
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

    for ℓ ∈ 0:ℓmax
        for m ∈ -ℓ:ℓ
            isnothing(predicate) || predicate(ℓ, m) || continue
            pℓ = findall(==(ℓ), ℓs)
            pm = findall(==(m), ms)
            p = pℓ ∩ pm ∩ pE
            isempty(p) && continue
            p = p[findall(>(0), ws[p])]
            x,y = fill_missing(is[p], ws[p], NaN)
            plot!(pa, x .+ ℓ, y;
                  label=L"$\ell = %$(ℓ)$, $m_\ell = %$(m)$",
                  markershape=:circle, markersize=3.0,
                  markerstrokewidth=0.0,
                  kwargs...)
        end
    end

    plot!(pa,
          title="Eigenstate-resolved population",
          xlabel=L"# + $\ell$", ylabel="Population",
          yaxis=:log10, size=(900,600),
          legend=:topright)
end

plot_eigen_decomposition(::Nothing; kwargs...) = nothing
plot_eigen_decomposition(results; kwargs...) =
    plot_eigen_decomposition!(plot(), results; kwargs...)

plot_photon_diagram(::Nothing; kwargs...) = nothing
function plot_photon_diagram(results; draw_photons=true)
    la = results.large_amplitudes
    isnothing(la) && return

    Es = real(la.E)
    ℓs = la.L
    ℓmax = maximum(ℓs)
    ms = la.M

    ℓ₀,m₀,i₀ = results.calc_params.initial_wfn_index
    E₀ = real(results.E₀)
    ω₀ = results.input_params.vp_param[1][2][1]

    pd = plot()

    for (E,ℓ) in zip(Es, ℓs)
        E > 0 && continue
        plot!(pd, ℓ .+ [0, 0.5], E*[1,1], linewidth=2.0, label=nothing)
    end

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
                plot!(pd, [a[1], a[1]+v[1]], [a[2], a[2]+v[2]],
                    arrow=true,
                    color=:black, linewidth=0.5,
                    label=nothing)
                push!(photon_stack, (n+1,ℓn))
            end
        end
    end

    hline!(pd, [0.0], color=:black, linewidth=3.0, label=nothing)

    plot!(pd, xlabel=L"\ell", ylabel=L"$E$ [Ha]", size=(700,600))
end
