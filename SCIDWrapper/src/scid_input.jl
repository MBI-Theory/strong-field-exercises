const DEFAULT_OPTIONS =
          (comment                = "Hydrogen 1S HHG, linear 800 nm 1e14 W/cm2, uniform grid",
           verbose                = 1,
           omp_num_threads        = 8,
           initial_wfn            = "single",
           initial_wfn_index      = (0, 0, 1),
           initial_wfn_energy     = -0.500+0.0im,
           sd_nradial             = 300,
           sd_lmax                =  20,
           sd_mmin                =   0,
           sd_mmax                =   0,
           field_unwrap           = true,
           rotation_mode          = "auto",
           sd_rgrid               = "uniform",
           sd_rgrid_zeta          = 1.0,
           sd_rgrid_dr            = 0.20,
           pot_name               = "hydrogenic",
           pot_param              = 1.0,
           task                   = "real time",
           dt                     = 0.01,
           timesteps              = 49000,
           dt_subdivision         = "off",
           vp_shape               = "z Gaussian",
           vp_scale               = 0.938582,
           vp_param               = [1:4 => (0.05695, 0.000000, 250.000,  200.0), 11:12 => (170.0, 245.0)],
           skip_tests             = true,
           output_each            = 200,
           composition_threshold  = 1e-5,
           initial_wf_dump_prefix = " ",
           final_wf_dump_prefix   = " ",
           field_preview          = " ",
           wt_atomic_cache_prefix = " ",
           detail_output          = "output.detail",
           detail_frequency       = 100,
           do_dipole_plasma       = true,
           sts_volkov_table       = "volkov.pes",
           sts_coulomb_table      = "coulomb.pes",
           sts_coulomb_waves      = "coulomb.waves",
           sts_verbose            = 0)

format_value(value::Tuple) = join(map(format_value, value), ", ")
format_value(value::AbstractString) = format("'{1:s}'", value)
format_value(value::Bool) = format(".{1:s}.", value ? "true" : "false")
format_value(value::Integer) = format("{1:s}", value)
format_value(value::AbstractFloat) = format("{1:.17g}", value)
format_value(value::Complex) = format("({1:s},{2:s})", format_value(real(value)), format_value(imag(value)))

format_param(param::Symbol, value::Tuple) = [(string(param), format_value(value))]

format_param(param::Symbol, value::Vector{<:Pair{<:UnitRange}}) =
    [(format("{1:s}({2:d}:{3:d})", param, first(k), last(k)), format_value(v))
    for (k,v) in value]

format_param(param::Symbol, value) =
    [(string(param), format_value(value))]

function generate_scid_input(params::NamedTuple)
    params_strs = Vector{Tuple{String,String}}()
    for (k,v) in pairs(params)
        s = format_param(k, v)
        append!(params_strs, s)
    end

    n = maximum(length, first.(params_strs))
    fmt = FormatExpr("   {1:$(n)s} = {2:s}")
    params_str = map(((k,v),) -> format(fmt, k, v), params_strs)

    " &sph_tdse\n"*join(params_str, ",\n")*"\n /"
end

function scid_params(F::ElectricFields.LinearField, δt;
                     auto_decrease_timestep=true, default_sampling_factor=10,
                     δtobs=1.0, verbosity=0, kwargs...)
    env = envelope(F)
    env isa ElectricFields.TruncatedGaussianEnvelope ||
    throw(ArgumentError("Only ElectricFields.TruncatedGaussianEnvelope supported for now"))

    δt > 0 ||
    throw(ArgumentError("Invalid time step δt = $(δt), must be larger than zero"))

    A₀ = vector_potential(F)
    ω = photon_energy(F)
    # ElectricFields.jl uses a sine carrier for the vector potential,
    # SCID a cosine one.
    ϕ = Float64(carrier(F).ϕ) - π/2

    s = span(F)
    a,b = s.left,s.right
    origin = (b-a)/2

    fs = 1/δt
    timesteps = steps(F, fs)
    max_fs = max(fs, 1/0.01, default_sampling_factor*austrip(max_frequency(F)))
    δtn = step(timeaxis(F, max_fs))
    if δtn < δt && auto_decrease_timestep
        timesteps = steps(F, max_fs)
        verbosity > 0 && @info "Decreased time step from $(δt) to $(δtn)" timesteps
        δt = δtn
    end
    τ = duration(F)

    δtobs = max(δtobs, δt)
    detail_frequency = ceil(Int, δtobs/δt)

    (dt=δt, timesteps=timesteps,
     vp_scale=A₀,
     vp_param=[1:4 => (ω, ϕ, origin, τ), 11:12 => (env.toff, env.tmax)],
     detail_frequency=detail_frequency)
end

function scid_params(p::PointCharge; kwargs...)
    Z = Float64(charge(p))
    (pot_name="hydrogenic", pot_param=Z, sd_rgrid_zeta=Z)
end

function scid_params(p::Yukawa; kwargs...)
    Z = Float64(charge(p))
    A = Float64(p.α*p.m)
    (pot_name="yukawa", pot_param=[1:2 => (Z,A)], sd_rgrid_zeta=Z)
end

function scid_params(; rmax=60, dr=0.2, ℓmax::Integer=20,
                     cap_kind="manolopoulos", cap_kmin=0.2, cap_delta=0.2,
                     composition_threshold=1e-5,
                     volkov_tsurff=false, volkov_isurff=false, coulomb_isurff=false,
                     surff_kmax=3.0, surff_nk=100, surff_nθ=15, surff_nϕ=30,
                     omp_num_threads=8,
                     kwargs...)
    rmax > 0 || throw(ArgumentError("Radial box size must be positive"))
    dr > 0 || throw(ArgumentError("Radial step must be positive"))
    ℓmax > 0 || throw(ArgumentError("ℓmax must be positive"))

    cap_kind = lowercase(cap_kind)
    cap_kind == "manolopoulos" || cap_kind == "none" ||
    throw(ArgumentError("Unknown CAP kind $(cap_kind)"))

    (sd_nradial=ceil(Int, rmax/dr), sd_rgrid_dr=dr, sd_lmax=ℓmax,

     cap_name=cap_kind, cap_param=[1:2 => (cap_kmin, cap_delta)],

     composition_threshold=composition_threshold,

     sts_volkov=volkov_tsurff, sts_volkov_atend=volkov_isurff,
     sts_coulomb_atend=coulomb_isurff,
     sts_kgrid_max=surff_kmax, sts_kgrid_count=surff_nk,
     sts_dgrid_ntheta=surff_nθ, sts_dgrid_nphi=surff_nϕ,

     omp_num_threads=omp_num_threads)
end

function generate_scid_input(sys::System; kwargs...)
    params = merge(DEFAULT_OPTIONS, scid_params(sys.F, sys.δt; kwargs...), scid_params(sys.potential; kwargs...),
                   scid_params(;kwargs...))
    generate_scid_input(params)
end

function generate_scid_input(filename::AbstractString, args...; kwargs...)
    if isdir(filename)
        filename = joinpath(filename, "run.inp")
        @info "Input will be saved to" filename
    end
    open(filename, "w") do file
        println(file, generate_scid_input(args...; kwargs...))
    end
end

const INT_PAT = int_pat = r"([+-]{0,1}[0-9]+)"
const FLOAT_PAT = r"([+-]{0,1}[0-9]+\.[0-9]*(?:[eE][+-]{0,1}[0-9]+){0,1})"
const FLOAT_OR_INT_PAT = r"((?:[+-]{0,1}[0-9]+\.[0-9]*(?:[eE][+-]{0,1}[0-9]+){0,1})|(?:[+-]{0,1}[0-9]+))"
const COMPLEX_PAT = r"\([ ]*"*FLOAT_OR_INT_PAT*r"[ ]*,[ ]*"*FLOAT_OR_INT_PAT*r"[ ]*\)"

function parse_value(v::AbstractString)
    m = match(r"^'(.*)'$", v)
    isnothing(m) || return String(strip(m[1]))
    m = match(r"^\"(.*)\"$", v)
    isnothing(m) || return String(strip(m[1]))

    isnothing(match(r"^\.true\.$", v)) || return true
    isnothing(match(r"^T$", v)) || return true
    isnothing(match(r"^\.false\.$", v)) || return false
    isnothing(match(r"^F$", v)) || return false

    m = match(r"^"*INT_PAT*r"$", v)
    isnothing(m) || return parse(Int, m[1])

    m = match(r"^"*FLOAT_PAT*r"$", v)
    isnothing(m) || return parse(Float64, m[1])

    m = match(r"^"*COMPLEX_PAT*r"$", v)
    isnothing(m) || return complex(parse_value(m[1]), parse_value(m[2]))

    # This will obviously not work with vectors/tuples of complex
    # numbers, since the string will be split in the wrong places.
    if ',' ∈ v || '*' ∈ v
        els = map(strip, split(v, ','))
        vv = map(els) do e
            if '*' ∈ e
                mult,ee = split(e, '*')
                fill(parse_value(ee), parse(Int, mult))
            else
                [parse_value(e)]
            end
        end
        return reduce(vcat, vv)
    end

    @warn "Could not parse value, returning string" v

    v
end

function load_input(io::IO)
    params = (;)
    pstart = r"^[ ]*&(?:sph_tdse)|(?:SPH_TDSE)[ ]*$"
    if isnothing(match(pstart, readline(io)))
        @error "Expected start of namelist &sph_tdse"
        return params
    end
    p = r"^[ ]*([^ ]+?)(?:\(([0-9]+):([0-9]+)\)){0,1}[ ]*=[ ]*(.+?),{0,1}$"
    pend = r"^[ ]*/[ ]*$"
    for line in eachline(io)
        isnothing(match(pend, line)) || break
        m = match(p, line)
        isnothing(m) && continue
        k = Symbol(lowercase(m[1]))
        v = parse_value(strip(m[4]))
        if isnothing(m[2])
            params = merge(params, (; k => v))
        else
            interval = parse(Int,m[2]):parse(Int,m[3])
            new_entry = [interval => v]
            v = if k ∈ keys(params)
                vcat(params[k], new_entry)
            else
                new_entry
            end
            params = merge(params, (; k => v))
        end
    end
    params
end

load_input(filename::AbstractString) = open(load_input, filename)
