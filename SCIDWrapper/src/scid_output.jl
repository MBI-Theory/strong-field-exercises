function load_calculation(dir)
    input_filename = joinpath(dir, "run.inp")
    isfile(input_filename) || return nothing
    input_params = load_input(input_filename)

    output_filename = joinpath(dir, "run.out")
    isfile(output_filename) || return nothing
    output = load_output(output_filename)

    detail_filename = joinpath(dir, output.calc_params.detail_output)
    detail = if isfile(detail_filename)
        load_detail(detail_filename, input_params, output.E₀)
    else
        (;)
    end

    volkov_filename = joinpath(dir, output.calc_params.sts_volkov_table)
    volkov_pes = if isfile(volkov_filename)
        load_pes(volkov_filename, output.calc_params)
    end

    coulomb_filename = joinpath(dir, output.calc_params.sts_coulomb_table)
    coulomb_pes = if isfile(coulomb_filename)
        load_pes(coulomb_filename, output.calc_params)
    end

    merge(detail, output, (volkov_pes=volkov_pes, coulomb_pes=coulomb_pes))
end

function load_detail(filename, input_params, E₀)
    data = readdlm(filename, comments=true)

    i = Int.(data[:,1])
    cplx(i,j) = data[:,i] + im*data[:,j]

    ω₀ = input_params.vp_param[1][2][1]
    A₀ = input_params.vp_scale
    Iₚ = -real(E₀)
    Uₚ = A₀^2/4
    cutoff = 3.17Uₚ + Iₚ
    # @info "Field params" ω₀ A₀ Iₚ Uₚ cutoff

    (i=i, t=data[:,2],
     vp_mag=data[:,3], vp_θ=data[:,3], vp_ϕ=data[:,3],
     norm=cplx(6,7), avg_h=cplx(8,9),
     x=cplx(12,13), y=cplx(14,15), z=cplx(16,17),
     ax=cplx(18,19), ay=cplx(20,21), az=cplx(22,23),
     vx=cplx(24,25), vy=cplx(26,27), vz=cplx(28,29),
     Ex=data[:,30], Ey=data[:,31], Ez=data[:,32],
     input_params=input_params,
     Iₚ=Iₚ, Uₚ=Uₚ, cutoff=cutoff)
end

function read_table(io::IO, header)
    data = readdlm(io)
    col_indices = UnitRange{Int}[]
    let i = 1
        for (_,T) in header
            if T <: Real
                push!(col_indices, i:i)
                i = i + 1
            elseif T <: Complex
                push!(col_indices, i:i+1)
                i = i + 2
            else
                @error "Do not know how to load type $(T)"
                i = i + 1
            end
        end
    end

    function mapT(T, d)
        m,n = size(d)
        [T(d[i,:]...) for i = 1:m]
    end

    DataFrame([name => mapT(T, data[:,is]) for ((name,T),is) in zip(header, col_indices)])
end

function load_output(filename)
    open(filename) do file
        function readlinesuntil(pred::Function)
            while !eof(file)
                line = readline(file)
                p = pred(line)
                isnothing(p) || return p
            end
            nothing
        end
        readlinesuntil(r::Regex) = readlinesuntil(Base.Fix1(match, r))
        readlinesuntil(s::AbstractString) = readlinesuntil(r""*s)

        function extract_table(header)
            readline(file)
            readline(file)
            buf = IOBuffer()
            while !eof(file)
                line = strip(readline(file))
                isempty(line) && break
                println(buf, line)
            end
            seek(buf, 0)
            read_table(buf, header)
        end

        readlinesuntil(r" ===== begin simulation parameters ===== ")
        calc_params = load_input(file)
        eng_line = readlinesuntil(r"Energy of the atomic solution is[ ]+"*FLOAT_PAT*r"[ ]*"*FLOAT_PAT*r"[ ]+Hartree")
        E₀ = if !isnothing(eng_line)
            parse(Float64, eng_line[1]) + im*parse(Float64, eng_line[2])
        end

        final_norm = if !isnothing(readlinesuntil("Final norm, by total angular momentum"))
            readline(file)
            extract_table([("L", Int), ("M", Int), ("norm", ComplexF64)])
        end

        final_population = if !isnothing(readlinesuntil("Final populations, by total angular momentum"))
            readline(file)
            extract_table([("L", Int), ("Total population", ComplexF64), ("Bound population", ComplexF64), ("Continuum population", ComplexF64)])
        end

        final_population_m_resolved = if !isnothing(readlinesuntil("Final populations, by total angular momentum and angular momentum projection"))
            readline(file)
            extract_table([("L", Int), ("M", Int), ("Total population", ComplexF64), ("Bound population", ComplexF64), ("Continuum population", ComplexF64)])
        end

        large_amplitudes = if !isnothing(readlinesuntil("Large amplitudes of individual field-free states"))
            readline(file)
            extract_table([("L", Int), ("M", Int), ("i", Int), ("E", ComplexF64), ("Wgt", ComplexF64), ("<I|W>", ComplexF64), ("<W|I>", ComplexF64)])
        end

        (calc_params=calc_params, E₀=E₀,
         final_norm=final_norm,
         final_population=final_population, final_population_m_resolved=final_population_m_resolved,
         large_amplitudes=large_amplitudes)
    end
end

function load_pes(filename, calc_params)
    nk = calc_params.sts_kgrid_count
    nθ = calc_params.sts_dgrid_ntheta
    nϕ = calc_params.sts_dgrid_nphi

    A = zeros(ComplexF64, nk, nθ, nϕ)

    k = zeros(Float64, nk)
    θ = zeros(Float64, nθ)
    ϕ = zeros(Float64, nϕ)

    open(filename) do file
        for _ = 1:7
            readline(file)
        end
        for iϕ = 1:nϕ
            for iθ = 1:nθ
                buf = IOBuffer()
                for ik = 1:nk
                    println(buf, strip(readline(file)))
                end
                readline(file)
                seek(buf, 0)
                data = readdlm(buf)
                if iϕ == iθ == 1
                    k .= data[:,4]
                end
                if iϕ == 1
                    θ[iθ] = data[1,5]
                end
                if iθ == 1
                    ϕ[iϕ] = data[1,6]
                end
                A[:,iθ,iϕ] .= data[:,7] .+ im .* data[:,8]
            end
        end
    end

    (A=A, k=k, θ=θ, ϕ=ϕ)
end

function load_scan(fun::Function, Is::CartesianIndices, scan_dir::AbstractString)
    scan_file = joinpath(scan_dir, "scan.csv")
    isfile(scan_file) && return DataFrame(CSV.File(scan_file))

    data = @withprogress name="Loading scan" begin
        N = length(Is)
        map(enumerate(Is)) do (i,I)
            run_dir = joinpath(scan_dir, string.(Tuple(I))...)
            data = load_calculation(run_dir)
            @logprogress i/N
            fun(I, data)::NamedTuple
        end
    end
    df = DataFrame(data)
    CSV.write(scan_file, df)
    df
end

load_scan(fun::Function, A::AbstractArray, args...; kwargs...) =
    load_scan((I,data) -> fun(A[I],data), CartesianIndices(A), args...; kwargs...)

load_scan(fun::Function, x::Base.Iterators.ProductIterator, args...; kwargs...) =
    load_scan(fun, collect(x), args...; kwargs...)
