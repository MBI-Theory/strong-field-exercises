function run_calc(sys::System, SCID_EXE::AbstractString;
                  run_dir="",
                  verbosity=0, overwrite=false, debug_plots=false,
                  generate_input_only=false,
                  omp_stacksize="65520k", kwargs...)
    if isempty(run_dir)
        run_dir = joinpath("runs", string(hash(sys)))
    end
    run_dir = abspath(run_dir)
    verbosity > 0 && @info "Will run calculation in subdirectory" run_dir
    mkpath(run_dir)

    input_filename = joinpath(run_dir, "run.inp")
    if !isfile(input_filename) || overwrite
        generate_scid_input(input_filename, sys;
            verbosity=verbosity-1, kwargs...)
        verbosity > 1 && @info "Input file created" input_filename
    end
    if generate_input_only
        @info "Input file generated" input_filename
        input_params = load_input(input_filename)
        @info "Parameters" input_params
        return run_dir
    end

    verbosity > 1 && @info "SCID command" SCID_EXE

    # SCID_EXE = abspath(SCID_EXE)
    # if !isfile(SCID_EXE)
    #     throw(ArgumentError("SCID executable path provided $(SCID_EXE) not present"))
    # end

    input_params = load_input(input_filename)
    verbosity > 10 && @info "Input parameters" input_params
    :timesteps ∈ keys(input_params) ||
    throw(ArgumentError("timesteps missing from input parameters"))
    timesteps = input_params.timesteps

    au2fs = auconvert(u"fs", 1)
    # ElectricFields.jl centres the pulse around zero, so we need to
    # subtract the origin in the plots below.
    t₀ = input_params.vp_param[1][2][3]

    if debug_plots
        F = sys.F
        t = timeaxis(F) .+ t₀
        Fv = field_amplitude(F, t .- t₀)
        Av = vector_potential(F, t .- t₀)
        tplot = au2fs*t

        pF = plot(tplot, Fv, ylabel=L"$F(t)$ [au]")
        pA = plot(tplot, Av, ylabel=L"$A(t)$ [au]")
        p = plot(pF, pA, layout=@layout([a;b]), xlabel=L"$t$ [fs]")

        @info "Proposed driving field" p
    end

    ENV["OMP_STACKSIZE"] = omp_stacksize

    cd(run_dir) do
        output_filename = joinpath(run_dir, "run.out")
        err_filename = joinpath(run_dir, "run.err")
        if !isfile(output_filename) || overwrite
            cmd = `$(split(SCID_EXE))`
            verbosity > 1 && @info "Calculation will commence" cmd output_filename

            # Remove old photoelectron spectra, if present, since if
            # this run does not request them, the old files may have
            # data of the wrong dimensions, leading to errors when
            # trying to load them.
            for filename in [:sts_volkov_table, :sts_coulomb_table, :sts_coulomb_waves]
                rm(joinpath(run_dir, input_params[filename]), force=true)
            end

            in_file = open(input_filename, "r")
            out_file = open(output_filename, "w")
            err_file = open(err_filename, "w")

            # https://discourse.julialang.org/t/avoiding-readavailable-when-communicating-with-long-lived-external-program/61611/21
            pout = Base.PipeEndpoint()
            proc = run(cmd, in_file, pout, err_file, wait = false)
            process_running(proc) || error("Could not start SCID: $(SCID_EXE)")

            progress_pat = r"^@[ ]*([0-9]+) "
            time_prop_finished_pat = r"Final left/right nradial"
            eigen_pat = r"Analyzing the final wavefunction in terms of field-free eigenstates"
            starting_time = now()
            try
                @withprogress name="Computing" begin
                    while !eof(pout) && process_running(proc)
                        line = readline(pout)
                        println(out_file, line)
                        m = match(progress_pat, line)
                        if !isnothing(m)
                            i = parse(Int, m[1])
                            @logprogress i/timesteps
                        elseif !isnothing(match(time_prop_finished_pat, line))
                            verbosity > 0 && @info "Time propagation finished in $(now()-starting_time)"
                        elseif !isnothing(match(eigen_pat, line))
                            verbosity > 0 && @info "Analyzing wavefunction in terms of eigenstates, hang tight"
                        end
                    end
                end
            catch e
                if e isa ProcessFailedException
                    println("External process failed:")
                    println(e)
                else
                    rethrow()
                end
            finally
                close(pout)
                close(in_file)
                close(out_file)
                close(err_file)
            end
            err = strip(read(err_filename, String))
            if !isempty(err)
                @error "SCID crashed:\n\n$(err)"
                return nothing
            end
            verbosity > 0 && @info "Calculation finished in $(now()-starting_time)"
        else
            verbosity > 0 && @warn "Calculation already exists, loading instead" output_filename
        end
    end

    all_runs_dir = abspath(joinpath(run_dir, ".."))
    cd(all_runs_dir) do
        du = read(`du -h`, String)
        verbosity > 0 && @info "Disk space usage:\n$(du)" all_runs_dir
    end

    if debug_plots
        data = load_calculation(run_dir)

        # Hack to please Plots.jl that refuses to deal with
        # non-positive values for logarithmix axes.
        nan_map(v) = v ≤ 0 ? NaN : v

        F = sys.F
        t = data.t
        Av = vector_potential(F, t .- t₀)
        tplot = au2fs*t

        ΔA = abs.(data.vp_mag-Av)

        pA = plot(tplot, [data.vp_mag Av],
                  label=[L"$A_{\mathrm{calc}}(t)$" L"$A_{\mathrm{prop}}(t)$"],
                  ylabel=L"$A(t)$ [au]")
        pΔA = plot(tplot, nan_map.(ΔA), ylabel=L"$\Delta A(t)$ [au]",
                   yscale=:log10, ylim=(1e-20, 10maximum(ΔA)))
        p = plot(pA, pΔA, layout=@layout([a;b]), xlabel=L"$t$ [fs]")
        @info "Actual vector potential used" p

        p = plot(tplot, [nan_map.(abs.(1.0 .- real(data.norm))) nan_map.(abs.(imag(data.norm)))],
            label=[L"|1 - \Re\{N\}|" L"|\Im\{N\}|"], xlabel=L"$t$ [fs]",
            yscale=:log10, ylim=(1e-14, 10))
        @info "Lost norm" p
    end

    run_dir
end

function run_scan(fun::Function, Is::CartesianIndices,
                  SCID_EXE::AbstractString,
                  scan_dir::AbstractString; kwargs...)
    @progress for I in Is
        run_dir = joinpath(scan_dir, string.(Tuple(I))...)
        mkpath(run_dir)
        sys = fun(I)
        # By default, we avoid detail output; the user may override
        # this.
        run_calc(sys, SCID_EXE, run_dir; detail_output=" ", kwargs...)
    end
    scan_dir
end

run_scan(fun::Function, A::AbstractArray, args...; kwargs...) =
    run_scan(I -> fun(A[I]), CartesianIndices(A), args...; kwargs...)

run_scan(fun::Function, x::Base.Iterators.ProductIterator, args...; kwargs...) =
    run_scan(fun, collect(x), args...; kwargs...)
