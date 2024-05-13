#!/usr/bin/env julia

function print_help_message()
    println("Usage: $(@__FILE__) [options]")
    println()
    println("Valid options are:")
    println("  -(n)lb : Do(n't) launch browser and navigate to the Pluto.jl notebook interface (default: true)")
    println("  -(n)ar : Do(n't) auto-reload notebook upon changes to the underlying Julia file (default: false)")
    println("  -(n)i  : Instantiate Julia environment (default: false)")
    println("  -h     : Print this help message")
end

function consume_negatable_flag(positive_flag, default)
    # We only allow -f or -nf (or neither, but not both).
    pf = "-"*positive_flag ∈ ARGS
    nf = "-n"*positive_flag ∈ ARGS
    if pf && nf
        error("-$(positive_flag) and -n$(positive_flag) are mutually exclusive")
    end
    deleteat!(ARGS, findall(==("-"*positive_flag), ARGS))
    deleteat!(ARGS, findall(==("-n"*positive_flag), ARGS))
    (pf | nf) ? pf : default
end

"SCID_EXE" ∈ keys(ENV) ||
throw(ArgumentError("You need to set the path to SCID in the environment variable SCID_EXE"))

if "-h" ∈ ARGS
    print_help_message()
    exit(0)
end

launch_browser = consume_negatable_flag("lb", true)
auto_reload_from_file = consume_negatable_flag("ar", false)
perform_instantiation = consume_negatable_flag("i", false)

if !isempty(ARGS)
    println("Invalid arguments: $(join(ARGS, ' '))")
    println()
    print_help_message()
    exit(1)
end

using Pkg
Pkg.activate(".")
perform_instantiation && Pkg.instantiate()

using Pluto
Pluto.run(host="127.0.0.1", port=1234, launch_browser=launch_browser, auto_reload_from_file=auto_reload_from_file)
