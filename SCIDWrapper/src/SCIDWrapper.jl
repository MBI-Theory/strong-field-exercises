module SCIDWrapper

using ElectricFields
using Unitful
using UnitfulAtomic

using AbstractFFTs
using FFTW
import DSP: hanning
using Statistics

using Format
using DelimitedFiles
using DataFrames
using CSV

using Dates
using ProgressLogging

using Plots
using LaTeXStrings

include("potentials.jl")

struct System{Field<:ElectricFields.AbstractField,T,
              Potential<:AbstractPotential}
    F::Field
    δt::T
    potential::Potential
end

# Base.hash(sys::System, h::UInt) = hash(sys.F, hash(sys.δt, hash(sys.potential, h)))
Base.hash(sys::System, h::UInt) = hash(sys.δt, hash(sys.potential, h))

include("scid_input.jl")
include("scid_output.jl")
include("run.jl")

include("dipole_analysis.jl")
include("population_analysis.jl")
include("pes_analysis.jl")

end
