abstract type AbstractPotential{T} end

include("table_of_elements.jl")

# Effective nuclear charge, subtracting those core electrons that are
# modelled by e.g. a pseudopotential.
effective_charge(p::AbstractPotential) = charge(p) - num_electrons(core(ground_state(p)))

isspinlocal(::AbstractPotential) = true

# * Point charge nucleus

struct PointCharge{T} <: AbstractPotential{T}
    Z::T
end

PointCharge(s::Symbol) = PointCharge(element_number(s))

Base.hash(p::PointCharge, h::UInt) = hash(p.Z, h)

macro pc_str(s)
    :(PointCharge(Symbol($s)))
end

function find_element(Z::Real)
    isinteger(Z) && (Z = Int(Z))
    Z ∈ eachindex(table_of_elements) || return nothing
    table_of_elements[Z]
end

function element_name(p::PointCharge)
    iszero(p.Z) && return "vacuum"
    element = find_element(p.Z)
    isnothing(element) ? "unknown" : element[2]
end

Base.show(io::IO, p::PointCharge{I}) where {I<:Integer} =
    write(io, "Z = $(p.Z) [$(element_name(p))]")

Base.show(io::IO, p::PointCharge) =
    write(io, "Z = $(p.Z)")

(p::PointCharge{T})( ::O, r::U) where {T,O,U} = -p.Z/r
(p::PointCharge{T})(o::O, r::VU) where {T,O,U,VU<:AbstractVector{U}} = p.(Ref(o), r)

charge(p::PointCharge) = p.Z
effective_charge(p::PointCharge) = charge(p)

islocal(::PointCharge) = true

Base.:(+)(p::PointCharge, q) =
    PointCharge(p.Z + q)

Base.:(-)(p::PointCharge, q) = p + (-q)

# * Yukawa

struct Yukawa{T} <: AbstractPotential{T}
    g::T
    α::T
    m::T
end

Base.hash(p::Yukawa, h::UInt) = hash(p.g, hash(p.α, hash(p.m, h)))

(p::Yukawa)( ::O, r::U) where {O,U} = -p.g^2*exp(-p.α*p.m*r)/r
(p::Yukawa)(o::O, r::VU) where {O,VU<:AbstractVector} = p.(Ref(o), r)

charge(p::Yukawa) = p.g^2
effective_charge(p::Yukawa) = charge(p)

islocal(::Yukawa) = true

Base.show(io::IO, p::Yukawa{T}) where T =
    printfmt(io, "Yukawa{{{1:s}}}({2:0.3f}, {3:0.3f}, {4:0.3f}): -{5:0.3f}*exp(-{6:0.3f}*r)/r",
             T, p.g, p.α, p.m,
             p.g^2, p.α*p.m)

# * Exports

export PointCharge, @pc_str, Yukawa
