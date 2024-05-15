module NotebookHelper

using DelimitedFiles

using Markdown
using InteractiveUtils

using PlutoUI
import PlutoUI: combine
using HypertextLiteral

using IntervalSets
using Statistics

using Unitful
using UnitfulAtomic
using ElectricFields

using JSON

export cached_calculation, notebook_styling,
    NumberIntervalField,
    Input, InputSection, format_input,
    run_dir_input, overwrite_input,
    presets_input, select_previous_runs,
    load_inputs, saved_inputs, save_inputs,
    get_electric_field

function cached_calculation(fun::Function, ::Type{T}, filename::AbstractString; overwrite::Bool=false) where  T
    if isfile(filename) && !overwrite
        @info "Loading cached results" filename
        return readdlm(filename,',',T,'\n')
    end
    res = fun()
    writedlm(filename, res, ',')
    res
end

function notebook_styling()
    # https://discourse.julialang.org/t/cell-width-in-pluto-notebook/49761/11
    html"""
<style>
  @media screen {
      main {
          margin: 0 auto;
          max-width: 2000px;
          padding-left: max(150px, 10%);
          padding-right: max(150px, 10%);
      }
      pluto-logs-container {
          max-height: unset;
      }
  }
</style>
"""
end


function closest(i::Interval{L,R,<:AbstractFloat}, v) where {L,R}
    l,r = i.left,i.right
    if v ∈ i
        v
    elseif v < l
        l + (L == :closed ? 0 : eps(l))
    elseif v > r
        r - (R == :closed ? 0 : eps(r))
    else
        # Empty interval
        l
    end
end

# * NumberIntervalField

"""A box where you can type in a number, within a specific interval.

## Examples
`@bind x NumberIntervalField(1.5..10.3)`

`@bind x NumberIntervalField(1.0..10.0; default=8)`

"""
struct NumberIntervalField{I<:Interval{<:Any,<:Any,<:AbstractFloat},N,Step}
    i::I
    default::N
    step::Step
end

function NumberIntervalField(i::Interval; default=missing, step="any")
    default = if default === missing
        mean(i)
    else
        closest(i, default)
    end
    NumberIntervalField(i, default, step)
end

function Base.show(io::IO, m::MIME"text/html", numberfield::NumberIntervalField)
    show(io, m, @htl("""<input $((
        type="number",
        min=closest(numberfield.i, numberfield.i.left),
        step=numberfield.step,
        max=closest(numberfield.i, numberfield.i.right),
        value=numberfield.default
    ))>"""))
end

Base.get(numberfield::NumberIntervalField) = numberfield.default
PlutoUI.BuiltinsNotebook.AbstractPlutoDingetjes.Bonds.initial_value(nf::NumberIntervalField) = nf.default
# PlutoUI.BuiltinsNotebook.AbstractPlutoDingetjes.Bonds.possible_values(nf::NumberIntervalField) = nf.range
PlutoUI.BuiltinsNotebook.AbstractPlutoDingetjes.Bonds.transform_value(nf::NumberIntervalField, val) = Base.convert(eltype(nf.i), val)
function PlutoUI.BuiltinsNotebook.AbstractPlutoDingetjes.Bonds.validate_value(nf::NumberIntervalField, val)
    val ∈ nf.i
end

# * Input

struct Input{Name,InputT,Label,Unit,Descr}
    name::Name
    input::InputT
    label::Label
    unit::Unit
    description::Descr
end
Input(name, input, label; unit=nothing, description=nothing) =
    Input(name, input, label, unit, description)

function mdjoin(v::AbstractVector{Markdown.MD})
    e = deepcopy(first(v))
    ec = e.content
    if length(ec) == 1 && first(ec) isa Markdown.Paragraph
        ec = first(ec).content
    end
    space = Markdown.parse("  ").content[1].content[1]
    for i = 2:length(v)
        c = v[i].content
        if length(c) == 1 && first(c) isa Markdown.Paragraph
            c = first(c).content
        end
        push!(ec, space)
        append!(ec, c)
    end
    e
end

function format_input(i::Input, Child, _, _)
    ustr = (isnothing(i.unit) ? "" : i.unit)
    gui_element = Child(i.name, i.input)
    label = Markdown.parse(i.label)
    entry = [mdjoin([label, Markdown.MD(gui_element), md"$(ustr)"])]
    isnothing(i.description) || push!(entry, Markdown.parse("$(i.description)"))
    entry
end

# * InputSection

struct InputSection{Name,Children}
    name::Name
    children::Children
    detail::Bool
end
InputSection(name, children; detail=false) = InputSection(name, children, detail)

function format_input(s::InputSection, Child, depth, scope)
    if scope == (:begin)
        heading = Markdown.parse(repeat("#", depth)*" "*s.name)
        s.detail ? [HTML("<details><summary>"), heading, HTML("</summary>")] : [heading]
    elseif scope == (:end)
        s.detail ? [HTML("</details>")] : []
    end
end

dfs(fun::Function, i, depth) = fun(i, depth, :middle)

function dfs(fun::Function, s::InputSection, depth=1)
    fun(s, depth, :begin)
    for c in s.children
        dfs(fun, c, depth+1)
    end
    fun(s, depth, :end)
end

# * Input tree

function format_input(input_tree)
    combine() do Child
        inputs = []

        dfs(input_tree) do e,depth,scope
            append!(inputs, format_input(e, Child, depth+1, scope))
        end

        md"""
        $(inputs)
        """
    end
end

# * Save/load inputs

function run_dir_input(default="")
    i = Input("run_dir", TextField(default=default), "Run directory")
    combine() do Child
        md"""
        $(format_input(i, Child, 0, 0))
        """
    end
end

function overwrite_input()
    i = Input("overwrite", CheckBox(default=false), "Overwrite previous run")
    combine() do Child
        md"""
        $(format_input(i, Child, 0, 0))
        """
    end
end

function list_previous_runs(;runs_dir="runs", extra_files=String[])
    dirs = filter(readdir(runs_dir)) do d
        dd = joinpath(runs_dir, d)
        isdir(dd) || return false
        for f in vcat("inputs.json", extra_files)
            isfile(joinpath(dd, f)) || return false
        end
        true
    end
    map(Base.Fix1(joinpath, runs_dir), dirs)
end

presets_input(;kwargs...) = Select(list_previous_runs(;kwargs...))

select_previous_runs(;kwargs...) = MultiSelect(list_previous_runs(;kwargs...))

function load_inputs(run_dir; verbosity=1)
    input_file = joinpath(run_dir, "inputs.json")
    !isfile(input_file) && return Dict()
    verbosity > 0 && @info "Loading inputs form $(input_file)"
    open(JSON.parse, input_file)
end

function saved_inputs(fun::Function, run_dir::AbstractString)
    stored_values = load_inputs(run_dir; verbosity=0)
    fun((k,d) -> get(stored_values, k, d))
end

function inputs_equal(a, b)
    bk = keys(b)
    Kt = eltype(bk)
    for k in keys(a)
        ka = Kt(k)
        ka ∈ bk || return false
        Vt = typeof(b[ka])
        Vt(a[k]) == b[ka] || return false
    end
    true
end

function save_inputs(inputs, run_dir::AbstractString; overwrite=false)
    run_dir = abspath(run_dir)
    if !isdir(run_dir)
        @warn "Directory $(run_dir) does not exist, please create it" run_dir
        return
    end
    input_file = joinpath(run_dir, "inputs.json")
    if isfile(input_file) && !overwrite
        inputs_equal(inputs, load_inputs(run_dir, verbosity=0)) ||
            @error "Inputs do not match those already stored in\n\n$(run_dir)\n\nplease delete run and try again, or choose another directory" inputs
        return
    end
    @info "Saving inputs to $(input_file)"
    open(input_file, "w") do file
        io = JSON.Writer.PrettyContext(file, 4)
        JSON.print(io, inputs)
    end
    input_file
end

# * Electric Fields

function get_electric_field(inputs; gauss_kind=:trunc_gauss, ramp=1.0)
    Kt = keytype(inputs)
    getk(k) = inputs[Kt(k)]

    params = Dict{Symbol,Any}()
    params[:I₀] = getk(:I₀)*u"TW/cm^2"

    τ = if Symbol(getk(:λorω)) == :λ
        λ = getk(:λ)*u"nm"
        params[:λ] = λ
        getk(:cycles)*u"fs"(λ/u"c")
    else
        params[:ω] = getk(:ω)
        getk(:cycles)*(2π/getk(:ω))
    end

    if Symbol(getk(:envelope)) == :gaussian
        params[:τ] = τ
        params[:σoff] = 3.0
        params[:σmax] = 4.0
        params[:env] = gauss_kind
    elseif Symbol(getk(:envelope)) == :tophat
        params[:flat] = float(getk(:cycles))
        params[:ramp] = ramp
        params[:env] = :tophat
        params[:ramp_kind] = :sin²
    else
        @error "Unknown envelope $(getk(:envelope))"
    end

    ElectricFields.make_field(params)
end

end # module NotebookHelper
