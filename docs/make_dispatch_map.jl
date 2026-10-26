# Writes `docs/dispatch-map.md` from the loaded package: the supertypes, the hooks each one's
# docstring lists under "# Methods a subtype provides", which built-in types implement them, and
# what each pipeline step reads and calls. The "Dispatch map" test fails when the file is stale.
#
#     julia --project docs/make_dispatch_map.jl

using FjordSim
using FjordSim: Configs, CLI, Simulations
using FjordSim.Configs: AbstractFjordConfig
using FjordSim.Validation: AbstractObservationConfig
using InteractiveUtils: subtypes

const DISPATCH_MAP_PATH = joinpath(@__DIR__, "dispatch-map.md")

dispatch_supertypes() = sort!(
    [
        [getfield(Configs, name) for name in names(Configs; all=true)
         if startswith(String(name), "Abstract") && getfield(Configs, name) isa Type]
        AbstractObservationConfig
    ];
    by = nameof,
)

function package_modules(root = FjordSim, found = Module[])
    push!(found, root)
    for name in names(root; all=true)
        isdefined(root, name) || continue
        value = getfield(root, name)
        if value isa Module && parentmodule(value) === root && value !== root && !(value in found)
            package_modules(value, found)
        end
    end
    return found
end

in_package(mod::Module) = Base.moduleroot(mod) === FjordSim

function implementations(supertype)
    found = Any[]
    for type in subtypes(supertype)
        if isabstracttype(Base.unwrap_unionall(type))
            append!(found, implementations(type))
        elseif in_package(parentmodule(type))
            push!(found, type)
        end
    end
    return sort!(found; by = nameof)
end

function docstring(type)
    binding = Base.Docs.Binding(parentmodule(type), nameof(type))
    return join(only(values(Base.Docs.meta(parentmodule(type))[binding].docs)).text)
end

# The bullets under "# Methods a subtype provides", continuation lines joined.
function hook_bullets(type)
    lines = split(docstring(type), '\n')
    start = findfirst(==("# Methods a subtype provides"), lines)
    isnothing(start) && error("`$(nameof(type))` has no \"# Methods a subtype provides\" docstring section")
    bullets = String[]
    for line in lines[start+1:end]
        if startswith(line, "- ")
            push!(bullets, line[3:end])
        elseif startswith(line, "  ") && !isempty(bullets)
            bullets[end] *= " " * strip(line)
        else
            break
        end
    end
    return bullets
end

# The names a bullet declares: the code spans before its first colon outside code spans.
function declared_names(bullet)
    in_code = false
    head = bullet
    for (index, character) in pairs(bullet)
        character == '`' && (in_code = !in_code)
        if character == ':' && !in_code
            head = bullet[1:prevind(bullet, index)]
            break
        end
    end
    spans = [match.captures[1] for match in eachmatch(r"`([^`]*)`", head)]
    return [first(split(span, '(')) for span in spans if span != "nothing"]
end

function hook_status(bullet)
    occursin("Required", bullet) && return "required"
    conditional = match(r"Only if ([^.]*)\.", bullet)
    isnothing(conditional) || return "only if " * conditional.captures[1]
    occursin(r"Optional|[Dd]efaults to", bullet) && return "optional"
    return ""
end

function resolve_function(name, modules)
    if startswith(name, "FjordSim.")
        return foldl(getfield, Symbol.(split(name, '.')[2:end]); init = FjordSim)
    end
    symbol = Symbol(name)
    candidates = unique([getfield(mod, symbol) for mod in modules
                         if isdefined(mod, symbol) && getfield(mod, symbol) isa Function])
    length(candidates) == 1 && return only(candidates)
    error("hook `$name` names $(length(candidates)) functions in FjordSim; expected one")
end

function argument_types(method)
    parameters = Base.unwrap_unionall(method.sig).parameters[2:end]
    return [Base.rewrap_unionall(parameter isa TypeVar ? parameter.ub : parameter, method.sig)
            for parameter in parameters]
end

# "own" for a method on `type` or a supertype below `supertype`, "default" for one on
# `supertype` itself or wider, "" for none.
function implementation(hook, type, supertype)
    found = ""
    for method in methods(hook)
        in_package(method.module) || continue
        for argument in argument_types(method)
            (argument === Any || !(type <: argument)) && continue
            supertype <: argument || return "own"
            found = "default"
        end
    end
    return found
end

function cell(kind, status)
    kind == "own" && return "✓"
    kind == "default" && return "default"
    return status == "required" ? "**missing**" : "—"
end

function hook_table(io, supertype, modules)
    types = implementations(supertype)
    print(io, "| Hook | Status |")
    foreach(type -> print(io, " `", nameof(type), "` |"), types)
    print(io, "\n|---|---|", "---|"^length(types), "\n")
    for bullet in hook_bullets(supertype), name in declared_names(bullet)
        hook = resolve_function(name, modules)
        print(io, "| `", name, "` | ", hook_status(bullet), " |")
        foreach(type -> print(io, " ", cell(implementation(hook, type, supertype), hook_status(bullet)), " |"), types)
        print(io, "\n")
    end
    isempty(types) && print(io, "\nNo built-in implementation.\n")
    return nothing
end

callee(code, target) = target isa Core.SSAValue ? callee(code, code[target.id]) :
    target isa GlobalRef ? getglobal(target.mod, target.name) : target

# What a step's one method reads from the setup, the function its result comes from, and whether
# that function is handed the whole setup. A step with keywords is followed into its body function.
function step_calls(step, accessors)
    info = only(Base.code_lowered(step, Tuple{AbstractFjordConfig}))
    setup = Base.method_argnames(only(methods(step)))[2]
    return lowered_calls(info, setup, accessors)
end

function lowered_calls(info, setup, accessors)
    code = info.code
    functions = [callee(code, statement.args[1]) for statement in code if Meta.isexpr(statement, :call)]
    reads = unique([nameof(f) for f in functions if f in accessors])
    # A splatted keyword set lowers to two returns, without and with keywords; the last has them.
    returned = code[last([statement for statement in code if statement isa Core.ReturnNode]).val.id]
    work = callee(code, returned.args[1])
    if work === Core.kwcall
        work = callee(code, returned.args[3])
    end
    if startswith(String(nameof(work)), '#')
        return lowered_calls(only(Base.code_lowered(work)), setup, accessors)
    end
    whole = any(argument -> argument isa Core.SlotNumber && info.slotnames[argument.id] == setup,
                returned.args[2:end])
    return reads, work, whole
end

function component_supertype(name, supertypes)
    stem = replace(String(name), r"_configs?$" => "")
    for candidate in (stem, chopsuffix(stem, "s"))
        camel = join(uppercasefirst.(split(candidate, '_')))
        index = findfirst(type -> nameof(type) == Symbol("Abstract", camel, "Config"), supertypes)
        isnothing(index) || return supertypes[index]
    end
    return nothing
end

function component_graph(io, supertypes, modules)
    edges = Tuple{Any,Symbol,Any}[]
    for supertype in supertypes, bullet in hook_bullets(supertype), name in declared_names(bullet)
        component = component_supertype(last(split(name, '.')), supertypes)
        isnothing(component) || component === supertype || push!(edges, (supertype, Symbol(name), component))
    end
    print(io, "```mermaid\nflowchart LR\n")
    for (from, accessor, to) in edges
        print(io, "    ", nameof(from), " -- ", accessor, " --> ", nameof(to), "\n")
    end
    print(io, "```\n")
    reached = Set(Any[AbstractFjordConfig; [to for (_, _, to) in edges]])
    unreached = [nameof(type) for type in supertypes if !(type in reached)]
    isempty(unreached) || print(io, "\nNot reached through an accessor: ",
                                join(["`$name`" for name in unreached], ", "), ".\n")
    return nothing
end

function dispatch_map()
    modules = package_modules()
    supertypes = dispatch_supertypes()
    accessors = [resolve_function(name, modules) for bullet in hook_bullets(AbstractFjordConfig)
                 for name in declared_names(bullet) if endswith(name, r"_configs?")]
    io = IOBuffer()
    print(io, """
        # Dispatch map

        Generated by `julia --project docs/make_dispatch_map.jl` from the loaded package; do not edit by
        hand. The "Dispatch map" test fails when this file no longer matches the code.

        A hook is a function a supertype's docstring lists under "# Methods a subtype provides". A cell
        says whether a built-in type has a method of its own (✓), inherits the supertype's (default), or
        has none (—, or **missing** for a required hook).

        ## Components

        Which supertype each accessor returns, by the accessor's name.

        """)
    component_graph(io, supertypes, modules)
    print(io, """

        ## Steps

        Each step has one method, on `AbstractFjordConfig`. It reads components through the accessors
        and hands them to the work function, which dispatches on them.

        | Step | Reads | Work function |
        |---|---|---|
        """)
    steps = [last.(CLI.SUBCOMMANDS); Simulations.build_simulation]
    for step in steps
        reads, work, whole = step_calls(step, accessors)
        read_list = join(["`$name`" for name in reads], ", ") * (whole ? " and the whole setup" : "")
        print(io, "| `", nameof(step), "` | ", read_list, " | `", nameof(work), "` |\n")
    end
    print(io, "\n## Hooks\n")
    for supertype in supertypes
        print(io, "\n### `", nameof(supertype), "`\n\n")
        hook_table(io, supertype, modules)
    end
    return String(take!(io))
end

if abspath(PROGRAM_FILE) == @__FILE__
    write(DISPATCH_MAP_PATH, dispatch_map())
end
