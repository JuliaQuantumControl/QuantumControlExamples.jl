"""
    clean([distclean=false])

Clean up build/doc/testing artifacts. Restore to clean checkout state
(distclean)
"""
function clean(; distclean=false, _exit=true)

    _exists(name) = isfile(name) || isdir(name)
    _push!(lst, name) = _exists(name) && push!(lst, name)

    function _glob(folder, ending)
        if !_exists(folder)
            return []
        end
        [name for name in readdir(folder; join=true) if (name |> endswith(ending))]
    end

    function _glob_star(folder; except=[])
        if !_exists(folder)
            return []
        end
        [
            joinpath(folder, name) for
            name in readdir(folder) if !(name |> startswith(".") || name ∈ except)
        ]
    end


    ROOT = @__DIR__

    ###########################################################################
    CLEAN = String[]
    _push!(CLEAN, joinpath(ROOT, "coverage"))
    _push!(CLEAN, joinpath(ROOT, "docs", "build"))
    append!(CLEAN, _glob(ROOT, ".info"))
    append!(CLEAN, _glob(joinpath(ROOT, ".coverage"), ".info"))
    ###########################################################################

    ###########################################################################
    DISTCLEAN = String[]
    _push!(DISTCLEAN, joinpath(ROOT, "docs", "src", "examples"))
    _push!(DISTCLEAN, joinpath(ROOT, "docs", "src", "tutorials"))
    _push!(DISTCLEAN, joinpath(ROOT, "Manifest.toml"))
    ###########################################################################

    for name in CLEAN
        @info "rm $name"
        rm(name, force=true, recursive=true)
    end
    clean_examples(joinpath(ROOT, "examples"))
    clean_examples(joinpath(ROOT, "tutorials"))
    if distclean
        for name in DISTCLEAN
            @info "rm $name"
            rm(name, force=true, recursive=true)
        end
        if _exit
            @info "Exiting"
            exit(0)
        end
    end

end

include("src/collect_examples.jl")
include("src/clean_examples.jl")
