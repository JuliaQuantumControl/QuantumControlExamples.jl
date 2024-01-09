using TOML
using Glob


# Assume that `run_example(folder)` has been executed and deploy the resulting
# files the the appropriate `docs/src` subfolder.
function deploy_example(folder::AbstractString,)
    if contains(abspath(folder), r"[/\\]tutorials[/\\]")
        target = normpath(@__DIR__, "..", "docs", "src", "tutorials")
    else
        target = normpath(@__DIR__, "..", "docs", "src", "examples")
    end
    example = TOML.parsefile(joinpath(folder, "metadata.toml"))
    rootname = splitext(example["file"])[1]
    example_md = rootname * ".md"
    example_jl = rootname * ".jl"
    example_ipynb = rootname * ".ipynb"
    files_to_copy = [example_md, example_jl, example_ipynb, "metadata.toml"]
    for pattern in get(example, "deploy", [])
        for path in glob(pattern, folder)
            push!(files_to_copy, relpath(path, folder))
        end
    end
    missing_files = filter(name -> !isfile(joinpath(folder, name)), files_to_copy)
    if isempty(missing_files)
        mkpath(joinpath(target, rootname))
        for name in files_to_copy
            target_name = name
            if name == example_md
                target_name = "index.md"
            end
            println(
                "cp $(joinpath(folder, name)) $(joinpath(target, rootname, target_name))"
            )
            cp(joinpath(folder, name), joinpath(target, rootname, target_name), force=true)
        end
    else
        @warn "Example in $folder has not been run." missing_files
    end
    return nothing
end


# Recursively deploy all examples in the `folder` to the docs directory
function deploy_examples(folder)
    # TODO: allow to filter, manage parallel execution in subprocesses
    for (subfolder, example) in collect_examples(folder)
        deploy_example(joinpath(folder, subfolder))
    end
end
