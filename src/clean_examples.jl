using TOML

# Clean up generate files in the given `folder`, which should be a subfolder of
# the root `examples` or `tutorials` directory
function clean_example(folder::AbstractString)
    cd(folder) do
        example = TOML.parsefile("metadata.toml")
        for name in readdir(folder; join=false, sort=false)
            (name == ".gitignore") && continue
            (name == "metadata.toml") && continue
            (name == example["file"]) && continue
            (name in get(example, "assets", [])) && continue
            # Note that `assets` are not glob patterns. This is because this
            # file must run in a plain Julia environment and can't use the Glob
            # dependency.
            path = joinpath(folder, name)
            println("rm $path")
            rm(path; recursive=true)
        end
    end
end


# Recursively clean up all examples in the `folder`, which should be the root
# `examples` or `tutorials` directory.
function clean_examples(folder::AbstractString)
    for (subfolder, _) in collect_examples(folder)
        clean_example(joinpath(folder, subfolder))
    end
end
