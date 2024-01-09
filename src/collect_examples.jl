using TOML

function collect_examples(folder::AbstractString)
    examples = []
    for name in readdir(folder; join=false, sort=false)
        isdir(joinpath(folder, name)) || continue
        example = TOML.parsefile(joinpath(folder, name, "metadata.toml"))
        push!(examples, name => example)
    end
    sort!(examples; by=(pair -> get(pair[2], "order", 0)))
    # TODO: filter option
    return examples
end
