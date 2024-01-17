using TOML

function _script_run_test(folder; quiet=false)
    folder = abspath(folder)
    example = TOML.parsefile(joinpath(folder, "metadata.toml"))
    file = example["file"]
    if quiet
        logfile = splitext(file)[1] * ".log"
        script = """
        ENV["GKSwstype"] = "100"
        ENV["JULIA_CAPTURE_COLOR"] = "0"
        using Plots
        using TOML
        using Test
        cd($(repr(folder))) do
            @info "writing log to `$(joinpath(folder, logfile))`."
            open($(repr(logfile)), "w") do io
                redirect_stdout(io) do
                    redirect_stderr(io) do
                        Plots.with(:unicodeplots) do
                            example = TOML.parsefile("metadata.toml")
                            @testset "$(example["title"])" begin
                                include($(repr(joinpath(folder, file))))
                            end
                        end
                    end
                end
            end
        end
        """
    else
        script = """
        ENV["GKSwstype"] = "100"
        ENV["JULIA_CAPTURE_COLOR"] = "0"
        using Plots
        using TOML
        using Test
        cd($(repr(folder))) do
            Plots.with(:unicodeplots) do
                example = TOML.parsefile("metadata.toml")
                @testset "$(example["title"])" begin
                    include($(repr(joinpath(folder, file))))
                end
            end
        end
        """
    end
    return script
end


function _script_run_markdown(folder)
    folder = abspath(folder)
    example = TOML.parsefile(joinpath(folder, "metadata.toml"))
    script = """
    ENV["GKSwstype"] = "100"
    ENV["JULIA_CAPTURE_COLOR"] = "0"
    using Literate
    cd($(repr(folder))) do
        Literate.markdown(
            $(repr(example["file"])), "."; execute=true, credit=false
        )
    end
    """
    return script
end


function _script_run_notebook(folder)
    folder = abspath(folder)
    example = TOML.parsefile(joinpath(folder, "metadata.toml"))
    script = """
    ENV["GKSwstype"] = "100"
    ENV["JULIA_CAPTURE_COLOR"] = "0"
    using Literate
    cd($(repr(folder))) do
        Literate.notebook(
            $(repr(example["file"])), "."; execute=true, credit=false
        )
    end
    """
    return script
end


function _run(code; threads=1)
    julia = Base.julia_cmd().exec[1]
    project = dirname(Base.active_project())
    run_jl = tempname()
    write(run_jl, code)
    cmd = [julia, "--project=$project", "--startup-file=no", "--threads=$threads", run_jl]
    run(Cmd(cmd); wait=true)
end


function _shrinkuser(path::String)
    home = homedir()
    if startswith(path, home)
        return replace(path, home => "~")
    end
    return path
end



function run_example(
    folder::AbstractString;
    run_test=true,
    run_markdown=true,
    run_notebook=true,
    quiet=false,
)
    example = TOML.parsefile(joinpath(folder, "metadata.toml"))
    threads = get(example, "threads", 1)
    file = abspath(joinpath(folder, example["file"]))
    if run_test
        up_to_date = false
        if quiet
            logfile = splitext(file)[1] * ".log"
            if mtime(logfile) > mtime(file)
                up_to_date = true
            end
        end
        if up_to_date
            @info "skip running `$(_shrinkuser(file))` as script (up to date)."
        else
            @info "running `$(_shrinkuser(file))` as script."
            _run(_script_run_test(folder; quiet); threads)
        end
    end
    if run_markdown
        mdfile = splitext(file)[1] * ".md"
        if mtime(mdfile) > mtime(file)
            @info "skip converting `$(_shrinkuser(file))` to .md file (up to date)."
        else
            _run(_script_run_markdown(folder); threads)
        end
    end
    if run_notebook
        nbfile = splitext(file)[1] * ".ipynb"
        if mtime(nbfile) > mtime(file)
            @info "skip converting `$(_shrinkuser(file))` to .ipynb file (up to date)."
        else
            _run(_script_run_notebook(folder); threads)
        end
    end
end


mutable struct TaskGroup
    tasks::Vector{Any}
    available_threads::Int64
    function TaskGroup(threads::Int64)
        new([], threads)
    end
end



function run_examples(folder; nthreads=Base.Threads.nthreads())
    if nthreads > Base.Threads.nthreads()
        @error "Cannot run with nthreads=$nthreads > $(Base.Threads.nthreads())"
    end
    @info "Running examples in `$folder` with $nthreads threads"
    if nthreads == 1
        for (subfolder, example) in collect_examples(folder)
            path = joinpath(folder, subfolder)
            run_example(path; quiet=true)
        end
    else
        groups = []
        for (subfolder, example) in collect_examples(folder)
            example_threads = get(example, "threads", 1)
            path = joinpath(folder, subfolder)
            scheduled = false
            for group in groups
                if group.available_threads >= example_threads
                    push!(group.tasks, path)
                    group.available_threads -= example_threads
                    scheduled = true
                    break
                end
            end
            if !scheduled
                group = TaskGroup(nthreads)
                push!(group.tasks, path)
                group.available_threads -= example_threads
                push!(groups, group)
            end
        end
        for group in groups
            @info "Running group of examples: $(group.tasks)"
            Base.Threads.@threads :static for path in group.tasks
                run_example(path; quiet=true)
            end
            @info "Done with group of examples: $(group.tasks)"
        end
    end
end
