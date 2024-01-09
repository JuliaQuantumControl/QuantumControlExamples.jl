# Source this script as e.g.
#
#     include("PATH/TO/devrepl.jl")
#
# from *any* Julia REPL or run it as e.g.
#
#     julia -i --banner=no PATH/TO/devrepl.jl
#
# from anywhere. This will change the current working directory and
# activate/initialize the correct Julia environment for you.
#
# You may also run this in vscode to initialize a development REPL
#
using Pkg

cd(@__DIR__)
Pkg.activate(".")
Pkg.instantiate()


using Revise
using JuliaFormatter
using LiveServer: LiveServer, serve, servedocs

using QuantumControlExamples:
    QuantumControlExamples,
    run_example,
    run_examples,
    clean_example,
    clean_examples,
    deploy_example,
    deploy_examples
using QuantumControlTestUtils: test as _test, show_coverage, generate_coverage_html


REPL_MESSAGE = """
*******************************************************************************
DEVELOPMENT REPL

Revise, JuliaFormatter, LiveServer are active.

* `help()` – Show this message
* `run_example(folder)` – Process an example folder with Literate
* `run_examples("examples") – Run `run_example` for all examples
* `run_examples("tutorials") – Run `run_example` for all tutorials
* `clean_example(folder)` – Remove generated files in example folder
* `include("docs/make.jl")` – Generate the documentation
* `test()` – Run the entire test suite in a subprocess with coverage
* `show_coverage()` – Print a tabular overview of coverage data
* `generate_coverage_html()` – Generate an HTML coverage report
* `format(".")` – Apply code formatting to all files
* `clean()` – Clean up build/doc/testing artifacts
* `distclean()` – Restore to a clean checkout state
*******************************************************************************
"""

"""Show help"""
help() = println(REPL_MESSAGE)


module Clean
include("clean.jl")
end

using .Clean: clean as _clean

function clean(; distclean=false)
    _clean(; distclean, _exit=false)
end

distclean() = clean(distclean=true)


function test(args...; kwargs...)
    return _test(args...; project=".", kwargs...)
end

if abspath(PROGRAM_FILE) == @__FILE__
    # REPL mode (as opposed to `include("devrepl.jl")`)
    help()
end
