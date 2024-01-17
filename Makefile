.PHONY: help all deploy test docs clean distclean devrepl codestyle servedocs
.DEFAULT_GOAL := help

JULIA ?= julia
PORT ?= 8000

define PRINT_HELP_JLSCRIPT
rx = r"^([a-z0-9A-Z_-]+):.*?##[ ]+(.*)$$"
for line in eachline()
    m = match(rx, line)
    if !isnothing(m)
        target, help = m.captures
        println("$$(rpad(target, 20)) $$help")
    end
end
endef
export PRINT_HELP_JLSCRIPT

Manifest.toml: Project.toml
	$(JULIA) --project=. --banner=no --startup-file=no -e 'using Pkg; Pkg.activate("."); Pkg.instantiate(); Pkg.resolve()'

help:  ## show this help
	@julia -e "$$PRINT_HELP_JLSCRIPT" < $(MAKEFILE_LIST)


all: Manifest.toml ## execute all examples and deploy to the docs folder
	$(JULIA) --project=. --banner=no --startup-file=no -e 'using QuantumControlExamples; QuantumControlExamples.run_examples("examples")'
	$(JULIA) --project=. --banner=no --startup-file=no -e 'using QuantumControlExamples; QuantumControlExamples.deploy_examples("examples")'
	$(JULIA) --project=. --banner=no --startup-file=no -e 'using QuantumControlExamples; QuantumControlExamples.run_examples("tutorials")'
	$(JULIA) --project=. --banner=no --startup-file=no -e 'using QuantumControlExamples; QuantumControlExamples.deploy_examples("tutorials")'

deploy: Manifest.toml ## deploy all run examples to the docs folder
	$(JULIA) --project=. --banner=no --startup-file=no -e 'using QuantumControlExamples; QuantumControlExamples.deploy_examples("examples")'
	$(JULIA) --project=. --banner=no --startup-file=no -e 'using QuantumControlExamples; QuantumControlExamples.deploy_examples("tutorials")'

test:  Manifest.toml ## Run the test suite
	$(JULIA) --project=. --banner=no --startup-file=yes -e 'include("devrepl.jl"); test()'
	@echo "Done. Consider using 'make devrepl'"

devrepl:  Manifest.toml ## Start an interactive REPL for testing and building documentation
	JULIA_NUM_THREADS=auto $(JULIA) --project=. --banner=no --startup-file=yes -i devrepl.jl

docs: Manifest.toml ## Build the documentation
	$(JULIA) --project=. docs/make.jl
	@echo "Done. Consider using 'make devrepl'"

servedocs: Manifest.toml ## Build (auto-rebuild) and serve documentation at PORT=8000
	$(JULIA) --project=. -e 'include("devrepl.jl"); servedocs(port=$(PORT), verbose=true)'

clean:  ## Clean up build/doc/testing artifacts
	$(JULIA) -e 'include("clean.jl"); clean()'

codestyle: Manifest.toml ## Apply the codestyle to the entire project
	$(JULIA) --project=. -e 'using JuliaFormatter; format(".", verbose=true)'
	@echo "Done. Consider using 'make devrepl'"

distclean:  ## Restore to a clean checkout state
	$(JULIA) -e 'include("clean.jl"); clean(distclean=true)'
