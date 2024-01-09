# JuliaQuantumControl Tutorials and Examples

Examples and Tutorials for the [JuliaQuantumControl organization](https://github.com/JuliaQuantumControl/)

The tutorials and examples are organized in the `tutorials` and `examples` subfolders. Each subfolder contains a script file that is runnable as a test and can be converted to `.md` and `.ipynb` files via [Literate.jl](https://github.com/fredrikekre/Literate.jl).

This conversion should happen on a dedicated multi-core workstation, not CI. This is to show the examples under realistic conditions. For deployment, the converted files will be copied into the `docs/src` folder, so that they can be build into HTML without rerunning them. Building and deploying the HTML happens in CI.

## Development

Run `make` after checking out this repository.

## Deployment

In order to deploy the tutorial and example to https://juliaquantumcontrol.github.io/QuantumControlExamples.jl/, do the following on a workstation:

```
git checkout -b 2024-01
JULIA_NUM_THREADS=auto make all
git add .
git commit -m "Build 2024-01"
git push origin 2024-01
```

This locally runs the examples and writes them to the `docs/src` folder. Committing and pushing the resulting files to a branch named according to `YYYY-MM` will then build and deploy the website automatically.
