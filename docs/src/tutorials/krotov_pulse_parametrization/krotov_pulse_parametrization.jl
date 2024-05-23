# # Pulse Parametrization for Krotov's Method

#md # ``\gdef\op#1{\hat{#1}}``
#md # ``\gdef\init{\text{init}}``
#md # ``\gdef\tgt{\text{tgt}}``

#nb # $
#nb # \newcommand{tr}[0]{\operatorname{tr}}
#nb # \newcommand{diag}[0]{\operatorname{diag}}
#nb # \newcommand{abs}[0]{\operatorname{abs}}
#nb # \newcommand{pop}[0]{\operatorname{pop}}
#nb # \newcommand{aux}[0]{\text{aux}}
#nb # \newcommand{opt}[0]{\text{opt}}
#nb # \newcommand{tgt}[0]{\text{tgt}}
#nb # \newcommand{init}[0]{\text{init}}
#nb # \newcommand{lab}[0]{\text{lab}}
#nb # \newcommand{rwa}[0]{\text{rwa}}
#nb # \newcommand{bra}[1]{\langle#1\vert}
#nb # \newcommand{ket}[1]{\vert#1\rangle}
#nb # \newcommand{Bra}[1]{\left\langle#1\right\vert}
#nb # \newcommand{Ket}[1]{\left\vert#1\right\rangle}
#nb # \newcommand{Braket}[2]{\left\langle #1\vphantom{#2}\mid{#2}\vphantom{#1}\right\rangle}
#nb # \newcommand{op}[1]{\hat{#1}}
#nb # \newcommand{Op}[1]{\hat{#1}}
#nb # \newcommand{dd}[0]{\,\text{d}}
#nb # \newcommand{Liouville}[0]{\mathcal{L}}
#nb # \newcommand{DynMap}[0]{\mathcal{E}}
#nb # \newcommand{identity}[0]{\mathbf{1}}
#nb # \newcommand{Norm}[1]{\lVert#1\rVert}
#nb # \newcommand{Abs}[1]{\left\vert#1\right\vert}
#nb # \newcommand{avg}[1]{\langle#1\rangle}
#nb # \newcommand{Avg}[1]{\left\langle#1\right\rangle}
#nb # \newcommand{AbsSq}[1]{\left\vert#1\right\vert^2}
#nb # \newcommand{Re}[0]{\operatorname{Re}}
#nb # \newcommand{Im}[0]{\operatorname{Im}}
#nb # $

#md import DisplayAs #hide

# This example illustrates the parametrization of control pulses as a
# form of amplitude constraint.

datadir(names...) = joinpath(@__DIR__, names...);
#-
using QuantumControl
using QuantumControl.Shapes: flattop
using QuantumControl.Generators
using QuantumControl.Controls
using QuantumControl.PulseParametrizations:
    SquareParametrization,
    TanhParametrization,
    TanhSqParametrization,
    LogisticParametrization,
    LogisticSqParametrization,
    ParametrizedAmplitude
using QuantumPropagators: ExpProp
using LinearAlgebra

#-
using Test #src
using Plots
Plots.default(
    linewidth               = 3,
    size                    = (550, 300),
    legend                  = :top,
    foreground_color_legend = nothing,
    background_color_legend = RGBA(1, 1, 1, 0.8),
)
#-

# ## Parametrizations

# ## Symmetric Bounded Controls

#-
include(joinpath(@__DIR__, "plots", "symmetric_parametrization_comparison.jl"))  #hide
fig = plot_symmetric_parametrization_comparison()  #hide
#md fig |> DisplayAs.PNG #hide
display(fig) #src

# ## Positive (Bounded) Controls

#-
include(joinpath(@__DIR__, "plots", "positive_parametrization_comparison.jl"))  #hide
fig = plot_positive_parametrization_comparison()  #hide
#md fig |> DisplayAs.PNG #hide
display(fig) #src

# ## Two-level Hamiltonian

# We consider the Hamiltonian $\op{H}_{0} = - \frac{\omega}{2} \op{\sigma}_{z}$, representing
# a simple qubit with energy level splitting $\omega$ in the basis
# $\{\ket{0},\ket{1}\}$. The control field $\epsilon(t)$ is assumed to couple via
# the Hamiltonian $\op{H}_{1}(t) = \epsilon(t) \op{\sigma}_{x}$ to the qubit,
# i.e., the control field effectively drives transitions between both qubit
# states.
#
# We we will use

ϵ(t) = 0.2 * flattop(t, T=5, t_rise=0.3, func=:blackman);


#-
"""Two-level-system Hamiltonian."""
function tls_hamiltonian(; Ω=1.0, ampl=ϵ)
    σ̂_z = ComplexF64[
        1  0
        0 -1
    ]
    σ̂_x = ComplexF64[
        0  1
        1  0
    ]
    Ĥ₀ = -0.5 * Ω * σ̂_z
    Ĥ₁ = σ̂_x
    return hamiltonian(Ĥ₀, (Ĥ₁, ampl))
end;
#-

H = tls_hamiltonian();
@test length(H.ops) == 2 #src

# The control field here switches on from zero at $t=0$ to it's maximum amplitude
# 0.2 within the time period 0.3 (the switch-on shape is half a [Blackman pulse](https://en.wikipedia.org/wiki/Window_function#Blackman_window)).
# It switches off again in the time period 0.3 before the
# final time $T=5$). We use a time grid with 500 time steps between 0 and $T$:

tlist = collect(range(0, 5, length=500));
nothing #hide


#-
function plot_amplitude(ampl, tlist)
    plot(tlist, discretize(ampl, tlist), xlabel="time", ylabel="amplitude", legend=false)
end

fig = plot_amplitude(ϵ, tlist)
#md fig |> DisplayAs.PNG #hide
display(fig) #src

# ## Optimization target

# The `Krotov` package requires the optimization to be expressed over a set of
# "trajectories". In this example, there is only a single
# trajectory: the state-to-state transfer from initial state $\ket{\Psi_{\init}} =
# \ket{0}$ to the target state $\ket{\Psi_{\tgt}} = \ket{1}$, under the dynamics
# of the Hamiltonian $\op{H}(t)$:

function ket(label)
    result = Dict("0" => Vector{ComplexF64}([1, 0]), "1" => Vector{ComplexF64}([0, 1]),)
    return result[string(label)]
end;

#-
@test dot(ket(0), ket(1)) ≈ 0 #src
#-

trajectories = [Trajectory(ket(0), H, target_state=ket(1))]

#-
@test length(trajectories) == 1 #src
#-


# ## Square-parametrization for positive pulses

a = ParametrizedAmplitude(
    ϵ,
    tlist;
    parametrization=SquareParametrization(),
    parameterize=true
)

#-
function plot_amplitude(ampl::ParametrizedAmplitude, tlist)
    plot(
        tlist,
        discretize(Array(ampl), tlist),
        xlabel="time",
        ylabel="amplitude",
        legend=false
    )
end

fig = plot_amplitude(a, tlist)
#md fig |> DisplayAs.PNG #hide
display(fig) #src

#-
problem = ControlProblem(
    trajectories=substitute(trajectories, IdDict(ϵ => a)),
    prop_method=ExpProp,
    lambda_a=5,
    update_shape=(t -> flattop(t, T=5, t_rise=0.3, func=:blackman)),
    tlist=tlist,
    iter_stop=50,
    J_T=QuantumControl.Functionals.J_T_ss,
    check_convergence=res -> begin
        ((res.J_T < 1e-3) && (res.converged = true) && (res.message = "J_T < 10⁻³"))
    end
);
#-
using Krotov
opt_result_positive = @optimize_or_load(
    datadir("parametrization#opt_result_positive.jld2"),
    problem;
    method=Krotov
);
#-
opt_result_positive

# We can plot the optimized field:

#-
fig = plot_amplitude(
    substitute(a, IdDict(a.control => opt_result_positive.optimized_controls[1])),
    tlist
)
display(fig) #src
#md fig |> DisplayAs.PNG #hide
#-

#-
amplitude = #src
    Array(substitute(a, IdDict(a.control => opt_result_positive.optimized_controls[1]))) #src
@test minimum(amplitude) ≥ 0.0 #src
@test minimum(amplitude) < 1e-16 #src
@test maximum(amplitude) > 0.0 #src
#-

# ## Tanh-Square-Parametrization for positive amplitude-constrained pulses

a = ParametrizedAmplitude(
    ϵ,
    tlist;
    parametrization=TanhSqParametrization(3),
    parameterize=true
)

problem_tanhsq = ControlProblem(
    trajectories=substitute(trajectories, IdDict(ϵ => a)),
    prop_method=ExpProp,
    lambda_a=10,
    update_shape=(t -> flattop(t, T=5, t_rise=0.3, func=:blackman)),
    tlist=tlist,
    iter_stop=50,
    J_T=QuantumControl.Functionals.J_T_ss,
    check_convergence=res -> begin
        ((res.J_T < 1e-3) && (res.converged = true) && (res.message = "J_T < 10⁻³"))
    end
);
#-
opt_result_tanhsq = @optimize_or_load(
    datadir("parametrization#opt_result_tanhsq.jld2"),
    problem_tanhsq;
    method=Krotov
);
#-
opt_result_tanhsq

# We can plot the optimized field:

#-
fig = plot_amplitude(
    substitute(a, IdDict(a.control => opt_result_tanhsq.optimized_controls[1])),
    tlist
)
display(fig) #src
#md fig |> DisplayAs.PNG #hide
#-

#-
amplitude = #src
    Array(substitute(a, IdDict(a.control => opt_result_tanhsq.optimized_controls[1]))) #src
@test minimum(amplitude) ≥ 0.0 #src
@test minimum(amplitude) < 1e-16 #src
@test maximum(amplitude) > 0.0 #src
@test maximum(opt_result_tanhsq.optimized_controls[1]) < 3.0 #src
#-

# ## Logistic-Square-Parametrization for positive amplitude-constrained pulses

a = ParametrizedAmplitude(
    ϵ,
    tlist;
    parametrization=LogisticSqParametrization(3, k=1.0),
    parameterize=true
)

problem_logisticsq = ControlProblem(
    trajectories=substitute(trajectories, IdDict(ϵ => a)),
    prop_method=ExpProp,
    lambda_a=1,
    update_shape=(t -> flattop(t, T=5, t_rise=0.3, func=:blackman)),
    tlist=tlist,
    iter_stop=50,
    J_T=QuantumControl.Functionals.J_T_ss,
    check_convergence=res -> begin
        ((res.J_T < 1e-3) && (res.converged = true) && (res.message = "J_T < 10⁻³"))
    end
);
#-
opt_result_logisticsq = @optimize_or_load(
    datadir("parametrization#opt_result_logisticsq.jld2"),
    problem_logisticsq;
    method=Krotov
);
# We can plot the optimized field:

#-
fig = plot_amplitude(
    substitute(a, IdDict(a.control => opt_result_logisticsq.optimized_controls[1])),
    tlist
)
display(fig) #src
#md fig |> DisplayAs.PNG #hide
#-

#-
amplitude = #src
    Array(substitute(a, IdDict(a.control => opt_result_logisticsq.optimized_controls[1]))) #src
@test minimum(amplitude) ≥ 0.0 #src
@test minimum(amplitude) < 1e-16 #src
@test maximum(amplitude) > 0.0 #src
@test maximum(amplitude) < 3.0 #src
#-
# ## Tanh-parametrization for amplitude-constrained pulses

a = ParametrizedAmplitude(
    ϵ,
    tlist;
    parametrization=TanhParametrization(-0.5, 0.5),
    parameterize=true
)

problem_tanh = ControlProblem(
    trajectories=substitute(trajectories, IdDict(ϵ => a)),
    prop_method=ExpProp,
    lambda_a=1,
    update_shape=(t -> flattop(t, T=5, t_rise=0.3, func=:blackman)),
    tlist=tlist,
    iter_stop=50,
    J_T=QuantumControl.Functionals.J_T_ss,
    check_convergence=res -> begin
        ((res.J_T < 1e-3) && (res.converged = true) && (res.message = "J_T < 10⁻³"))
    end
);
#-
opt_result_tanh = @optimize_or_load(
    datadir("parametrization#opt_result_tanh.jld2"),
    problem_tanh;
    method=Krotov
);
#-
fig = plot_amplitude(
    substitute(a, IdDict(a.control => opt_result_tanh.optimized_controls[1])),
    tlist
)
display(fig) #src
#md fig |> DisplayAs.PNG #hide
#-
amplitude = Array(substitute(a, IdDict(a.control => opt_result_tanh.optimized_controls[1]))) #src
@test minimum(amplitude) > -0.5 #src
@test maximum(amplitude) < 0.5 #src
