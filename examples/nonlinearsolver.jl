using ModelingToolkit, NonlinearSolve

@variables x y z
@parameters σ ρ β

# Define a nonlinear system
eqs = [0 ~ σ * (y - x), 0 ~ x * (ρ - z) - y, 0 ~ x * y - β * z]
@mtkbuild ns = NonlinearSystem(eqs, [x, y, z], [σ, ρ, β])

u0 = [x => 1.0, y => 0.0, z => 0.0]

ps = [σ => 10.0
      ρ => 26.0
      β => 8 / 3]

prob = NonlinearProblem(ns, u0, ps)
sol = solve(prob, NewtonRaphson())