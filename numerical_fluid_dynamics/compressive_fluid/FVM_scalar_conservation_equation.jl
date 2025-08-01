# finite volume method
# scalar conservation equation: ∂u/∂t + ∂f(u)/∂x = 0

#=
 |----------|----------|---        ---|----------|---       ---|----------|----------|
    cell_1     cell_2                    cell_i                   cell_N-1   cell_N
x[1]       x[2]       x[3]           x[i]       x[i+1]        x[N-1]     x[N]       x[N+1]
                       Xl                                      Xr

     u[1]       u[2]                      u[i]                     u[N-1]     u[N]
flux
=#

mutable struct FVM
    f
    df
    N::Int
    Xl::Float64
    Xr::Float64
    dx::Float64
    dt::Float64

    x               # coordinate of boundary
    u               # mean value of u
    u_init          # initial value of u
    u_left          # value of left boundary of u
    u_right         # value of right boundary of u
    flux            #
    t::Float64      # time
    n::Int          # time step

    # constructor
    function FVM(f, df, init; N, Xl, Xr)
        dx = (Xr - Xl) / (N - 4)
        dt = 0.2 * dx
        x = [Xl + i * dx for i = -2 : N - 3]
        u_init = [init(x[i]) for i = 1 : N]
        u = copy(u_init)
        u_left = zeros(N+1)
        u_right = zeros(N+1)
        flux = zeros(N+1)
        t = 0.0
        n = 0
        new(f, df, N, Xl, Xr, dx, dt, x, u, u_init, u_left, u_right, flux, t, n)
    end
end

# methods
function show_parameters(fvm)
    println("N = $(fvm.N)\n[Xl, Xr] = [$(fvm.Xl), $(fvm.Xr)]\ndx = $(fvm.dx), dt = $(fvm.dt)\nt = $(fvm.t), steps = $(fvm.n)")
end

function show_graph(fvm)
    plot(fvm.x, fvm.u, xlabel="x", label="u")
    plot!(fvm.x, fvm.u_init, xlabel="x", label="u_init")
end

# spacial reconstruction methods
function reconstruction_linear!(fvm)
    for i = 2 : fvm.N - 2
        fvm.u_left[i+1] = fvm.u[i]
        fvm.u_right[i+1] = fvm.u[i+1]
    end
end

function reconstruction_Lax_Wendroff!(fvm)
    for i = 2 : fvm.N - 2 
        δ = 0.5 * (fvm.u[i+1] - fvm.u[i])
        v = _velocity(fvm, i)
        fvm.u_left[i+1] = fvm.u[i] + (1.0 - v * fvm.dt / fvm.dx) * δ
        fvm.u_right[i+1] = fvm.u[i+1] - (1.0 + v * fvm.dt / fvm.dx) * δ
    end
end

function reconstruction_Beam_Warming!(fvm)
    for i = 2 : fvm.N - 2 
        v = _velocity(fvm, i)
        fvm.u_left[i+1] = fvm.u[i] + (1.0 - v * fvm.dt / fvm.dx) * 0.5 * (fvm.u[i] - fvm.u[i-1])
        fvm.u_right[i+1] = fvm.u[i+1] - (1.0 + v * fvm.dt / fvm.dx) * 0.5 * (fvm.u[i+2] - fvm.u[i+1])
    end
end

function reconstruction_Fromm!(fvm)
    for i = 2 : fvm.N - 2
        v = _velocity(fvm, i)
        fvm.u_left[i+1] = fvm.u[i] + (1.0 - v * fvm.dt / fvm.dx) * 0.25 * (fvm.u[i+1] - fvm.u[i-1])
        fvm.u_right[i+1] = fvm.u[i+1] - (1.0 + v * fvm.dt / fvm.dx) * 0.25 * (fvm.u[i+2] - fvm.u[i])
    end
end

# TVD spacial reconstruction methods 
function _rL(fvm, i)
    if abs(fvm.u[i+1] - fvm.u[i]) < 0.00001
        return fvm.df(fvm.u[i])
    else
        return (fvm.u[i] - fvm.u[i-1]) / (fvm.u[i+1] - fvm.u[i])
    end
end

function _rR(fvm, i)
    if abs(fvm.u[i+1] - fvm.u[i]) < 0.00001
        return fvm.df(fvm.u[i+1])
    else
        return (fvm.u[i+2] - fvm.u[i+1]) / (fvm.u[i+1] - fvm.u[i])
    end
end

function reconstruction_minmod!(fvm)
    for i = 2 : fvm.N - 2
        δ = 0.5 * (fvm.u[i+1] - fvm.u[i])
        minmod = (x -> max(0.0, min(1.0, x)))
        v = _velocity(fvm, i)
        fvm.u_left[i+1] = fvm.u[i] + (1.0 - v * fvm.dt / fvm.dx) * δ * minmod(_rL(fvm, i))
        fvm.u_right[i+1] = fvm.u[i+1] - (1.0 + v * fvm.dt / fvm.dx) * δ * minmod(_rR(fvm, i))
    end
end

function reconstruction_superbee!(fvm)
    for i = 2 : fvm.N - 2
        δ = 0.5 * (fvm.u[i+1] - fvm.u[i])
        superbee = (x -> max(0.0, min(1.0, 2.0 * x), min(2.0, x)))
        v = _velocity(fvm, i)
        fvm.u_left[i+1] = fvm.u[i] + (1.0 - v * fvm.dt / fvm.dx) * δ * superbee(_rL(fvm, i))
        fvm.u_right[i+1] = fvm.u[i+1] - (1.0 + v * fvm.dt / fvm.dx) * δ * superbee(_rR(fvm, i))
    end
end

function reconstruction_van_Leer!(fvm)
    for i = 2 : fvm.N - 2
        δ = 0.5 * (fvm.u[i+1] - fvm.u[i])
        van_Leer = (x -> (x + abs(x)) / (1.0 + abs(x)))
        v = _velocity(fvm, i)
        fvm.u_left[i+1] = fvm.u[i] + (1.0 - v * fvm.dt / fvm.dx) * δ * van_Leer(_rL(fvm, i))
        fvm.u_right[i+1] = fvm.u[i+1] - (1.0 + v * fvm.dt / fvm.dx) * δ * van_Leer(_rR(fvm, i))
    end
end

function reconstruction_van_Albada!(fvm)
    for i = 2 : fvm.N - 2
        δ = 0.5 * (fvm.u[i+1] - fvm.u[i])
        van_Albada = (x -> (x + x^2) / (1.0 + x^2))
        v = _velocity(fvm, i)
        fvm.u_left[i+1] = fvm.u[i] + (1.0 - v * fvm.dt / fvm.dx) * δ * van_Albada(_rL(fvm, i))
        fvm.u_right[i+1] = fvm.u[i+1] - (1.0 + v * fvm.dt / fvm.dx) * δ * van_Albada(_rR(fvm, i))
    end
end

# Riemann solver
function _velocity(fvm, i)
    if abs(fvm.u_right[i] - fvm.u_left[i]) < 0.000001
        return fvm.df(fvm.u[i])
    else
        return (fvm.f(fvm.u_right[i]) - fvm.f(fvm.u_left[i])) / (fvm.u_right[i] - fvm.u_left[i])
    end
end

function Roe!(fvm)
    for i = 3 : fvm.N - 1
        fvm.flux[i] = 0.5 * (fvm.f(fvm.u_right[i]) + fvm.f(fvm.u_left[i]) - abs(_velocity(fvm, i)) * (fvm.u_right[i] - fvm.u_left[i]))
    end
end

function Local_Lax_Friedrichs!(fvm)
    for i = 3 : fvm.N - 1
        velocity = max(abs(fvm.df(fvm.u_right[i])), abs(fvm.df(fvm.u_left[i])))
        fvm.flux[i] = 0.5 * (fvm.f(fvm.u_right[i]) + fvm.f(fvm.u_left[i]) - velocity * (fvm.u_right[i] - fvm.u_left[i]))
    end
end

function Harten!(fvm; eps=0.25)      # eps: entropy modification: 0 < eps < 0.5
    for i = 3 : fvm.N - 1
        velocity = abs(_velocity(fvm, i))
        if velocity < 2.0 * eps 
            velocity = (velocity * velocity / (4.0 * eps)) + eps 
        end
        fvm.flux[i] = 0.5 * (fvm.f(fvm.u_right[i]) + fvm.f(fvm.u_left[i]) - velocity * (fvm.u_right[i] - fvm.u_left[i]))
    end
end

# Lax_Wendroff_Richtmyer method 
function Richtmyer!(fvm)
    for i = 3 : fvm.N 
        fvm.flux[i] = fvm.f(0.5 * ((fvm.u[i] + fvm.u[i-1]) - (fvm.u[i] - fvm.u[i-1]) * fvm.dt / fvm.dx))
    end
end

# time evolution
function update!(fvm)
    for i = 3 : fvm.N - 2
        fvm.u[i] -= fvm.dt / fvm.dx * (fvm.flux[i+1] - fvm.flux[i])
    end
end

# cyclic boundary condition
function cyclic_boundary_condition!(fvm)
    fvm.u[1] = fvm.u[fvm.N - 3]
    fvm.u[2] = fvm.u[fvm.N - 2]
    fvm.u[fvm.N - 1] = fvm.u[3]
    fvm.u[fvm.N] = fvm.u[4]
end

# solver selection
function _reconstructor_selection(rec, solver)
    if solver == "Richtmyer"
        return (fvm -> fvm)
    elseif solver == "Roe" || solver == "LLF" || solver == "Harten"
        if rec == "linear"
            return reconstruction_linear!
        elseif rec == "Lax_Wendroff"
            return reconstruction_Lax_Wendroff!
        elseif rec == "Beam_Warming"
            return reconstruction_Beam_Warming!
        elseif rec == "Fromm"
            return reconstruction_Fromm!
        elseif rec == "minmod"
            return reconstruction_minmod!
        elseif rec == "superbee"
            return reconstruction_superbee!
        elseif rec == "van_Leer"
            return reconstruction_van_Leer!
        elseif rec == "van_Albada"
            return reconstruction_van_Albada!
        else
            error("invalid reconstruction method: $rec")
        end
    else
        error("invalid solver: $solver")
    end
end

function _solver_selection(solver)
    if solver == "Richtmyer"
        return Richtmyer! 
    elseif solver == "Roe"
        return Roe!
    elseif solver == "LLF"
        return Local_Lax_Friedrichs!
    elseif solver == "Harten"
        return Harten!
    else 
        error("invalid solver: $solver")
    end
end

# main loop
function solve!(fvm, t_stop; rec="linear", solver="Roe")
    cyclic_boundary_condition!(fvm)
    spacial_reconstructor = _reconstructor_selection(rec, solver)
    Riemann_solver = _solver_selection(solver)

    while fvm.t <= t_stop
        fvm.t += fvm.dt
        fvm.n += 1
        spacial_reconstructor(fvm)
        Riemann_solver(fvm)
        update!(fvm)
        cyclic_boundary_condition!(fvm)
    end
    show_parameters(fvm)
    return fvm.x, fvm.u
end