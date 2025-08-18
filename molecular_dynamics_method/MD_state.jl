# Molecular Dynamics method
# 2-dim, rigid body molecules 
# solid state, liquid state or gas state 

using Plots 
using Random
using Statistics

mutable struct Molecular_Dynamics
    # scale
    xL::Float64          # size
    yL::Float64
    time::Float64
    iter::Int
    dt::Float64           

    # physical quantities of molecules
    x           # position: -L/2 < x < L/2 
    y           # position: -L/2 < y < L/2 
    v_x         # velocity of molecules
    v_y         

    diameter    # diameter of molecules

    # thermodynamical parameters
    T::Float64              # temperature
    N::Int                  # total number of molecules
    density::Float64        # number density of molecules
    area_fraction::Float64  

    # collision data 
    collision_data 
    collision_molecules

    # subroutines
    _new_collision_data

    function _init_size_and_position(N, density; ϵ=0.01)
        # set molecules at close-packed lattice points
        lattice_constant = sqrt(2.0 / sqrt(3.0) / density)
        scale = Int(round(sqrt(N / 4.0)))
        lattice_x = sqrt(3.0) * lattice_constant
        lattice_y = 2.0 * lattice_constant 
        xL = lattice_x * scale
        yL = lattice_y * scale 

        x = []
        y = []
        for layer = 1 : 4
            if layer % 2 == 1
                x0 = ϵ 
                y0 = ϵ + ((layer - 1) * 0.5) * lattice_constant
            elseif layer % 2 == 0
                x0 = ϵ + 0.5 * lattice_x
                y0 = ϵ + ((layer - 1) * 0.5) * lattice_constant
            end
            
            for j = 1 : scale
                y_temp = lattice_y * (j - 1) + y0
                if y_temp >= yL
                    break
                else
                    for i = 1 : scale
                        x_temp = lattice_x * (i - 1) + x0
                        if x_temp >= xL
                            break 
                        else
                            push!(x, x_temp)
                            push!(y, y_temp)
                        end
                    end
                end
            end
        end

        return xL, yL, x, y
    end

    function _init_velocity(N, T; criterion=3.5)
        # set v_x and v_y according to Maxwell-Boltzmann distribution via Box-Muller algorithm
        v_x = []
        v_y = []
        i = 1
        while i <= N
            box_muller_1 = sqrt(2.0 * T * (-log(rand(Float64))))
            box_muller_2 = 2.0 * π * rand(Float64)
            v_x_candidate = box_muller_1 * cos(box_muller_2)
            v_y_candidate = box_muller_1 * sin(box_muller_2)

            # reject too fast molecules 
            if (v_x_candidate ^ 2 + v_y_candidate ^ 2) > T * criterion ^ 2 
                continue 
            else
                push!(v_x, v_x_candidate)
                push!(v_y, v_y_candidate)
                i += 1
            end
        end

        # correct total momentum to 0
        momentum_x = mean(v_x)
        momentum_y = mean(v_y)
        for i = 1 : N
            v_x[i] -= momentum_x
            v_y[i] -= momentum_y
        end

        # scale velocity to keep formula: T = 1/2 * v^2
        tmp = 0.0
        for i = N 
            tmp += v_x[i] ^ 2 + v_y[i] ^ 2
        end
        velocity_temperature = tmp / (2.0 * N)

        scale = sqrt(T / velocity_temperature)
        v_x = v_x .* scale 
        v_y = v_y .* scale 

        return v_x, v_y
    end

    function _new_collision_data(N, diameter, i, xL, yL, x, y, v_x, v_y)
        time = Inf64
        partner = N
        for j = 1 : N 
            if j == i 
                continue 
            end

            distance_x = x[i] - x[j]
            distance_y = y[i] - y[j]
            relative_velocity_x = v_x[i] - v_x[j]
            relative_velocity_y = v_y[i] - v_y[j]
            # cyclic boundary condition
            distance_x -= round(distance_x / xL) * xL
            distance_y -= round(distance_y / yL) * yL

            inner_product = distance_x * relative_velocity_x + distance_y * relative_velocity_y
            distance_squared = distance_x ^ 2 + distance_y ^ 2
            relative_velocity_squared = relative_velocity_x ^ 2 + relative_velocity_y ^ 2
            if inner_product < 0.0
                discriminator = inner_product ^ 2 - relative_velocity_squared * (distance_squared - diameter ^ 2)
                if discriminator > 0.0
                    collision_time_ij = (-inner_product - sqrt(discriminator)) / relative_velocity_squared
                    if collision_time_ij < time 
                        time = collision_time_ij
                        partner = j 
                    end
                end
            end
        end

        return time, partner
    end

    function _init_collision_data(N, diameter, xL, yL, x, y, v_x, v_y)
        collision_data = []
        for i = 1 : N 
            collision_time, collision_partner = _new_collision_data(N, diameter, i, xL, yL, x, y, v_x, v_y)
            push!(collision_data, [collision_time, collision_partner])
        end

        return collision_data
    end

    function Molecular_Dynamics(T, N, area_fraction; diameter=1.0)
        density = area_fraction * 4.0 / π
        iter = 0
        time = 0.0
        dt = 0.0

        xL, yL, x, y = _init_size_and_position(N, density)
        v_x, v_y = _init_velocity(N, T)
        collision_data = _init_collision_data(N, diameter, xL, yL, x, y, v_x, v_y)
        collision_molecules = (0, 0)

        new(xL, yL, time, iter, dt, x, y, v_x, v_y, diameter, T, N, density, area_fraction, collision_data, collision_molecules, _new_collision_data)        
    end
end

function show_parameters(md)
    println("size_x: $(md.xL), size_y: $(md.yL)")
    println("molecular diameter: $(md.diameter)")
    println("thermodynamical parameters:")
    println("temperature: $(md.T), # of molecules: $(md.N), density: $(md.density), area_fraction: $(md.area_fraction)")
end

function solve!(md, iter_max)
    show_parameters(md)
    data = init_data(md, iter_max)

    while md.iter <= iter_max
        dt_and_collision_molecules_update!(md)
        time_evolution!(md)
        velocity_update!(md)
        collision_data_update!(md)
        save_position_data!(md, data)
    end
    
    plot_data(md, data)
    return data
end

function init_data(md, iter_max)
    data = zeros(md.N, iter_max+2, 2)
    for i = 1 : md.N 
        data[i, 1, 1] = md.x[i]
        data[i, 1, 2] = md.y[i]
    end
    return data
end

function dt_and_collision_molecules_update!(md)
    dt = Inf64
    molecule = 0
    for i = 1 : md.N 
        if md.collision_data[i][1] < dt 
            dt = md.collision_data[i][1] 
            molecule = i 
        end
    end

    md.dt = dt
    partner = Int(round(md.collision_data[molecule][2]))
    md.collision_molecules = (molecule, partner) 
end

function time_evolution!(md)
    md.iter += 1
    md.time += md.dt
    for i = 1 : md.N 
        md.collision_data[i][1] -= md.dt 

        x_temp = md.x[i] + md.v_x[i] * md.dt
        y_temp = md.y[i] + md.v_y[i] * md.dt
        x_temp -= round(x_temp / md.xL - 0.5) * md.xL       # cyclic boundary condition 
        y_temp -= round(y_temp / md.yL - 0.5) * md.yL
        md.x[i] = x_temp
        md.y[i] = y_temp
    end
end

function velocity_update!(md)
    (mol_1, mol_2) = md.collision_molecules

    distance_x = md.x[mol_1] - md.x[mol_2]
    distance_y = md.y[mol_1] - md.y[mol_2]
    distance_x -= round(distance_x / md.xL) * md.xL         # cyclic boundary condition
    distance_y -= round(distance_y / md.yL) * md.yL 
    distance = sqrt(distance_x ^ 2 + distance_y ^ 2)
    relative_velocity_x = md.v_x[mol_1] - md.v_x[mol_2]
    relative_velocity_y = md.v_y[mol_1] - md.v_y[mol_2]
    factor = (distance_x * relative_velocity_x + distance_y * relative_velocity_y) / md.diameter
    Δv_x = factor * distance_x / distance
    Δv_y = factor * distance_y / distance

    md.v_x[mol_1] -= Δv_x
    md.v_x[mol_2] += Δv_x
    md.v_y[mol_1] -= Δv_y
    md.v_y[mol_2] += Δv_y
end

function collision_data_update!(md)
    (mol_1, mol_2) = md.collision_molecules
    _collision_data_update!(md, mol_1)
    _collision_data_update!(md, mol_2)

    for i = 1 : md.N 
        partner = Int(round(md.collision_data[i][2]))
        if (partner == mol_1) || (partner == mol_2)
            _collision_data_update!(md, i)
        end
    end
end

function _collision_data_update!(md, molecule)
    time, partner = md._new_collision_data(md.N, md.diameter, molecule, md.xL, md.yL, md.x, md.y, md.v_x, md.v_y)
    md.collision_data[molecule] = [time, partner]
end

function save_position_data!(md, data)
    for i = 1 : md.N 
        data[i, md.iter+1, 1] = md.x[i]
        data[i, md.iter+1, 2] = md.y[i]
    end
end

function plot_data(md, data)
    plt = scatter()
    for i = 1 : md.N 
        scatter!(plt, data[i, :, 1], data[i, :, 2], size=(1050, 900), markeralpha=0.7, markerstrokewidth=0)
    end
    plot(plt)
end