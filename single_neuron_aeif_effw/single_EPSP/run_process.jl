using Statistics

@with_kw mutable struct the_whole_system
    dt::Float = 1e-4
    E::NeuronArray
end

function run_process!(system::the_whole_system, rcd, frame)
    t = frame * system.dt
    update_neurons!(system.E, system.dt, t)
    MOD_single_neuron_aeif.take_record!(frame, system, rcd)
end

function initialize_system!(system::the_whole_system)
    initialize_neurons!(system.E);
end

function initialize_test_system!(system::the_whole_system)
    initialize_neurons!(system.E);
end
