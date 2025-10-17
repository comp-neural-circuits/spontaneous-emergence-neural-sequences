using Statistics

@with_kw mutable struct the_whole_system
    dt::Float = 1e-4
    E::NeuronArray
    I::NeuronArray
    EE::Network
    EI::Network
    IE::Network
    II::Network
    buffer_frames_x::Int = Int(round(1/dt))
    buffer_frames_w::Int = Int(round(100/dt))
end

function run_process!(system::the_whole_system, rcdx, rcdw, frame)
    t = frame * system.dt
    update_network!(system.EE, system.dt, t)
    update_network!(system.IE, system.dt, t)
    update_network!(system.EI, system.dt, t)
    update_network!(system.II, system.dt, t)
    update_neurons!(system.E, system.dt, t)
    update_neurons!(system.I, system.dt, t)
    MOD_aeif.take_record_x!(frame, system, rcdx, system.buffer_frames_x)
    MOD_aeif.take_record_w!(frame, system, rcdw, system.buffer_frames_w)
end

function initialize_system!(system::the_whole_system)
    initialize_network!(system.EE);
    initialize_network!(system.EI);
    initialize_network!(system.IE);
    initialize_network!(system.II);
    initialize_test_neurons!(system.E);
    initialize_test_neurons!(system.I);
end
