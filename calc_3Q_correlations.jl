# 3Q Generation Refactored using CoreSimulation OOP and multiple threads
include("src/CoreSimulation.jl")
using .CoreSimulation
using Sunny

function main()
    # 1. Hamiltonian Parameters: J1-J2-Jc1 Minimal Model
    J1  = 1.9        # meV (Nearest neighbor)
    Jc1 = 0.5 * J1   # meV (Interlayer AF)
    J2  = 0.3 * J1   # meV (Long-range intra-layer)
    K_bq = 0.06 * J1 # meV (Biquadratic interaction) - set to stabilize 3Q
    
    params = SimulationParams(J1, Jc1, J2, K_bq)
    config = SystemConfig((20, 20, 6), 2000, 5, 5.0)

    H_range = -1.0:0.01:1.0
    K_raw_range = -1.0:0.01:1.0 
    qs = [ [0, 0, 1] + H * [1, 0, 0] + K * [-0.5, 1.0, 0.0] for H = H_range, K = K_raw_range ]

    runner = SimulationRunner(params, config)
    
    find_3Q_domain(runner)
    
    fname = "CoNbS_result_3Q_state.jld2"
    # Spawns parallel calculations of intensity for 5 samples
    run_parallel_simulation!(runner, qs, fname)
    
    println("3Q domain generation and saving complete.")
end

main()
