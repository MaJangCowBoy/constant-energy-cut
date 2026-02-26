module CoreSimulation

using Sunny, Statistics, LinearAlgebra, FFTW
using ProgressBars, JLD2
using Base.Threads
using Random

export SimulationParams, SystemConfig, SimulationRunner
export setup_system, find_1Q_domain, find_3Q_domain
export run_parallel_simulation!

# --- Data Structures (OOP equivalent) ---

"""
    SimulationParams

Data struct holding physical parameters for the Hamiltonian
"""
struct SimulationParams
    J1::Float64
    Jc1::Float64
    J2::Float64
    K_bq::Float64
end

"""
    SystemConfig

Data struct configuring the system properties
"""
struct SystemConfig
    L_size::Tuple{Int, Int, Int}
    n_thermal::Int
    n_samples::Int
    kT::Float64
end

"""
    SimulationRunner

Main struct that encapsulates the initial System and the generated SampledCorrelations
"""
mutable struct SimulationRunner
    params::SimulationParams
    config::SystemConfig
    crystal::Crystal
    base_sys::System
    
    function SimulationRunner(params::SimulationParams, config::SystemConfig)
        CoNbS = Crystal("CoTaS.cif", symprec = 1e-3)
        crystal = subcrystal(CoNbS, "Co1")
        
        # Init simple 1x1x1 system
        sys = System(crystal, [1 => Moment(s=3/2, g=2)], :dipole)
        
        # Set Hamiltonian based on SimulationParams
        if abs(params.K_bq) > 10 * eps()
            set_pair_coupling!(sys, (Si, Sj) -> params.J1*(Si'*Sj) + params.K_bq*(Si'*Sj)^2, Bond(1, 1, [1, 0, 0]))
        else
            set_exchange!(sys, params.J1, Bond(1, 1, [1, 0, 0]))
        end
        set_exchange!(sys, params.Jc1, Bond(1, 2, [0, 0, 0])) 
        set_exchange!(sys, params.J2, Bond(1, 1, [1, 2, 0]))
        
        new(params, config, crystal, sys)
    end
end

# --- Methods acting on Objects ---

"""
    find_1Q_domain(runner::SimulationRunner, target::String)

Deterministically constructs the target Q-domain states mapped onto the specific lattice layout.
The structural mapping constraint specifies that Spin(x+1, y+1, z, 2) = Spin(x, y, z, 1),
which simplifies into Spin(x, y, z, 2) = Spin(x-1, y-1, z, 1).
"""
function find_1Q_domain(runner::SimulationRunner, target::String)
    println("\n============================================")
    println("▶ Generating deterministic domain: ", target)
    println("============================================")
    
    # Directly build the desired simulation size system
    runner.base_sys = repeat_periodically(runner.base_sys, runner.config.L_size)
    
    Lx, Ly, Lz = runner.config.L_size
    S0 = [3/2, 0.0, 0.0] # Initial moment aligned on x-axis
    
    for z in 1:Lz, y in 1:Ly, x in 1:Lx
        # Based on defined properties for wavevectors:
        # Q1 spans x, Q2 spans y, Q3 spans x+y
        if target == "Q1"
            s1 = S0 * (-1)^x
            s2 = -s1
        elseif target == "Q2"
            s1 = S0 * (-1)^y
            s2 = -s1
        elseif target == "Q3"
            s1 = S0 * (-1)^(x+y)
            s2 = -s1
        else
            error("Unknown target $target")
        end
        
        runner.base_sys.dipoles[x,y,z,1] = s1;
        runner.base_sys.dipoles[x,y,z,2] = s2;
    end
    
    # Let the exact system relax energetically based on the assigned starting configurations
    for _ in 1:100  minimize_energy!(runner.base_sys) end
    
    println("✔ Successfully generated and minimized deterministic domain $target!")
end

"""
    find_3Q_domain(runner::SimulationRunner)

Creates the 3Q domain using random initialization.
"""
function find_3Q_domain(runner::SimulationRunner)
    println("\n============================================")
    println("▶ Generating 3Q Domain")
    println("============================================")
    
    sys_small = repeat_periodically(runner.base_sys, (2, 2, 1))
    
    randomize_spins!(sys_small)
    for _ in 1:40  randomize_spins!(sys_small) end
    for _ in 1:60  minimize_energy!(sys_small) end
    
    println("Domain minimized.")
    runner.base_sys = repeat_periodically(sys_small, runner.config.L_size)
end

"""
    run_parallel_simulation!(runner::SimulationRunner, qs, fname::String)

Executes parallel SampledCorrelations over the configured number of samples.
It utilizes multiple threads. Each thread copies the base system, thermalizes, and calculates intensities.
The results are averaged and saved at the end.
"""
function run_parallel_simulation!(runner::SimulationRunner, qs, fname::String)
    println("\n============================================")
    println("▶ Starting Parallel Simulation")
    println("  Available Threads: ", Threads.nthreads())
    println("============================================")
    
    # 1. Base Thermalization of the global system
    dt = 0.1 / (2 * runner.params.J1)
    langevin_base = Langevin(dt; damping=0.1, kT=runner.config.kT * meV_per_K)
    
    println("Initial base thermalization...")
    # NOTE: using simple progress bar logic because ProgressBar package inside threads might glitch
    for i in 1:(runner.config.n_thermal)
        step!(runner.base_sys, langevin_base)
    end
    
    dt_measure = 2 * langevin_base.dt
    energies = range(0, 20, 101)
    formfactors = [1 => FormFactor("Co1"; g_lande=2)]
    
    # Pre-allocate array for thread results
    n_samples = runner.config.n_samples
    intensities_results = Vector{Any}(undef, n_samples)
    sc_results = Vector{SampledCorrelations}(undef, n_samples)
    
    println("Sampling ", n_samples, " times in parallel...")
    
    # Note: Threads.@threads provides simple multi-threading in Julia
    Threads.@threads for i in 1:n_samples
        println("  Thread $(Threads.threadid()): Starting Sample $i/$n_samples")
        
        # Clone the completely thermalized system to branch off a new uncorrelated trajectory
        local_sys = clone_system(runner.base_sys)
        
        # Inject an independent random seed to ensure Langevin dynamics explores a different path 
        # for each thread. Uses current time, threadid, and sample index to guarantee uniqueness.
        Random.seed!(local_sys.rng, hash(time_ns()) + hash(Threads.threadid()) + i)
        
        local_langevin = Langevin(dt; damping=0.1, kT=runner.config.kT * meV_per_K)
        
        # Decorrelate
        for _ in 1:(runner.config.n_thermal)
            step!(local_sys, local_langevin)
        end
        
        # Setup Sc
        local_sc = SampledCorrelations(local_sys; dt=dt_measure, energies, measure=ssf_perp(local_sys; formfactors))
        add_sample!(local_sc, local_sys)
        
        # Compute Intensity for this sample and store
        # This keeps memory footprint bounded by collecting raw data right away rather than huge SC objects.
        local_res = intensities(local_sc, qs[:]; energies, kT=local_langevin.kT)
        intensities_results[i] = local_res
        sc_results[i] = local_sc
        
        println("  Thread $(Threads.threadid()): Finished Sample $i")
    end
    
    println("Aggregating parallel results...")
    
    # Merge SampledCorrelations objects securely
    merged_sc = merge_correlations(sc_results)
    
    # res.data contains complex or float matrices. Average them out.
    # Take base structure from the first task
    final_res = intensities_results[1]
    
    if n_samples > 1
        for i in 2:n_samples
            final_res.data .+= intensities_results[i].data
        end
        # Compute mean
        final_res.data ./= n_samples
    end
    
    println("Saving results to $fname")
    
    J1, Jc1, J2, K = runner.params.J1, runner.params.Jc1, runner.params.J2, runner.params.K_bq
    L_size = runner.config.L_size
    H_range = -1.0:0.01:1.0
    K_raw_range = -1.0:0.01:1.0 
    
    jldsave(fname; res=final_res, sc=merged_sc, L_size, qs, H_range, K_range=K_raw_range, J1, Jc1, J2, K)
    
    println("Save complete.")
end

end # module CoreSimulation
