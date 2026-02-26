using Sunny, CairoMakie, JLD2, ImageFiltering, LinearAlgebra

# --- OOP Declarations ---

"""
    PlotConfig

Defines the structure for drawing plots. 
"""
struct PlotConfig
    E_ranges::Vector{Tuple{Float64, Float64}}
    H_range::StepRangeLen{Float64}
    K_raw_range::StepRangeLen{Float64}
    b1::Vector{Float64}
    b2::Vector{Float64}
    grid_x::Matrix{Float64}
    grid_y::Matrix{Float64}
    bz_x::Vector{Float64}
    bz_y::Vector{Float64}
    
    function PlotConfig()
        E_ranges = [(0.8, 1.4), (1.1, 1.7), (2.2, 2.7)]
        H_range = -1.0:0.01:1.0
        K_raw_range = -1.0:0.01:1.0
        
        b1 = [1.0, 0.0]
        b2 = [0.5, sqrt(3)/2]
        
        grid_x = [h * b1[1] + k * b2[1] for h in H_range, k in K_raw_range]
        grid_y = [h * b1[2] + k * b2[2] for h in H_range, k in K_raw_range]
        
        k_points = [(2*b1 - b2)/3, (b1+b2)/3, (2*b2-b1)/3, (b2-2*b1)/3, (-b1-b2)/3, (b1-2*b2)/3]
        push!(k_points, k_points[1])
        bz_x = [p[1] for p in k_points]
        bz_y = [p[2] for p in k_points]
        
        new(E_ranges, H_range, K_raw_range, b1, b2, grid_x, grid_y, bz_x, bz_y)
    end
end

"""
    IntensityData

Data wrapper to encapsulate loaded intensity metrics
"""
struct IntensityData
    tag::String
    data::Array{Float64, 2}
    energies::Vector{Float64}
end

# --- Methods ---

function apply_mirror_sym(img::Matrix{Float64})
    img_symx = (img + reverse(img, dims=1)) / 2.0
    img_symy = (img_symx + reverse(img_symx, dims=2)) / 2.0
    img_symorigin = (img_symy + reverse(reverse(img_symy, dims=1), dims=2)) / 2.0
    return img_symorigin
end

# Multi-dispatch logic is handled via function parameters based on Tag Strings or sub-types
function load_intensity_data(tag::String)::IntensityData
    if tag == "3Q"
        fname = "CoNbS_result_3Q_state.jld2"
        println("Loading 3Q state: $fname")
        res = load(fname)["res"]
        return IntensityData("3Q", Float64.(res.data), collect(res.energies))
        
    elseif tag == "1Q"
        combined_data = nothing
        fnames = ["CoNbS_result_1Q_state_Q1.jld2", 
                  "CoNbS_result_1Q_state_Q2.jld2", 
                  "CoNbS_result_1Q_state_Q3.jld2"]
        
        sample_energies = nothing
        for (i, f) in enumerate(fnames)
            if !isfile(f) error("File not found: $f") end
            println("Merging 1Q state ($i/3): $f")
            
            l = load(f)
            d = Float64.(l["res"].data)
            
            if isnothing(combined_data)
                combined_data = copy(d)
                sample_energies = l["res"].energies
            else
                combined_data .+= d  # 강도 동등하게 합산
            end
        end
        return IntensityData("1Q", combined_data, collect(sample_energies))
    else
        error("Unsupported tag: $tag")
    end
end

function plot_structure!(config::PlotConfig, item::IntensityData)
    println("==========================================")
    println("Processing $(item.tag) state...")
    
    for (E_min, E_max) in config.E_ranges
        E_indices = findall(e -> E_min <= e <= E_max, item.energies)
        integrated_intensity = sum(item.data[E_indices, :], dims=1) |> vec
        intensity_map = reshape(integrated_intensity, length(config.H_range), length(config.K_raw_range))

        # 3. Gaussian Broadening (Isotropic 보정 적용)
        base_sigma = 2.5
        sigma_spatial = (base_sigma, base_sigma / (sqrt(3)/2)) # Cartesian isotropic 보정
        smoothed = imfilter(intensity_map, Kernel.gaussian(sigma_spatial))
        
        smoothed_map_final = smoothed; # to check comparison with/without symmetry
        # smoothed_map_final = apply_mirror_sym(smoothed)
        
        fig = Figure(size=(800, 700))
        ax = Axis(fig[1, 1], title = "CoNbS $(item.tag) state ($(E_min)-$(E_max) meV)", 
                  aspect = DataAspect(), xlabel="Qx", ylabel="Qy")
        
        K_range = √3/2 * (-1.0:0.01:1.0)
        
        # Color Scale setting
        max_val = maximum(smoothed_map_final)
        cmax = max(max_val * 0.8, 1.0) # avoid getting 0
        
        hm = surface!(ax, config.H_range, K_range, zeros(size(config.grid_x)); 
                      color = smoothed_map_final, colormap = :magma, 
                      shading = NoShading, colorrange = (0, cmax))
                      
        for n in -2:2, m in -2:2
            shift = n*config.b1 + m*config.b2
            dist = sqrt(n^2 + m^2 + n*m)
            l_alpha = dist == 0 ? 1.0 : (dist <= 1.1 ? 0.5 : 0.2)
            lines!(ax, config.bz_x .+ shift[1], config.bz_y .+ shift[2], color = (:white, l_alpha), linewidth = 2)
        end
        xlims!(ax, -0.9, 0.9); ylims!(ax, -0.8, 0.8)
        Colorbar(fig[1, 2], hm)
        
        # 6. 저장
        out_name = "CoNbS_$(item.tag)_Integrated_$(E_min)_$(E_max).png"
        out_name = replace(replace(out_name, "." => "_"), "_png" => ".png")
        save(out_name, fig, px_per_unit = 2.0)
        println("Saved: $out_name")
    end
end

function main()
    cfg = PlotConfig()
    
    for tag in ["1Q", "3Q"]
        data_item = load_intensity_data(tag)
        plot_structure!(cfg, data_item)
    end
end

main()
