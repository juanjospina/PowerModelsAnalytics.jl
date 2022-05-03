"Plots branch impedances"
function plot_branch_impedance(data::Dict{String,Any}; branch_key::Any="branch", resistance_key::String="br_r", reactance_key::String="br_x")::Vega.VGSpec
    spec = deepcopy(default_branch_impedance_spec)

    scatter_data = [Dict("resistance" => sum(sum(branch["br_r"])), "reactance" => sum(sum(branch["br_x"])), "id" => id) for (id, branch) in data["branch"]]

    sort!(scatter_data; by=x -> parse(Int, x["id"]))

    @set! spec.data = []

    pushfirst!(spec.data, Dict("name" => "branch-impedances", "values" => scatter_data))

    resistance = Float64[x["resistance"] for x in scatter_data]
    reactance = Float64[x["reactance"] for x in scatter_data]

    for (field, _data) in zip(["resistance", "reactance"], [resistance, reactance])
        push!(spec.data, Dict(
            "name" => "binned-$field",
            "source" => "branch-impedances",
            "transform" => [
                Dict(
                    "type" => "bin", "field" => field,
                    "extent" => [minimum(_data), maximum(_data)],
                    "anchor" => mean(_data),
                    "step" => std(_data) / 2,
                    "nice" => true
                ),
                Dict(
                    "type" => "aggregate",
                    "key" => "bin0",
                    "groupby" => ["bin0", "bin1"],
                    "fields" => ["bin0"],
                    "ops" => ["count"],
                    "as" => ["count"]
                )
            ]
        )
        )
    end

    return spec
end


"""
    `plot_load_summary(file, result, case; kwargs...)`

    Plots total generation, total load served, and total forecasted load for a given `case` and `result`, saving to `file`

    Arguments:

    `file::String`: file path to saved figure
    `result::Dict{String,Any}`: multinetwork solution data (contains load statuses)
    `case::Dict{String,Any}`: Original case file (without calcuated loads) for forecasted loads
    `log::Bool`: If `true`, plots y-axis on log scale
    `intermediate::Bool`: If `true`, plots intermediate steps of plot (for animations).
    `legend_position::Symbol`: Position of legend, accepts the following symbols: `:right`, `:left`, `:top`, `:bottom`, `:inside`, `:best`, `:legend`, `:topright`, `:topleft`, `:bottomleft`, `:bottomright`
"""
function plot_load_summary(file::String, result::Dict{String,Any}, case::Dict{String,Any}; log::Bool=false, intermediate::Bool=false, legend_position::Symbol=:best)::Vega.VGSpec
    @assert Int(get(case, "data_model", 1)) == 1 && get(case, "per_unit", true) "This function only supports plotting MATHEMATICAL data models in per-unit representation"

    spec = Vega.loadvgspec("src/vega/load_summary.json")

    x = 0:length(result["nw"]) - 1
    generation = [x for (n, x) in sort([(parse(Int, n), sum(sum(_replace_nan(gen["pg"])) * nw["baseMVA"] for (i, gen) in nw["gen"])) for (n, nw) in result["nw"]]; by=x -> x[1])]
    storage = [x for (n, x) in sort([(parse(Int, n), sum(sum(_replace_nan(strg["ps"])) * nw["baseMVA"] for (i, strg) in nw["storage"])) for (n, nw) in result["nw"]]; by=x -> x[1])]
    total_generated = generation .+ storage
    total_load_served = [x for (n, x) in sort([(parse(Int, n), sum(sum(_replace_nan(load["status"] * case["nw"]["$n"]["load"]["$i"]["pd"])) * nw["baseMVA"] for (i, load) in nw["load"])) for (n, nw) in result["nw"]]; by=x -> x[1])]
    total_load_forecast = [x for (n, x) in sort([(parse(Int, n), sum(sum(_replace_nan(load["pd"])) * case["nw"]["$n"]["baseMVA"] for (i, load) in nw["load"])) for (n, nw) in case["nw"]]; by=x -> x[1])]

    max_digits = max_digits = maximum([length("$n") for n in x])
    @debug "" total_generated total_load_served total_load_forecast

    spec = deepcopy(default_total_source_demand_summary_spec)

    power_summary_data = [
        Dict(
            "x" => x[i],
            "y" => y[i],
            "c" => c - 1,
        ) for (c, y) in enumerate([total_generated, total_load_served, total_load_forecast]) for i in 1:length(x)
    ]

    @set! spec.data = [
        Dict(
            "name" => "table",
            "values" => power_summary_data,
            "transform" => [
                Dict(
                    "type" => "stack",
                    "groupby" => ["x"],
                    "sort" => Dict(
                        "field" => "c"
                    ),
                    "field" => "y"
                )
            ]
        )
    ]

    @set! spec.axes[2]["title"] = "Power (MW)"

    if log
        @set! spec.scales[2]["type"] = "log"
    end

    if intermediate
        _tmp_data = []
        for (i, _data) in enumerate(eachrow(reshape(power_summary_data, :, 3)))
            append!(_tmp_data, _data)
            @set! spec.data = [
                Dict(
                    "name" => "table",
                    "values" => _tmp_data,
                    "transform" => [
                        Dict(
                            "type" => "stack",
                            "groupby" => ["x"],
                            "sort" => Dict(
                                "field" => "c"
                            ),
                            "field" => "y"
                        )
                    ]
                )
            ]

            filename_parts = split(file, ".")
            filename = join(filename_parts[1:end-1], ".")
            ext = filename_parts[end]

            _fileout = "$(filename)_$(lpad(i, max_digits, "0")).$(ext)"

            Vega.save(_fileout, spec)
        end
    else
        Vega.save(file, spec)
    end

    return spec
end


"""
    `plot_source_demand_summary(file::String, mn_case::Dict{String,<:Any}; kwargs...)`

    Plots the total delivery from sources (generation) and total receipts by demands (load)

    Arguments:

    `fileout::String`: path to file where plot will be saved
    `mn_case::Dict{String,<:Any}`: a multinetwork case
    `yscale::Symbol`: To set log scale, `:log10`, else `:identity`
    `save_intermediate_frames::Bool`: if `true`, each frame of the multinetwork will be saved separately
    `legend_position::Symbol`: Position of legend, accepts the following symbols: `:right`, `:left`, `:top`, `:bottom`, `:inside`, `:best`, `:legend`, `:topright`, `:topleft`, `:bottomleft`, `:bottomright`
    `sources::Dict{String,<:Any}`: information about sources (e.g. generators)
    `demands::Dict{String,<:Any}`: information about demands (e.g. loads)
    `totals::Symbol`: Choose `:real`, `:imaginary`, `:complex`
"""
function plot_source_demand_summary(fileout::String, mn_case::Dict{String,<:Any};
    yscale::Symbol=:identity,
    save_intermediate_frames::Bool=false,
    legend_position::Symbol=:best,
    sources::Dict{String,<:Any}=default_sources_eng,
    demands::Dict{String,<:Any}=default_demands_eng,
    totals::Symbol=:real,
    )::Vega.VGSpec

    x = 1:length(mn_case["nw"])
    total_generated = Vector{Real}(undef, length(x))

    for (n, nw) in get(mn_case, "nw", Dict())
        value = Complex(0.0, 0.0)

        for (type, settings) in sources
            real_key = get(settings, "inactive_real", "" => 0)[1]
            imag_key = get(settings, "inactive_imaginary", "" => 0)[1]

            for (_,obj) in get(nw, type, Dict())
                if totals == :real
                    v_real = get(obj, real_key, 0.0)
                    v_imag = 0.0
                elseif totals == :imaginary
                    v_real = 0.0
                    v_imag = get(obj, imag_key, 0.0)
                else
                    v_real = get(obj, real_key, 0.0)
                    v_imag = get(obj, imag_key, 0.0)
                end

                value += sum(Complex.(v_real, v_imag))
            end
        end

        total_generated[parse(Int, n)] = norm(value)
    end

    total_demand_served = Vector{Real}(undef, length(x))
    total_demand_forecast = Vector{Real}(undef, length(x))
    for (n, nw) in get(mn_case, "nw", Dict())
        served_value = Complex(0.0, 0.0)
        forecast_value = Complex(0.0, 0.0)

        for (type, settings) in demands
            served_real_key = get(settings, "inactive_real", "" => 0)[1]
            served_imag_key = get(settings, "inactive_imaginary", "" => 0)[1]

            forecast_real_key = get(settings, "original_demand_real", "")
            forecast_imag_key = get(settings, "original_demand_imaginary", "")

            for (_,obj) in get(nw, type, Dict())
                if totals == :real
                    v_served_real = get(obj, served_real_key, 0.0)
                    v_served_imag = 0.0
                    v_forecast_real = get(obj, forecast_real_key, 0.0)
                    v_forecast_imag = 0.0
                elseif totals == :imaginary
                    v_served_real = 0.0
                    v_served_imag = get(obj, served_imag_key, 0.0)
                    v_forecast_real = 0.0
                    v_forecast_imag = get(obj, forecast_imag_key, 0.0)
                else
                    v_served_real = get(obj, served_real_key, 0.0)
                    v_served_imag = get(obj, served_imag_key, 0.0)
                    v_forecast_real = get(obj, forecast_real_key, 0.0)
                    v_forecast_imag = get(obj, forecast_imag_key, 0.0)
                end

                served_value += sum(Complex.(v_served_real, v_served_imag))
                forecast_value += sum(Complex.(v_forecast_real, v_forecast_imag))
            end
        end

        total_demand_served[parse(Int, n)] = norm(served_value)
        total_demand_forecast[parse(Int, n)] = norm(forecast_value)
    end

    max_digits = maximum([length(n) for (n,_) in mn_case["nw"]])

    power_scale_factor = mn_case["nw"]["1"]["settings"]["sbase"]
    units_str = power_scale_factor == 1.0 ? "W" : power_scale_factor == 1e3 ? "kW" : power_scale_factor == 1e6 ? "MW" : "$power_scale_factor W"

    spec = deepcopy(default_total_source_demand_summary_spec)

    power_summary_data = [
        Dict(
            "x" => x[i],
            "y" => y[i],
            "c" => c - 1,
        ) for (c, y) in enumerate([total_generated, total_demand_served, total_demand_forecast]) for i in 1:length(x)
    ]

    @set! spec.data = [
        Dict(
            "name" => "table",
            "values" => power_summary_data,
            "transform" => [
                Dict(
                    "type" => "stack",
                    "groupby" => ["x"],
                    "sort" => Dict("field" => "c"),
                    "field" => "y"
                )
            ]
        )
    ]

    @set! spec.axes[2]["title"] = "Power ($units_str)"

    if save_intermediate_frames
        _tmp_data = []
        for (i, _data) in enumerate(eachrow(reshape(power_summary_data, :, 3)))
            append!(_tmp_data, _data)
            @set! spec.data = [
                Dict(
                    "name" => "table",
                    "values" => _tmp_data,
                    "transform" => [
                        Dict(
                            "type" => "stack",
                            "groupby" => ["x"],
                            "sort" => Dict("field" => "c"),
                            "field" => "y"
                        )
                    ]
                )
            ]

            filename_parts = split(fileout, ".")
            filename = join(filename_parts[1:end-1], ".")
            ext = filename_parts[end]

            _fileout = "$(filename)_$(lpad(i, max_digits, "0")).$(ext)"

            Vega.save(_fileout, spec)
        end
    else
        Vega.save(fileout, spec)
    end

    return spec
end

# --------------------------------------------------------------------------------------

"""
    `function plot_distribution_total_source_demand_summary(
        file::String,
        mn_solution::Dict{String,<:Any};
        kwargs...
    )`

    Plots the total generation from sources and total load demands of the distribution system.

    Arguments:

    `fileout::String`: path to file where plot will be saved
    `mn_solution::Dict{String,<:Any}`: a multinetwork solution from PMD or PMITD.
    `yscale::Symbol`: To set log scale, `:log10`, else `:identity`
    `save_intermediate_frames::Bool`: if `true`, each frame of the multinetwork will be saved separately
    `legend_position::Symbol`: Position of legend, accepts the following symbols: `:right`, `:left`, `:top`, `:bottom`, `:inside`, `:best`, `:legend`, `:topright`, `:topleft`, `:bottomleft`, `:bottomright`
    `sources::Dict{String,<:Any}`: information about sources (e.g. generators)
    `demands::Dict{String,<:Any}`: information about demands (e.g. loads)
    `totals::Symbol`: Choose `:real`, `:imaginary`, `:complex`
"""
function plot_distribution_total_source_demand_summary(fileout::String, mn_solution::Dict{String,<:Any};
    yscale::Symbol=:identity,
    save_intermediate_frames::Bool=false,
    legend_position::Symbol=:best,
    sources::Dict{String,<:Any}=default_sources_eng,
    demands::Dict{String,<:Any}=default_demands_eng,
    totals::Symbol=:real,
    )::Vega.VGSpec

    # assign correct data based on multiinfrastructure condition
    if get(mn_solution, "multiinfrastructure", false) != false
        # Extract pmd solution only
        mn_solution = mn_solution["it"]["pmd"]
    end

    x = 1:length(mn_solution["nw"])
    total_generated = Vector{Real}(undef, length(x))

    for (n, nw) in get(mn_solution, "nw", Dict())
        value = Complex(0.0, 0.0)

        for (type, settings) in sources
            real_key = get(settings, "inactive_real", "" => 0)[1]
            imag_key = get(settings, "inactive_imaginary", "" => 0)[1]

            for (_,obj) in get(nw, type, Dict())
                if totals == :real
                    v_real = get(obj, real_key, 0.0)
                    v_imag = 0.0
                elseif totals == :imaginary
                    v_real = 0.0
                    v_imag = get(obj, imag_key, 0.0)
                else
                    v_real = get(obj, real_key, 0.0)
                    v_imag = get(obj, imag_key, 0.0)
                end

                value += sum(Complex.(v_real, v_imag))
            end
        end

        total_generated[parse(Int, n)] = norm(value)
    end

    total_demand_served = Vector{Real}(undef, length(x))
    for (n, nw) in get(mn_solution, "nw", Dict())
        served_value = Complex(0.0, 0.0)

        for (type, settings) in demands
            served_real_key = get(settings, "inactive_real", "" => 0)[1]
            served_imag_key = get(settings, "inactive_imaginary", "" => 0)[1]

            for (_,obj) in get(nw, type, Dict())
                if totals == :real
                    v_served_real = get(obj, served_real_key, 0.0)
                    v_served_imag = 0.0
                elseif totals == :imaginary
                    v_served_real = 0.0
                    v_served_imag = get(obj, served_imag_key, 0.0)
                else
                    v_served_real = get(obj, served_real_key, 0.0)
                    v_served_imag = get(obj, served_imag_key, 0.0)
                end

                served_value += sum(Complex.(v_served_real, v_served_imag))
            end
        end

        total_demand_served[parse(Int, n)] = norm(served_value)
    end

    max_digits = maximum([length(n) for (n,_) in mn_solution["nw"]])

    power_scale_factor = mn_solution["nw"]["1"]["settings"]["sbase"]

    if totals == :real
        units_str = power_scale_factor == 1e2 ? "W" : power_scale_factor == 1e5 ? "kW" : power_scale_factor == 1e10 ? "MW" : "$power_scale_factor W"
    elseif totals == :imaginary
        units_str = power_scale_factor == 1e2 ? "VAR" : power_scale_factor == 1e5 ? "kVAR" : power_scale_factor == 1e10 ? "MVAR" : "$power_scale_factor VAR"
    else
        units_str = power_scale_factor == 1e2 ? "VA" : power_scale_factor == 1e5 ? "kVA" : power_scale_factor == 1e10 ? "MVA" : "$power_scale_factor VA"
    end

    spec = deepcopy(default_total_source_demand_summary_spec)

    power_summary_data = [
        Dict(
            "x" => x[i],
            "y" => y[i],
            "c" => c - 1,
        ) for (c, y) in enumerate([total_generated, total_demand_served]) for i in 1:length(x)
    ]

    @set! spec.data = [
        Dict(
            "name" => "table",
            "values" => power_summary_data,
            "transform" => [
                Dict(
                    "type" => "stack",
                    "groupby" => ["x"],
                    "sort" => Dict("field" => "c"),
                    "field" => "y"
                )
            ]
        )
    ]

    @set! spec.axes[2]["title"] = "Power ($units_str)"

    if save_intermediate_frames
        _tmp_data = []
        for (i, _data) in enumerate(eachrow(reshape(power_summary_data, :, 3)))
            append!(_tmp_data, _data)
            @set! spec.data = [
                Dict(
                    "name" => "table",
                    "values" => _tmp_data,
                    "transform" => [
                        Dict(
                            "type" => "stack",
                            "groupby" => ["x"],
                            "sort" => Dict("field" => "c"),
                            "field" => "y"
                        )
                    ]
                )
            ]

            filename_parts = split(fileout, ".")
            filename = join(filename_parts[1:end-1], ".")
            ext = filename_parts[end]

            _fileout = "$(filename)_$(lpad(i, max_digits, "0")).$(ext)"

            Vega.save(_fileout, spec)
        end
    else
        Vega.save(fileout, spec)
    end

    return spec
end


"""
    `function plot_transmission_total_source_demand_summary(
        file::String,
        mn_solution::Dict{String,<:Any};
        kwargs...
    )`

    Plots the total generation from sources of the transmission system (load demands are not available from PM).

    Arguments:

    `fileout::String`: path to file where plot will be saved
    `mn_solution::Dict{String,<:Any}`: a multinetwork solution from PM or PMITD.
    `yscale::Symbol`: To set log scale, `:log10`, else `:identity`
    `save_intermediate_frames::Bool`: if `true`, each frame of the multinetwork will be saved separately
    `legend_position::Symbol`: Position of legend, accepts the following symbols: `:right`, `:left`, `:top`, `:bottom`, `:inside`, `:best`, `:legend`, `:topright`, `:topleft`, `:bottomleft`, `:bottomright`
    `sources::Dict{String,<:Any}`: information about sources (e.g. generators)
    `totals::Symbol`: Choose `:real`, `:imaginary`, `:complex`
"""
function plot_transmission_total_source_demand_summary(fileout::String, mn_solution::Dict{String,<:Any};
    yscale::Symbol=:identity,
    save_intermediate_frames::Bool=false,
    legend_position::Symbol=:best,
    sources::Dict{String,<:Any}=default_sources_math,
    totals::Symbol=:real,
    )::Vega.VGSpec

    # assign correct data based on multiinfrastructure condition
    if get(mn_solution, "multiinfrastructure", false) != false
        # Extract pm solution only
        mn_solution = mn_solution["it"]["pm"]
    end

    x = 1:length(mn_solution["nw"])
    total_generated = Vector{Real}(undef, length(x))

    for (n, nw) in get(mn_solution, "nw", Dict())
        value = Complex(0.0, 0.0)

        for (type, settings) in sources

            real_key = get(settings, "inactive_real", "" => 0)[1]
            imag_key = get(settings, "inactive_imaginary", "" => 0)[1]

            for (_,obj) in get(nw, type, Dict())
                if totals == :real
                    v_real = get(obj, real_key, 0.0)
                    v_imag = 0.0
                elseif totals == :imaginary
                    v_real = 0.0
                    v_imag = get(obj, imag_key, 0.0)
                else
                    v_real = get(obj, real_key, 0.0)
                    v_imag = get(obj, imag_key, 0.0)
                end

                value += sum(Complex.(v_real, v_imag))
            end
        end

        total_generated[parse(Int, n)] = norm(value)
    end

    max_digits = maximum([length(n) for (n,_) in mn_solution["nw"]])

    power_scale_factor = mn_solution["nw"]["1"]["baseMVA"]

    if totals == :real
        units_str = power_scale_factor == 1e-6 ? "W" : power_scale_factor == 0.001 ? "kW" : power_scale_factor == 100 ? "MW" : "$power_scale_factor W"
    elseif totals == :imaginary
        units_str = power_scale_factor == 1e-6 ? "VAR" : power_scale_factor == 0.001 ? "kVAR" : power_scale_factor == 100 ? "MVAR" : "$power_scale_factor VAR"
    else
        units_str = power_scale_factor == 1e-6 ? "VA" : power_scale_factor == 0.001 ? "kVA" : power_scale_factor == 100 ? "MVA" : "$power_scale_factor VA"
    end

    spec = deepcopy(default_total_source_demand_summary_spec)

    power_summary_data = [
        Dict(
            "x" => x[i],
            "y" => y[i],
            "c" => c - 1,
        ) for (c, y) in enumerate([total_generated]) for i in 1:length(x)
    ]

    @set! spec.data = [
        Dict(
            "name" => "table",
            "values" => power_summary_data,
            "transform" => [
                Dict(
                    "type" => "stack",
                    "groupby" => ["x"],
                    "sort" => Dict("field" => "c"),
                    "field" => "y"
                )
            ]
        )
    ]

    @set! spec.axes[2]["title"] = "Power ($units_str)"

    if save_intermediate_frames
        _tmp_data = []
        for (i, _data) in enumerate(eachrow(reshape(power_summary_data, :, 3)))
            append!(_tmp_data, _data)
            @set! spec.data = [
                Dict(
                    "name" => "table",
                    "values" => _tmp_data,
                    "transform" => [
                        Dict(
                            "type" => "stack",
                            "groupby" => ["x"],
                            "sort" => Dict("field" => "c"),
                            "field" => "y"
                        )
                    ]
                )
            ]

            filename_parts = split(fileout, ".")
            filename = join(filename_parts[1:end-1], ".")
            ext = filename_parts[end]

            _fileout = "$(filename)_$(lpad(i, max_digits, "0")).$(ext)"

            Vega.save(_fileout, spec)
        end
    else
        Vega.save(fileout, spec)
    end

    return spec
end


"""
    `function plot_distribution_source_demand_summary(
        file::String,
        mn_solution::Dict{String,<:Any};
        kwargs...
    )`

    Plots the individual generation from sources and load demands of the distribution system.

    Arguments:

    `fileout::String`: path to file where plot will be saved
    `mn_solution::Dict{String,<:Any}`: a multinetwork solution from PMD or PMITD.
    `yscale::Symbol`: To set log scale, `:log10`, else `:identity`
    `legend_position::Symbol`: Position of legend, accepts the following symbols: `:right`, `:left`, `:top`, `:bottom`, `:inside`, `:best`, `:legend`, `:topright`, `:topleft`, `:bottomleft`, `:bottomright`
    `sources::Dict{String,<:Any}`: information about sources (e.g. generators)
    `demands::Dict{String,<:Any}`: information about demands (e.g. loads)
    `totals::Symbol`: Choose `:real` or `:imaginary` (`:complex` not available for this function)
"""
function plot_distribution_source_demand_summary(fileout::String, mn_solution::Dict{String,<:Any};
    yscale::Symbol=:identity,
    legend_position::Symbol=:best,
    sources::Dict{String,<:Any}=default_sources_eng,
    demands::Dict{String,<:Any}=default_demands_eng,
    totals::Symbol=:real,
    )::Vega.VGSpec

    # assign correct data based on multiinfrastructure condition
    if get(mn_solution, "multiinfrastructure", false) != false
        # Extract pmd solution only
        mn_solution = mn_solution["it"]["pmd"]

        # Get total number of gens, solar sotrage
        gens_number = length(get(mn_solution["nw"]["1"], "generator", Dict()))
        gens_number += length(get(mn_solution["nw"]["1"], "solar", Dict()))
        gens_number += length(get(mn_solution["nw"]["1"], "storage", Dict()))

    else
        # Get total number of gens, solar sotrage + 'voltage_source (slack gen)'
        gens_number = length(get(mn_solution["nw"]["1"], "generator", Dict()))
        gens_number += length(get(mn_solution["nw"]["1"], "solar", Dict()))
        gens_number += length(get(mn_solution["nw"]["1"], "storage", Dict()))
        gens_number += length(get(mn_solution["nw"]["1"], "voltage_source", Dict()))
    end

    x = 1:length(mn_solution["nw"]) # number of nws
    gen_names = Vector{String}()  # for the labels of the graph

    # Create arrays to store values (depending on total type)
    if totals == :real
        active_gen = Array{Real}(undef, length(x), gens_number)
    elseif totals == :imaginary
        reactive_gen = Array{Real}(undef, length(x), gens_number)
    else
        active_gen = Array{Real}(undef, length(x), gens_number)
        reactive_gen = Array{Real}(undef, length(x), gens_number)
    end

    for (n, nw) in get(mn_solution, "nw", Dict())
        gen_counter = 0 # counts the number of gen. devices (gens, solar, storage).

        for (type, settings) in sources
            real_key = get(settings, "inactive_real", "" => 0)[1]
            imag_key = get(settings, "inactive_imaginary", "" => 0)[1]

            for (name,obj) in get(nw, type, Dict())
                gen_counter += 1

                if totals == :real
                    v_real = get(obj, real_key, 0.0)
                    v_imag = 0.0
                    active_gen[parse(Int, n),gen_counter] = sum(v_real)
                elseif totals == :imaginary
                    v_real = 0.0
                    v_imag = get(obj, imag_key, 0.0)
                    reactive_gen[parse(Int, n),gen_counter] = sum(v_imag)
                else
                    v_real = get(obj, real_key, 0.0)
                    v_imag = get(obj, imag_key, 0.0)
                    active_gen[parse(Int, n),gen_counter] = sum(v_real)
                    reactive_gen[parse(Int, n),gen_counter] = sum(v_imag)
                end

                # only in the first nw (save names of the gens)
                if (parse(Int, n) == 1)
                    push!(gen_names, name)
                end

            end
        end
    end

    load_names = Vector{String}()  # For the labels of the graph
    loads_number = length(get(mn_solution["nw"]["1"], "load", Dict())) # Get total number of loads

    # Create arrays to store values (depending on total type)
    if totals == :real
        active_load = Array{Real}(undef, length(x), loads_number)
    elseif totals == :imaginary
        reactive_load = Array{Real}(undef, length(x), loads_number)
    else
        active_load = Array{Real}(undef, length(x), loads_number)
        reactive_load = Array{Real}(undef, length(x), loads_number)
    end


    for (n, nw) in get(mn_solution, "nw", Dict())
        load_counter = 0 # counts the number of loads

        for (type, settings) in demands
            served_real_key = get(settings, "inactive_real", "" => 0)[1]
            served_imag_key = get(settings, "inactive_imaginary", "" => 0)[1]

            for (name,obj) in get(nw, type, Dict())
                load_counter += 1

                if totals == :real
                    v_served_real = get(obj, served_real_key, 0.0)
                    v_served_imag = 0.0
                    active_load[parse(Int, n),load_counter] = sum(v_served_real)
                elseif totals == :imaginary
                    v_served_real = 0.0
                    v_served_imag = get(obj, served_imag_key, 0.0)
                    reactive_load[parse(Int, n),load_counter] = sum(v_served_imag)
                else
                    v_served_real = get(obj, served_real_key, 0.0)
                    v_served_imag = get(obj, served_imag_key, 0.0)
                    active_load[parse(Int, n),load_counter] = sum(v_served_real)
                    reactive_load[parse(Int, n),load_counter] = sum(v_served_imag)
                end

                # only in the first nw (save names of the loads)
                if (parse(Int, n) == 1)
                    push!(load_names, name)
                end

            end
        end

    end

    # Assign active/reactive power for plotting
    if totals == :real
        gen_power = active_gen
        load_power =  active_load
    elseif totals == :imaginary
        gen_power = reactive_gen
        load_power =  reactive_load
    else
        error("Can't plot active and reactive powers at the same time for this function.")
    end

    power_scale_factor = mn_solution["nw"]["1"]["settings"]["sbase"]

    if totals == :real
        units_str = power_scale_factor == 1e2 ? "W" : power_scale_factor == 1e5 ? "kW" : power_scale_factor == 1e10 ? "MW" : "$power_scale_factor W"
    elseif totals == :imaginary
        units_str = power_scale_factor == 1e2 ? "VAR" : power_scale_factor == 1e5 ? "kVAR" : power_scale_factor == 1e10 ? "MVAR" : "$power_scale_factor VAR"
    else
        units_str = power_scale_factor == 1e2 ? "VA" : power_scale_factor == 1e5 ? "kVA" : power_scale_factor == 1e10 ? "MVA" : "$power_scale_factor VA"
    end

    spec = deepcopy(default_source_demand_summary_spec)

    power_summary_data = []
    for i in 1:length(x)
        elem_count = 0
        for (c, y) in enumerate([gen_power[i,:], load_power[i,:]])
            for (v,b) in enumerate(y)
                elem_count += 1
                push!(power_summary_data, Dict("x" => x[i], "y" => b, "c" => elem_count-1))
            end
        end
    end

    @set! spec.data = [
        Dict(
            "name" => "table",
            "values" => power_summary_data,
            "transform" => [
                Dict(
                    "type" => "stack",
                    "groupby" => ["x"],
                    "sort" => Dict("field" => "c"),
                    "field" => "y"
                )
            ]
        )
    ]

    @set! spec.scales[4]["range"] = [gen_names; load_names] # set legends
    @set! spec.axes[2]["title"] = "Power ($units_str)"      # set title

    Vega.save(fileout, spec)

    return spec
end


"""
    `function plot_distribution_source_perphase_summary(
        file::String,
        mn_solution::Dict{String,<:Any};
        kwargs...
    )`

    Plots the individual per-phase generation from sources of the distribution system.

    Arguments:

    `fileout::String`: path to file where plot will be saved
    `mn_solution::Dict{String,<:Any}`: a multinetwork solution from PMD or PMITD.
    `yscale::Symbol`: To set log scale, `:log10`, else `:identity`
    `legend_position::Symbol`: Position of legend, accepts the following symbols: `:right`, `:left`, `:top`, `:bottom`, `:inside`, `:best`, `:legend`, `:topright`, `:topleft`, `:bottomleft`, `:bottomright`
    `sources::Dict{String,<:Any}`: information about sources (e.g. generators)
    `totals::Symbol`: Choose `:real` or `:imaginary` (`:complex` not available for this function)
"""
function plot_distribution_source_perphase_summary(fileout::String, mn_solution::Dict{String,<:Any};
    yscale::Symbol=:identity,
    legend_position::Symbol=:best,
    sources::Dict{String,<:Any}=default_sources_eng,
    totals::Symbol=:real,
    )::Vega.VGSpec

    # assign correct data based on multiinfrastructure condition
    if get(mn_solution, "multiinfrastructure", false) != false
        # Extract pmd solution only
        mn_solution = mn_solution["it"]["pmd"]

        # Get total number of gens, solar sotrage
        gens_number = length(get(mn_solution["nw"]["1"], "generator", Dict()))
        gens_number += length(get(mn_solution["nw"]["1"], "solar", Dict()))
        gens_number += length(get(mn_solution["nw"]["1"], "storage", Dict()))

    else
        # Get total number of gens, solar sotrage + 'voltage_source (slack gen)'
        gens_number = length(get(mn_solution["nw"]["1"], "generator", Dict()))
        gens_number += length(get(mn_solution["nw"]["1"], "solar", Dict()))
        gens_number += length(get(mn_solution["nw"]["1"], "storage", Dict()))
        gens_number += length(get(mn_solution["nw"]["1"], "voltage_source", Dict()))
    end

    x = 1:length(mn_solution["nw"]) # number of nws
    gen_names = Vector{String}()  # for the labels of the graph
    phases = 3 # number of multiconductor phases

     # Create arrays to store values (depending on total type)
     if totals == :real
        active_gen = Array{Real}(undef, length(x), gens_number*phases)
    elseif totals == :imaginary
        reactive_gen = Array{Real}(undef, length(x), gens_number*phases)
    else
        active_gen = Array{Real}(undef, length(x), gens_number*phases)
        reactive_gen = Array{Real}(undef, length(x), gens_number*phases)
    end

    for (n, nw) in get(mn_solution, "nw", Dict())
        gen_counter = 0 # counts the number of gen. devices (gens, solar, storage).

        for (type, settings) in sources
            real_key = get(settings, "inactive_real", "" => 0)[1]
            imag_key = get(settings, "inactive_imaginary", "" => 0)[1]

            for (name,obj) in get(nw, type, Dict())
                gen_counter += 1

                if totals == :real
                    v_real = get(obj, real_key, 0.0)
                    v_imag = 0.0
                elseif totals == :imaginary
                    v_real = 0.0
                    v_imag = get(obj, imag_key, 0.0)
                else
                    v_real = get(obj, real_key, 0.0)
                    v_imag = get(obj, imag_key, 0.0)
                end

                # Loop through phases, save names in vector and values in arrays
                for p in phases:-1:1
                    # only in the first nw (save names of the gens)
                    if (parse(Int, n) == 1)
                        push!(gen_names, name*"_$p")
                    end

                    # check that real or imag values are not zero
                    if totals == :real
                        active_gen[parse(Int, n),(gen_counter*phases)-p+1] = v_real[p]
                    elseif totals == :imaginary
                        reactive_gen[parse(Int, n),(gen_counter*phases)-p+1] = v_imag[p]
                    else
                        active_gen[parse(Int, n),(gen_counter*phases)-p+1] = v_real[p]
                        reactive_gen[parse(Int, n),(gen_counter*phases)-p+1] = v_imag[p]
                    end
                end
            end
        end
    end

    # Assign active/reactive power for plotting
    if totals == :real
        gen_power = active_gen
    elseif totals == :imaginary
        gen_power = reactive_gen
    else
        error("Can't plot active and reactive powers at the same time for this function.")
    end

    power_scale_factor = mn_solution["nw"]["1"]["settings"]["sbase"]

    if totals == :real
        units_str = power_scale_factor == 1e2 ? "W" : power_scale_factor == 1e5 ? "kW" : power_scale_factor == 1e10 ? "MW" : "$power_scale_factor W"
    elseif totals == :imaginary
        units_str = power_scale_factor == 1e2 ? "VAR" : power_scale_factor == 1e5 ? "kVAR" : power_scale_factor == 1e10 ? "MVAR" : "$power_scale_factor VAR"
    else
        units_str = power_scale_factor == 1e2 ? "VA" : power_scale_factor == 1e5 ? "kVA" : power_scale_factor == 1e10 ? "MVA" : "$power_scale_factor VA"
    end

    spec = deepcopy(default_source_demand_summary_spec)

    power_summary_data = []
    for i in 1:length(x)
        for (c, y) in enumerate([gen_power[i,:]])
            for (v,b) in enumerate(y)
                push!(power_summary_data, Dict("x" => x[i], "y" => b, "c" => v))
            end
        end
    end

    @set! spec.data = [
        Dict(
            "name" => "table",
            "values" => power_summary_data,
            "transform" => [
                Dict(
                    "type" => "stack",
                    "groupby" => ["x"],
                    "sort" => Dict("field" => "c"),
                    "field" => "y"
                )
            ]
        )
    ]

    @set! spec.scales[4]["range"] = gen_names # set legends
    @set! spec.axes[2]["title"] = "Power ($units_str)"

    Vega.save(fileout, spec)

    return spec
end



"""
    `function plot_transmission_source_demand_summary(
        file::String,
        mn_solution::Dict{String,<:Any};
        kwargs...
    )`

    Plots the individual generation from sources of the transmission system (load demands are not available from PM).

    Arguments:

    `fileout::String`: path to file where plot will be saved
    `mn_solution::Dict{String,<:Any}`: a multinetwork solution from PM or PMITD.
    `yscale::Symbol`: To set log scale, `:log10`, else `:identity`
    `legend_position::Symbol`: Position of legend, accepts the following symbols: `:right`, `:left`, `:top`, `:bottom`, `:inside`, `:best`, `:legend`, `:topright`, `:topleft`, `:bottomleft`, `:bottomright`
    `sources::Dict{String,<:Any}`: information about sources (e.g. generators)
    `totals::Symbol`: Choose `:real` or `:imaginary` (`:complex` not available for this function)
"""
function plot_transmission_source_demand_summary(fileout::String, mn_solution::Dict{String,<:Any};
    yscale::Symbol=:identity,
    legend_position::Symbol=:best,
    sources::Dict{String,<:Any}=default_sources_math,
    totals::Symbol=:real,
    )::Vega.VGSpec

    # assign correct data based on multiinfrastructure condition
    if get(mn_solution, "multiinfrastructure", false) != false
        # Extract pm solution only
        mn_solution = mn_solution["it"]["pm"]
    end

    x = 1:length(mn_solution["nw"]) # number of nws
    gen_names = Vector{String}()  # for the labels of the graph

    # Get total number of gens, solar sotrage
    gens_number = length(get(mn_solution["nw"]["1"], "gen", Dict()))

    # Create arrays to store values (depending on total type)
    if totals == :real
        active_gen = Array{Real}(undef, length(x), gens_number)
    elseif totals == :imaginary
        reactive_gen = Array{Real}(undef, length(x), gens_number)
    else
        active_gen = Array{Real}(undef, length(x), gens_number)
        reactive_gen = Array{Real}(undef, length(x), gens_number)
    end

    for (n, nw) in get(mn_solution, "nw", Dict())
        gen_counter = 0 # counts the number of gen. devices (gens, solar, storage).

        for (type, settings) in sources
            real_key = get(settings, "inactive_real", "" => 0)[1]
            imag_key = get(settings, "inactive_imaginary", "" => 0)[1]

            for (name,obj) in get(nw, type, Dict())
                gen_counter += 1

                if totals == :real
                    v_real = get(obj, real_key, 0.0)
                    v_imag = 0.0
                    active_gen[parse(Int, n),gen_counter] = sum(v_real)
                elseif totals == :imaginary
                    v_real = 0.0
                    v_imag = get(obj, imag_key, 0.0)
                    reactive_gen[parse(Int, n),gen_counter] = sum(v_imag)
                else
                    v_real = get(obj, real_key, 0.0)
                    v_imag = get(obj, imag_key, 0.0)
                    active_gen[parse(Int, n),gen_counter] = sum(v_real)
                    reactive_gen[parse(Int, n),gen_counter] = sum(v_imag)
                end

                # only in the first nw (save names of the gens)
                if (parse(Int, n) == 1)
                    push!(gen_names, name)
                end

            end
        end
    end

    # Assign active/reactive power for plotting
    if totals == :real
        gen_power = active_gen
    elseif totals == :imaginary
        gen_power = reactive_gen
    else
        error("Can't plot active and reactive powers at the same time for this funciton.")
    end

    power_scale_factor = mn_solution["nw"]["1"]["baseMVA"]

    if totals == :real
        units_str = power_scale_factor == 1e-6 ? "W" : power_scale_factor == 0.001 ? "kW" : power_scale_factor == 100 ? "MW" : "$power_scale_factor W"
    elseif totals == :imaginary
        units_str = power_scale_factor == 1e-6 ? "VAR" : power_scale_factor == 0.001 ? "kVAR" : power_scale_factor == 100 ? "MVAR" : "$power_scale_factor VAR"
    else
        units_str = power_scale_factor == 1e-6 ? "VA" : power_scale_factor == 0.001 ? "kVA" : power_scale_factor == 100 ? "MVA" : "$power_scale_factor VA"
    end

    spec = deepcopy(default_source_demand_summary_spec)

    power_summary_data = []
    counter = 0
    for i in 1:length(x)
        elem_count = 0
        for (c, y) in enumerate([gen_power[i,:]])
            for (v,b) in enumerate(y)
                elem_count += 1
                push!(power_summary_data, Dict("x" => x[i], "y" => b, "c" => elem_count-1))
            end
        end
    end

    @set! spec.data = [
        Dict(
            "name" => "table",
            "values" => power_summary_data,
            "transform" => [
                Dict(
                    "type" => "stack",
                    "groupby" => ["x"],
                    "sort" => Dict("field" => "c"),
                    "field" => "y"
                )
            ]
        )
    ]

    @set! spec.scales[4]["range"] = gen_names # set legends
    @set! spec.axes[2]["title"] = "Power ($units_str)"      # set title

    Vega.save(fileout, spec)

    return spec
end



"""
    `function plot_boundary_summary(
        file::String,
        mn_solution::Dict{String,<:Any};
        kwargs...
    )`

    Plots the power flowing at the boundary(ies) of the T&D system. to (distribution) and from (transmission).

    Arguments:

    `fileout::String`: path to file where plot will be saved
    `mn_solution::Dict{String,<:Any}`: a multinetwork solution from PMITD.
    `yscale::Symbol`: To set log scale, `:log10`, else `:identity`
    `legend_position::Symbol`: Position of legend, accepts the following symbols: `:right`, `:left`, `:top`, `:bottom`, `:inside`, `:best`, `:legend`, `:topright`, `:topleft`, `:bottomleft`, `:bottomright`
    `totals::Symbol`: Choose `:real` or `:imaginary` (`:complex` not available for this function)
"""
function plot_boundary_summary(fileout::String, mn_solution::Dict{String,<:Any};
    yscale::Symbol=:identity,
    legend_position::Symbol=:best,
    totals::Symbol=:real,
    )::Vega.VGSpec

    BOUNDARY_NUMBER = 100000

    # Needed just to get the basis of the units
    mn_pmd = mn_solution["it"]["pmd"]

    # Extract pmitd solution only
    mn_solution = mn_solution["it"]["pmitd"]

    x = 1:length(mn_solution["nw"]) # number of nws
    phases = 3

    # Get total number of gens, solar sotrage
    boundaries_number = length(mn_solution["nw"]["1"]["boundary"])/2 # divide by two
    boundaries_number = Int(boundaries_number)
    boundary_names = Vector{String}(undef, (boundaries_number*phases)+boundaries_number)  # for the labels of the graph

    # Create arrays to store values (depending on total type)
    if totals == :real
        active_bound = Array{Real}(undef, length(x), (boundaries_number*phases)+boundaries_number)
    elseif totals == :imaginary
        reactive_bound = Array{Real}(undef, length(x), (boundaries_number*phases)+boundaries_number)
    else
        active_bound = Array{Real}(undef, length(x), (boundaries_number*phases)+boundaries_number)
        reactive_bound = Array{Real}(undef, length(x), (boundaries_number*phases)+boundaries_number)
    end

    for (n, nw) in get(mn_solution, "nw", Dict())
        for itd in nw["boundary"]
            boundary_counter = eval(Meta.parse(itd[1])) # parse tuple from string
            boundary_counter = boundary_counter[1]-BOUNDARY_NUMBER # get boundary number - BIAS
            # distribution system boundaries
            for boundary in itd[2]
                if totals == :real
                    if boundary[1]=="pbound_to"
                        for p in phases:-1:1
                            if (parse(Int, n) == 1)
                                boundary_names[(boundary_counter*phases)-p+1] = string(itd[1])*"_$p"
                            end
                            active_bound[parse(Int, n), (boundary_counter*phases)-p+1] = boundary[2][p]
                        end
                    elseif boundary[1]=="pbound_fr"
                        # only add 1 if it is a multi-system problem
                        if (parse(Int, n) == 1)
                            boundary_names[(boundaries_number*phases)+boundary_counter] = string(itd[1])
                        end
                        active_bound[parse(Int, n), (boundaries_number*phases)+boundary_counter] = boundary[2][1]
                    end
                elseif totals == :imaginary
                    if boundary[1]=="qbound_to"
                        for p in phases:-1:1
                            if (parse(Int, n) == 1)
                                boundary_names[(boundary_counter*phases)-p+1] = string(itd[1])*"_$p"
                            end
                            reactive_bound[parse(Int, n), (boundary_counter*phases)-p+1] = boundary[2][p]
                        end
                    elseif boundary[1]=="qbound_fr"
                        # only add 1 if it is a multi-system problem
                        if (parse(Int, n) == 1)
                            boundary_names[(boundaries_number*phases)+boundary_counter] = string(itd[1])
                        end
                        reactive_bound[parse(Int, n), (boundaries_number*phases)+boundary_counter] = boundary[2][1]
                    end
                else
                    error("Can't plot active and reactive powers at the same time.")
                end
            end
        end
    end

    # Assign active/reactive power for plotting
    if totals == :real
        boundary_power = active_bound
    elseif totals == :imaginary
        boundary_power = reactive_bound
    else
        error("Can't plot active and reactive powers at the same time for this function.")
    end

    power_scale_factor = mn_pmd["nw"]["1"]["settings"]["sbase"]

    if totals == :real
        units_str = power_scale_factor == 1e2 ? "W" : power_scale_factor == 1e5 ? "kW" : power_scale_factor == 1e10 ? "MW" : "$power_scale_factor W"
    elseif totals == :imaginary
        units_str = power_scale_factor == 1e2 ? "VAR" : power_scale_factor == 1e5 ? "kVAR" : power_scale_factor == 1e10 ? "MVAR" : "$power_scale_factor VAR"
    else
        units_str = power_scale_factor == 1e2 ? "VA" : power_scale_factor == 1e5 ? "kVA" : power_scale_factor == 1e10 ? "MVA" : "$power_scale_factor VA"
    end

    spec = deepcopy(default_source_demand_summary_spec)

    power_summary_data = []
    for i in 1:length(x)
        for (c, y) in enumerate([boundary_power[i,:]])
            for (v,b) in enumerate(y)
                push!(power_summary_data, Dict("x" => x[i], "y" => b, "c" => v))
            end
        end
    end

    @set! spec.data = [
        Dict(
            "name" => "table",
            "values" => power_summary_data,
            "transform" => [
                Dict(
                    "type" => "stack",
                    "groupby" => ["x"],
                    "sort" => Dict("field" => "c"),
                    "field" => "y"
                )
            ]
        )
    ]

    @set! spec.scales[4]["range"] = boundary_names # set legends
    @set! spec.axes[2]["title"] = "Power ($units_str)"      # set title

    Vega.save(fileout, spec)

    return spec

end


"""
    `function plot_transmission_voltage_magnitudes(
        file::String,
        mn_solution::Dict{String,<:Any};
        kwargs...
    )`

    Plots the voltage magnitude of buses in the transmission system.

    Arguments:

    `fileout::String`: path to file where plot will be saved
    `mn_solution::Dict{String,<:Any}`: a multinetwork solution from PM or PMITD.
    `nw::Int=0`: multinetwork network (step) to be plotted. If not given, mn_solution must come from a non-multinetwork solution.
"""
function plot_transmission_voltage_magnitudes(fileout::String, mn_solution::Dict{String,<:Any};
    nw::Int=0,
    )::Vega.VGSpec

    # convert nw to str
    nw_str = string(nw)

    # assign correct data based on multiinfrastructure condition
    if get(mn_solution, "multiinfrastructure", false) != false
        # Get the respective nw to plot
        if nw==0
            # Extract pm solution only
            mn_pm = mn_solution["it"]["pm"]
        else
            # Extract pm solution only
            mn_pm = mn_solution["it"]["pm"]["nw"][nw_str]
        end
    else
        # Get the respective nw to plot
        if nw!=0
            mn_pm = mn_solution["nw"][nw_str]
        else
            mn_pm = mn_solution
        end
    end

    # Check for a possible usual error.
    # Error: User passing a multinetwork solution and not providing the nw.
    if nw==0
        try
            temp_assig = mn_pm["bus"]
        catch e
            error("User not providing the 'nw' (step) for the multinetwork solution. Please provide one.")
        end
    end

    # ---- Transmission plot ----
    tbus_names = Vector{String}()
    tbus_voltages = Vector{Float64}()

    for bus in mn_pm["bus"]
        push!(tbus_names, string(bus[1]))
        push!(tbus_voltages, bus[2]["vm"])
    end

    spec = deepcopy(default_voltage_magnitude_spec)

    trans_voltage_data = []
    for (i,vmag) in enumerate(tbus_voltages)
        push!(trans_voltage_data, Dict("voltage" => vmag, "busname" => tbus_names[i]))
    end

    @set! spec.data = [
        Dict(
            "name" => "source",
            "values" => trans_voltage_data
        )
    ]

    Vega.save(fileout, spec)

    return spec

end


"""
    `function plot_distribution_voltage_magnitudes(
        file::String,
        mn_solution::Dict{String,<:Any};
        kwargs...
    )`

    Plots the per-phase voltage magnitude of buses in the distribution system.

    Arguments:

    `fileout::String`: path to file where plot will be saved
    `mn_solution::Dict{String,<:Any}`: a multinetwork solution from PMD or PMITD.
    `nw::Int=0`: multinetwork network (step) to be plotted. If not given, mn_solution must come from a non-multinetwork solution.
"""
function plot_distribution_voltage_magnitudes(fileout::String, mn_solution::Dict{String,<:Any};
    nw::Int=0,
    )::Vega.VGSpec

    # convert nw to str
    nw_str = string(nw)

    # assign correct data based on multiinfrastructure condition
    if get(mn_solution, "multiinfrastructure", false) != false
        # Get the respective nw to plot
        if nw==0
            # Extract pmd solution only
            mn_pmd = mn_solution["it"]["pmd"]
        else
            # Extract pmd solution only
            mn_pmd = mn_solution["it"]["pmd"]["nw"][nw_str]
        end
    else
        # Get the respective nw to plot
        if nw!=0
            mn_pmd = mn_solution["nw"][nw_str]
        else
            mn_pmd = mn_solution
        end
    end

    # Check for a possible usual error.
    # Error: User passing a multinetwork solution and not providing the nw.
    if nw==0
        try
            temp_assig = mn_pmd["bus"]
        catch e
            error("User not providing the 'nw' (step) for the multinetwork solution. Please provide one.")
        end
    end


    # ---- Distribution plot ----
    dbus_names = Vector{String}()
    dbus_voltages = Vector{Float64}()

    for bus in mn_pmd["bus"]

        # get the number of phases of the bus
        phases = length(bus[2]["vm"])

        if phases == 1
            push!(dbus_names, string(bus[1])*"_1")
            push!(dbus_voltages, bus[2]["vm"][1])
        elseif phases == 2
            for j in 1:1:2
                push!(dbus_names, string(bus[1])*"_$j")
                push!(dbus_voltages, bus[2]["vm"][j])
            end
        elseif phases == 3
            for j in 1:1:3
                push!(dbus_names, string(bus[1])*"_$j")
                push!(dbus_voltages, bus[2]["vm"][j])
            end
        else
            @warn "There are more than 3 phases! Plot may be wrong."
        end
    end

    spec = deepcopy(default_voltage_magnitude_spec)

    dist_voltage_data = []
    for (i,vmag) in enumerate(dbus_voltages)
        push!(dist_voltage_data, Dict("voltage" => vmag, "busname" => dbus_names[i]))
    end

    @set! spec.data = [
        Dict(
            "name" => "source",
            "values" => dist_voltage_data
        )
    ]

    Vega.save(fileout, spec)

    return spec

end


"""
    `function plot_transmission_duals(
        file::String,
        mn_solution::Dict{String,<:Any};
        kwargs...
    )`

    Plots the duals of buses in the transmission system.

    Arguments:

    `fileout::String`: path to file where plot will be saved
    `mn_solution::Dict{String,<:Any}`: a multinetwork solution from PM or PMITD.
    `nw::Int=0`: multinetwork network (step) to be plotted. If not given, mn_solution must come from a non-multinetwork solution.
"""
function plot_transmission_duals(fileout::String, mn_solution::Dict{String,<:Any};
    nw::Int=0,
    )::Vega.VGSpec

    # convert nw to str
    nw_str = string(nw)

    # assign correct data based on multiinfrastructure condition
    if get(mn_solution, "multiinfrastructure", false) != false
        # Get the respective nw to plot
        if nw==0
            # Extract pm solution only
            mn_pm = mn_solution["it"]["pm"]
        else
            # Extract pm solution only
            mn_pm = mn_solution["it"]["pm"]["nw"][nw_str]
        end
    else
        # Get the respective nw to plot
        if nw!=0
            mn_pm = mn_solution["nw"][nw_str]
        else
            mn_pm = mn_solution
        end
    end

    # Check for a possible usual error.
    # Error: User passing a multinetwork solution and not providing the nw.
    if nw==0
        try
            temp_assig = mn_pm["bus"]
        catch e
            error("User not providing the 'nw' (step) for the multinetwork solution. Please provide one.")
        end
    end

    # ---- Transmission plot ----
    tbus_names = Vector{String}()
    tbus_duals = Vector{Float64}()

    for bus in mn_pm["bus"]
        # TODO: probably a horrible way to avoid a single entry that does not exists (slack bus).
        try
            push!(tbus_duals, bus[2]["lam_kcl_r"])
            push!(tbus_names, string(bus[1]))
        catch e
            continue
        end
    end

    spec = deepcopy(default_duals_spec)

    trans_duals_data = []
    for (i,dual) in enumerate(tbus_duals)
        push!(trans_duals_data, Dict("duals" => dual, "busname" => tbus_names[i]))
    end

    @set! spec.data = [
        Dict(
            "name" => "source",
            "values" => trans_duals_data
        )
    ]

    Vega.save(fileout, spec)

    return spec

end


"""
    `function plot_histogram_transmission_duals(
        file::String,
        mn_solution::Dict{String,<:Any};
        kwargs...
    )`

    Plots histogram counting the duals of buses in the transmission system.

    Arguments:

    `fileout::String`: path to file where plot will be saved
    `mn_solution::Dict{String,<:Any}`: a multinetwork solution from PM or PMITD.
"""
function plot_histogram_transmission_duals(fileout::String, mn_solution::Dict{String,<:Any}
    )::Vega.VGSpec

    # assign correct data based on multiinfrastructure condition
    if get(mn_solution, "multiinfrastructure", false) != false
        mn_pm = mn_solution["it"]["pm"]
    else
        mn_pm = mn_solution
    end


    # # ---- Transmission plot ----
    tbus_duals = Vector{Float64}()

    # check if the problem is multinetwork
    if get(mn_pm, "multinetwork", false) != false
        # run for all nws
        for (nw, nw_data) in mn_pm["nw"]
            for bus in nw_data["bus"]
                try
                    push!(tbus_duals, bus[2]["lam_kcl_r"])
                catch e
                    continue
                end
            end
        end

    else

        for bus in mn_pm["bus"]
            try
                push!(tbus_duals, bus[2]["lam_kcl_r"])
            catch e
                continue
            end
        end

    end

    # get the spec
    spec = deepcopy(default_duals_histogram_spec)

    # add to dictionary using specific format
    trans_duals_data = []
    for (i,dual) in enumerate(tbus_duals)
        push!(trans_duals_data, Dict("duals" => dual))
    end

    @set! spec.data = [
        Dict(
            "name" => "table",
            "values" => trans_duals_data,
            "transform" => [
            Dict(
                "type" => "extent",
                "field" => "duals",
                "signal" => "extent"
            ),
            Dict(
                "type" => "bin",
                "signal" => "bins",
                "field" => "duals",
                "extent" => Dict("signal" => "extent"),
                "maxbins" => Dict("signal" => "maxbins")
            )
        ]),
        Dict(
            "name" => "counts",
            "source" => "table",
            "transform" => [
            Dict(
                "type" => "filter",
                "expr" => "datum['duals'] != null"
            ),
            Dict(
                "type" => "aggregate",
                "groupby" => ["bin0", "bin1"]
            )
        ])
    ]

    Vega.save(fileout, spec)

    return spec

end


"""
    `function plot_histogram_distribution_duals(
        file::String,
        mn_solution::Dict{String,<:Any};
        kwargs...
    )`

    Plots histogram counting the duals of buses in the distribution system.

    Arguments:

    `fileout::String`: path to file where plot will be saved
    `mn_solution::Dict{String,<:Any}`: a multinetwork solution from PM or PMITD.
"""
function plot_histogram_distribution_duals(fileout::String, mn_solution::Dict{String,<:Any}
    )::Vega.VGSpec

    # assign correct data based on multiinfrastructure condition
    if get(mn_solution, "multiinfrastructure", false) != false
        mn_pmd = mn_solution["it"]["pmd"]
    else
        mn_pmd = mn_solution
    end

    # # ---- Transmission plot ----
    dbus_duals = Vector{Float64}()

    # check if the problem is multinetwork
    if haskey(mn_pmd, "nw")
        # run for all nws
        for (nw, nw_data) in mn_pmd["nw"]

            # scale LMPs - TODO: Temporary, since this should happen in PMD and PMITD
            if get(nw_data, "per_unit", false) != true
                power_scale_factor = nw_data["settings"]["sbase"]
            else
                power_scale_factor = 1
            end

            for bus in nw_data["bus"]
                # get the number of phases of the bus
                phases = length(bus[2]["lam_kcl_r"])

                # loop through phases
                for phase in 1:1:phases
                    push!(dbus_duals, bus[2]["lam_kcl_r"][phase]/power_scale_factor)
                end
            end
        end
    else
        # scale LMPs - TODO: Temporary, since this should happen in PMD and PMITD
        if get(mn_pmd, "per_unit", false) != true
            power_scale_factor = mn_pmd["settings"]["sbase"]
        else
            power_scale_factor = 1
        end

        for bus in mn_pmd["bus"]
            # get the number of phases of the bus
            phases = length(bus[2]["lam_kcl_r"])

            # loop through phases
            for phase in 1:1:phases
                push!(dbus_duals, bus[2]["lam_kcl_r"][phase]/power_scale_factor)
            end
        end
    end

    # get the spec
    spec = deepcopy(default_duals_histogram_spec)

    # add to dictionary using specific format
    dist_duals_data = []
    for (i,dual) in enumerate(dbus_duals)
        push!(dist_duals_data, Dict("duals" => dual))
    end

    @set! spec.data = [
        Dict(
            "name" => "table",
            "values" => dist_duals_data,
            "transform" => [
            Dict(
                "type" => "extent",
                "field" => "duals",
                "signal" => "extent"
            ),
            Dict(
                "type" => "bin",
                "signal" => "bins",
                "field" => "duals",
                "extent" => Dict("signal" => "extent"),
                "maxbins" => Dict("signal" => "maxbins")
            )
        ]),
        Dict(
            "name" => "counts",
            "source" => "table",
            "transform" => [
            Dict(
                "type" => "filter",
                "expr" => "datum['duals'] != null"
            ),
            Dict(
                "type" => "aggregate",
                "groupby" => ["bin0", "bin1"]
            )
        ])
    ]

    Vega.save(fileout, spec)

    return spec

end


"""
    `function plot_histogram_transmission_voltage_magnitudes(
        file::String,
        mn_solution::Dict{String,<:Any};
        kwargs...
    )`

    Plots histogram counting the voltage magnitudes of buses in the transmission system.

    Arguments:

    `fileout::String`: path to file where plot will be saved
    `mn_solution::Dict{String,<:Any}`: a multinetwork solution from PM or PMITD.
"""
function plot_histogram_transmission_voltage_magnitudes(fileout::String, mn_solution::Dict{String,<:Any}
    )::Vega.VGSpec

    # assign correct data based on multiinfrastructure condition
    if get(mn_solution, "multiinfrastructure", false) != false
        mn_pm = mn_solution["it"]["pm"]
    else
        mn_pm = mn_solution
    end


    # # ---- Transmission plot ----
    tbus_vmags = Vector{Float64}()

    # check if the problem is multinetwork
    if get(mn_pm, "multinetwork", false) != false
        # run for all nws
        for (nw, nw_data) in mn_pm["nw"]
            for bus in nw_data["bus"]
                try
                    push!(tbus_vmags, bus[2]["vm"])
                catch e
                    continue
                end
            end
        end

    else

        for bus in mn_pm["bus"]
            try
                push!(tbus_vmags, bus[2]["vm"])
            catch e
                continue
            end
        end

    end

    # get the spec
    spec = deepcopy(default_voltage_histogram_spec)

    # add to dictionary using specific format
    trans_vmags_data = []
    for (i,vmag) in enumerate(tbus_vmags)
        push!(trans_vmags_data, Dict("vmags" => vmag))
    end

    @set! spec.data = [
        Dict(
            "name" => "table",
            "values" => trans_vmags_data,
            "transform" => [
            Dict(
                "type" => "extent",
                "field" => "vmags",
                "signal" => "extent"
            ),
            Dict(
                "type" => "bin",
                "signal" => "bins",
                "field" => "vmags",
                "extent" => Dict("signal" => "extent"),
                "maxbins" => Dict("signal" => "maxbins")
            )
        ]),
        Dict(
            "name" => "counts",
            "source" => "table",
            "transform" => [
            Dict(
                "type" => "filter",
                "expr" => "datum['vmags'] != null"
            ),
            Dict(
                "type" => "aggregate",
                "groupby" => ["bin0", "bin1"]
            )
        ])
    ]

    Vega.save(fileout, spec)

    return spec

end


"""
    `function plot_histogram_distribution_voltage_magnitudes(
        file::String,
        mn_solution::Dict{String,<:Any};
        kwargs...
    )`

    Plots histogram counting the voltage magnitudes of buses in the distribution system.

    Arguments:

    `fileout::String`: path to file where plot will be saved
    `mn_solution::Dict{String,<:Any}`: a multinetwork solution from PM or PMITD.
"""
function plot_histogram_distribution_voltage_magnitudes(fileout::String, mn_solution::Dict{String,<:Any}
    )::Vega.VGSpec

    # assign correct data based on multiinfrastructure condition
    if get(mn_solution, "multiinfrastructure", false) != false
        mn_pmd = mn_solution["it"]["pmd"]
    else
        mn_pmd = mn_solution
    end

    # # ---- Transmission plot ----
    dbus_vmags = Vector{Float64}()

    # check if the problem is multinetwork
    if haskey(mn_pmd, "nw")
        # run for all nws
        for (nw, nw_data) in mn_pmd["nw"]
            for bus in nw_data["bus"]
                # get the number of phases of the bus
                phases = length(bus[2]["vm"])

                # loop through phases
                for phase in 1:1:phases
                    push!(dbus_vmags, bus[2]["vm"][phase])
                end
            end
        end
    else
        for bus in mn_pmd["bus"]
            # get the number of phases of the bus
            phases = length(bus[2]["vm"])

            # loop through phases
            for phase in 1:1:phases
                push!(dbus_vmags, bus[2]["vm"][phase])
            end
        end
    end

    # get the spec
    spec = deepcopy(default_voltage_histogram_spec)

    # add to dictionary using specific format
    dist_vmags_data = []
    for (i,vmag) in enumerate(dbus_vmags)
        push!(dist_vmags_data, Dict("vmags" => vmag))
    end

    @set! spec.data = [
        Dict(
            "name" => "table",
            "values" => dist_vmags_data,
            "transform" => [
            Dict(
                "type" => "extent",
                "field" => "vmags",
                "signal" => "extent"
            ),
            Dict(
                "type" => "bin",
                "signal" => "bins",
                "field" => "vmags",
                "extent" => Dict("signal" => "extent"),
                "maxbins" => Dict("signal" => "maxbins")
            )
        ]),
        Dict(
            "name" => "counts",
            "source" => "table",
            "transform" => [
            Dict(
                "type" => "filter",
                "expr" => "datum['vmags'] != null"
            ),
            Dict(
                "type" => "aggregate",
                "groupby" => ["bin0", "bin1"]
            )
        ])
    ]

    Vega.save(fileout, spec)

    return spec

end


"""
    `function plot_distribution_duals(
        file::String,
        mn_solution::Dict{String,<:Any};
        kwargs...
    )`

    Plots the per_phase duals of buses in the distribution system.

    Arguments:

    `fileout::String`: path to file where plot will be saved
    `mn_solution::Dict{String,<:Any}`: a multinetwork solution from PMD or PMITD.
    `nw::Int=0`: multinetwork network (step) to be plotted. If not given, mn_solution must come from a non-multinetwork solution.
"""
function plot_distribution_duals(fileout::String, mn_solution::Dict{String,<:Any};
    nw::Int=0,
    )::Vega.VGSpec

    # convert nw to str
    nw_str = string(nw)

    # assign correct data based on multiinfrastructure condition
    if get(mn_solution, "multiinfrastructure", false) != false
        # Get the respective nw to plot
        if nw==0
            # Extract pmd solution only
            mn_pmd = mn_solution["it"]["pmd"]
        else
            # Extract pmd solution only
            mn_pmd = mn_solution["it"]["pmd"]["nw"][nw_str]
        end
    else
        # Get the respective nw to plot
        if nw!=0
            mn_pmd = mn_solution["nw"][nw_str]
        else
            mn_pmd = mn_solution
        end
    end

    # Check for a possible usual error.
    # Error: User passing a multinetwork solution and not providing the nw.
    if nw==0
        try
            temp_assig = mn_pmd["bus"]
        catch e
            error("User not providing the 'nw' (step) for the multinetwork solution. Please provide one.")
        end
    end

    # ---- Distribution plot ----
    dbus_names = Vector{String}()
    dbus_duals = Vector{Float64}()

    for bus in mn_pmd["bus"]

        # get the number of phases of the bus
        phases = length(bus[2]["lam_kcl_r"])

        if phases == 1
            push!(dbus_names, string(bus[1])*"_1")
            push!(dbus_duals, bus[2]["lam_kcl_r"][1])
        elseif phases == 2
            for j in 1:1:2
                push!(dbus_names, string(bus[1])*"_$j")
                push!(dbus_duals, bus[2]["lam_kcl_r"][j])
            end
        elseif phases == 3
            for j in 1:1:3
                push!(dbus_names, string(bus[1])*"_$j")
                push!(dbus_duals, bus[2]["lam_kcl_r"][j])
            end
        else
            @warn "There are more than 3 phases! Plot may be wrong."
        end
    end

    spec = deepcopy(default_duals_spec)

    dist_duals_data = []
    for (i,dual) in enumerate(dbus_duals)
        push!(dist_duals_data, Dict("duals" => dual, "busname" => dbus_names[i]))
    end

    @set! spec.data = [
        Dict(
            "name" => "source",
            "values" => dist_duals_data
        )
    ]

    Vega.save(fileout, spec)

    return spec

end

