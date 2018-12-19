
function plot_branch_impedance(data::Dict{String,Any})
    r = [branch["br_r"] for (i,branch) in data["branch"]]
    x = [branch["br_x"] for (i,branch) in data["branch"]]

    s = scatter(r, x, xlabel="resistance (p.u.)", ylabel="reactance (p.u.)", label="")
    r_h = histogram(r, xlabel="resistance (p.u.)", ylabel="branch count", label="", reuse=false)
    x_h = histogram(x, xlabel="reactance (p.u.)", ylabel="branch count", label="", reuse=false)
end
