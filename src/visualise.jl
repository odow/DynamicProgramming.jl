#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
const ASSET_DIR = joinpath(dirname(dirname(@__FILE__)), "assets")
const HTML_FILE = joinpath(ASSET_DIR, "visualise.html")
const ASSETS = ["d3.v3.min.js", "visualise.js", "visualise.css"]
const PLOT_DATA = Dict{String,Any}(
    "cumulative" => false,
    "title" => "",
    "ylabel" => "",
    "xlabel" => "Stages",
    "interpolate" => "linear",
    "ymin" => "",
    "ymax" => "",
)

function add_to_output!(output::Dict{String,Any}, sym, value)
    if haskey(output, string(sym))
        output[string(sym)] = value
    else
        error("Keyword $(sym)=$value not recognised in @visualise.")
    end
end

function adddata!(
    output::Dict{String,Any},
    plot_data::Function,
    results::Dict{Symbol,Any},
)
    output["data"] = Vector{Float64}[]
    for i in 1:results[:n_replications]
        push!(output["data"], Float64[])
        y = 0.0
        for j in 1:results[:n_stages]
            if output["cumulative"]
                y += plot_data(j, i, results)
            else
                y = plot_data(j, i, results)
            end
            push!(output["data"][i], y)
        end
    end
end

function launch_plot(html_file)
    if is_windows()
        run(`$(ENV["COMSPEC"]) /c start $(html_file)`)
    elseif is_apple()
        run(`open $(html_file)`)
    elseif is_linux() || is_bsd()
        run(`xdg-open $(html_file)`)
    end
end

macro visualise(results, stage, replication, block)
    @assert block.head == :block || error("Invalid syntax for @visualise")
    kw = Expr(:tuple, stage, replication, results)
    code = quote
        plot_list = Dict{String,Any}[]
    end
    for it in block.args
        Base.Meta.isexpr(it, :line) && continue
        output = copy(PLOT_DATA)
        if it.head == :tuple
            if length(it.args) == 2
                if it.args[2].head == :tuple
                    for arg in it.args[2].args
                        if arg.head != :(=)
                            error(
                                "Must be a keyword argument in @visualise: $(arg)",
                            )
                        end
                        add_to_output!(output, arg.args[1], arg.args[2])
                    end
                elseif it.args[2].head == :(=)
                    add_to_output!(
                        output,
                        it.args[2].args[1],
                        it.args[2].args[2],
                    )
                end
                f = Expr(:->, kw, Expr(:block, it.args[1]))
            elseif length(it.args) == 1
                f = Expr(:->, kw, Expr(:block, it.args))
            else
                error("Unknown arguments in @visualise")
            end
        else
            f = Expr(:->, kw, Expr(:block, it))
        end
        push!(code.args, quote
            adddata!($output, $(esc(f)), $(esc(results)))
            push!(plot_list, $output)
        end)
    end
    push!(
        code.args,
        quote
            temporary_html_file = replace(tempname(), ".tmp", ".html")
            html_string = readstring($HTML_FILE)
            for asset in $ASSETS
                cp(
                    joinpath(ASSET_DIR, asset),
                    joinpath(dirname(temporary_html_file), asset),
                    remove_destination = true,
                )
            end
            json_data = json(plot_list)
            open(temporary_html_file, "w") do f
                return write(f, replace(html_string, "<!--DATA-->", json_data))
            end
            launch_plot(temporary_html_file)
        end,
    )
    return code
end
