"""
    treepositions(n::Node)

Get the coordinates for plotting a tree using Plots.jl. Only for bifurcating
trees! Based on the dendrogram implementatoin in StatsPlots.jl.

```julia
x, y = treepositions(n)
plot(y, x, color=:black, framestyle=:none, legend=false)
```
"""
function treepositions(tr, orientation=:horizontal)
    # we need two matrices, x and y
    # col i of x has the x coordinates of the four points of an open rectangle
    # col i of y has the corresponding y coordinates in the following order:
    # 2--3
    # |  |
    # 1  4
    # TODO: multifurcating trees require a different approach... (a real tree
    # layout algorithm probably)
    @assert isbifurcating(tr) "Only implemented for bifurcating trees..."
    o = postwalk(tr)
    n = length(o)
    h = getheights(tr)
    nodepos = Dict()
    i = 0
    for n in o
        isleaf(n) && (i += 1)
        nodepos[id(n)] = (i, h[id(n)])
    end
    x = zeros(4, n)
    y = zeros(4, n)
    function walk(n)
        isleaf(n) && return
        for c in children(n)
            walk(c)
        end
        x1, y1 = nodepos[id(n[1])]
        x2, y2 = nodepos[id(n[2])]
        xn = (x1 + x2)/2
        yn = nodepos[id(n)][2]
        nodepos[id(n)] = (xn, yn)
        x[:,id(n)] = [x1, x1, x2, x2]
        y[:,id(n)] = [y1, yn, yn, y2]
    end
    walk(tr)
    if orientation == :horizontal
        return y, x
    else
        y = abs.(y .- maximum(y))
        return x, y
    end
end

# TODO: taxon labels for horizontal tree
@recipe function f(n::Node; orientation=:vertical)
    nleaves = length(getleaves(n))
    legend --> false
    linecolor --> :black
    x, y = treepositions(n, orientation)
    taxlabels = name.(getleaves(n))
    if orientation == :horizontal
        yforeground_color_axis --> :white
        ygrid --> false
        ylims --> (0.5, nleaves + 0.5)
        yticks --> false
        yshowaxis --> false
        xlims --> (0, maximum(x))
    else
        xshowaxis --> false
        xforeground_color_axis --> :white
        xgrid --> false
        xlims --> (0.5, nleaves + 0.5)
        ylims --> (0, maximum(y))
        xticks --> (1:nleaves, taxlabels)
    end
    x, y
end
