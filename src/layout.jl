indexmap(xs) = Dict(x=>i for (i,x) in enumerate(xs))

"""
    treepositions(n::Node)

Get the coordinates for plotting a tree using Plots.jl. Only for bifurcating
trees! Based on the dendrogram implementation in StatsPlots.jl.

```julia
x, y = treepositions(n)
plot(y, x, color=:black, framestyle=:none, legend=false)
```
"""
function treepositions(tr, usedistance=true)
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
    h = usedistance ? getheights(tr) : zeros(length(o))
    nodepos = Dict()
    idx = indexmap([id(n) for n in o])
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
        yn = usedistance ? nodepos[id(n)][2] : min(y1, y2) - 1
        nodepos[id(n)] = (xn, yn)
        x[:,idx[id(n)]] = [x1, x1, x2, x2]
        y[:,idx[id(n)]] = [y1, yn, yn, y2]
    end
    walk(tr)
    # if the root branch has non-zero length we add this:
    if !isnan(distance(tr))
        xr, yr = nodepos[id(tr)]
        d = distance(tr)
        x = hcat(x, [xr, xr, xr, xr])
        y = hcat(y, [yr, yr, yr-d, yr-d])
    end
    y = abs.(y .- maximum(y))
    return x, y, nodepos
end

# TODO: taxon labels for horizontal tree
@recipe function f(n::Node; orientation=:v, usedistance=true)
    nleaves = length(getleaves(n))
    legend --> false
    linecolor --> :black
    grid --> false
    yshowaxis --> false
    xshowaxis --> false
    x, y, nodepos = treepositions(n, usedistance)
    taxlabels = name.(getleaves(n))
    if orientation == :h
        yforeground_color_axis --> :white
        ygrid --> false
        ylims --> (0.5, nleaves + 0.5)
        yticks --> (1:nleaves, taxlabels)
        xlims --> (0, maximum(y))
        xticks --> false
        x, y = y, x
    else
        xforeground_color_axis --> :white
        xgrid --> false
        xlims --> (0.5, nleaves + 0.5)
        ylims --> (0, maximum(y))
        xticks --> (1:nleaves, taxlabels)
        yticks --> false
    end
    x, y
end
