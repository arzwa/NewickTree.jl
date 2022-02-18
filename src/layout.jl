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
function treepositions(tr, transform=false)
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
    maxh = maximum(values(h))
    nodepos = Dict()
    idx = indexmap(id.(o))
    i = 0
    for n in o
        if isleaf(n) 
            i += 1
            hn = transform ? maxh : h[id(n)]
        else
            hn = h[id(n)]
        end
        nodepos[id(n)] = (i, hn)
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
    return x, y, nodepos
end

# very inelegant, but works
@recipe function f(n::Node; orientation=1, pad=1., internal=false, transform=false)
    nleaves = length(getleaves(n))
    legend --> false
    linecolor --> :black
    grid --> false
    yshowaxis --> false
    xshowaxis --> false
    xticks --> false
    yticks --> false
    x, y, nodepos = treepositions(n, transform)
    if internal
        l = Dict(id(n)=>name(n) for n in postwalk(n))
    else
        l = Dict(id(n)=>name(n) for n in getleaves(n))
    end
    fontfamily --> "helvetica oblique"
    if orientation == 1
        yforeground_color_axis --> :white
        ylims --> (0.5, nleaves + 0.5)
        xlims --> (0, maximum(y) + pad)
        x, y = y, x
        anns = [(y, x, (" " * l[i], :left)) for (i,(x,y)) in nodepos if haskey(l, i)]
    elseif orientation == 2
        xforeground_color_axis --> :white
        xlims --> (0.5, nleaves + 0.5)
        ylims --> (0, maximum(y) + pad)
        anns = [(x, y, (l[i] * "\n", :bottom)) for (i,(x,y)) in nodepos if haskey(l, i)]
    elseif orientation == 3
        xforeground_color_axis --> :white
        xlims --> (-maximum(y) - pad, 0)
        ylims --> (0, nleaves + 0.5)
        x, y = -y, x
        anns = [(-y, x, (l[i] * " ", :right)) for (i,(x,y)) in nodepos if haskey(l, i)]
    elseif orientation == 4
        xforeground_color_axis --> :white
        xgrid --> false
        xlims --> (0.5, nleaves + 0.5)
        ylims --> (-maximum(y) - pad, 0)
        x, y = x, -y
        anns = [(x, -y, (l[i], :top)) for (i,(x,y)) in nodepos if haskey(l, i)]
    end
    @series begin
        annotations := anns
        x, y
    end
end
