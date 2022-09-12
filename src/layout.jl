"""
    treepositions(tree, transform=false)

Compute the coordinates for the nodes of a tree. `transform=true` will
transform the terminal branches so that the tips are all equidistant from the
root (as in FigTree).
"""
function treepositions(tr::N, transform=false) where N
    o = postwalk(tr)
    n = length(o)
    h = getheights(tr)
    maxh = maximum(values(h))
    nodepos = Dict{N,Tuple{Float64,Float64}}()
    i = 0.
    for n in o
        if isleaf(n) 
            i += 1.
            hn = transform ? maxh : h[id(n)]
        else
            hn = h[id(n)]
        end
        nodepos[n] = (hn, i)
    end
    function walk(n)
        isleaf(n) && return nodepos[n]
        cs = map(walk, children(n))
        xn = nodepos[n][1]
        yn = sum(last.(cs)) / length(cs)
        nodepos[n] = (xn, yn)
    end
    walk(tr)
    # if the root branch has non-zero length we add this:
    return nodepos
end

@recipe function f(tree::Node; sq=true, fs=9, internal=true, transform=false, scalebar=0.)
    d = treepositions(tree, transform)
    framestyle --> :none
    grid --> false
    legend --> false
    for n in prewalk(tree)
        if isroot(n) 
            if isfinite(distance(n))
                (x, y) = d[n]
                @series begin
                    seriestype := :path
                    seriescolor --> :black
                    [(x, y-distance(n)), (x,y)]
                end
            end
        else
            @series begin
                seriestype := :path
                seriescolor --> :black
                (x1, y1) = d[n]
                (x2, y2) = d[parent(n)]
                if sq
                    [(x2, y2), (x2, y1), (x1, y1)]
                else
                    [(x2, y2), (x1, y1)]
                end
            end
        end
    end
    if scalebar > 0
        @series begin
            color --> :black
            linewidth --> 1
            marker --> 3
            markershape --> :vline
            x1, x2 = extrema(first.(collect(values(d))))
            re = (x2 - x1) % scalebar
            xs = collect(x2:-scalebar:x1-re) 
            @show (x2 - x1) ÷ scalebar
            ys = repeat([0.5], length(xs))
            xs, ys
        end
    end
    ns = internal ? prewalk(tree) : getleaves(tree)
    anns = [(d[n]..., (" " * name(n), fs, :left)) for n in ns]
    @series begin
        annotations := anns
        [], []
    end
end
