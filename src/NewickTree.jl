module NewickTree

using AbstractTrees
export Node, NewickData
export isroot, isleaf, postwalk, prewalk, children
export getroot, getlca, getleaves, nodefilter
export insertnode!, print_tree, readnw, writenw
export distance, name, id, nwstr, degree

# defaults
name(x) = string(x)
distance(x) = NaN
support(x) = NaN

"""
    NewickData{T,S<:AbstractString}

A simple container for the alowed fields in a newick tree. Those are
the `distance` (expected number of substitutions, time, what have you),
`support` (e.g. bootstrap support value, posterior clade probability)
and `name` (leaf names). Note that internal node labels are not allowed.
"""
mutable struct NewickData{T,S}
    distance::T
    support ::T
    name    ::S
end

NewickData(; d=NaN, s=NaN, n="") = NewickData(promote(d, s)..., n)
name(n::NewickData) = n.name
support(n::NewickData) = n.support
distance(n::NewickData) = n.distance
setname!(n::NewickData, s) = n.name = string(s)
setsupport!(n::NewickData, x) = n.support = x
setdistance!(n::NewickData, x) = n.distance = x

"""
    Node{I,T}

A node of a (phylogenetic) tree.
"""
mutable struct Node{I,T}
    id      ::I
    data    ::T
    parent  ::Node{I,T}
    children::Vector{Node{I,T}}

    Node(id::I, data::T) where {I,T} = new{I,T}(id, data)
    Node(id::I, data::T, p::Nothing) where {I,T} = new{I,T}(id, data)
    function Node(id::I, data::T, p) where {I,T}
        n = new{I,T}(id, data, p)
        push!(p, n)
        return n
    end
end

Node(i; kwargs...) = Node(i, NewickData(; kwargs...))
Node(i, p::Node; kwargs...) = Node(i, NewickData(; kwargs...), p)

id(n::Node) = n.id
data(n::Node) = n.data
name(n::Node) = name(n.data)
support(n::Node) = support(n.data)
distance(n::Node) = distance(n.data)

setname!(n::Node, s) = setname!(n.data, s)
setsupport!(n::Node, s) = setsupport(n.data, s)
setdistance!(n::Node, d) = setdistance!(n.data, d)

isroot(n::Node) = !isdefined(n, :parent)
isleaf(n::Node) = !isdefined(n, :children) || length(n) == 0
degree(n::Node) = isleaf(n) ? 0 : length(n.children)
nv(n::Node) = length(prewalk(n))

nodefilter(n::Node{I}, i::I) where I = filter(x->id(x) == i, prewalk(n))
nodefilter(n::Node, f::Base.Callable) = filter(f, prewalk(n))

Base.pop!(n::Node) = pop!(n.children)
Base.last(n::Node) = last(children(n))
Base.first(n::Node) = first(children(n))
Base.length(n::Node) = length(n.children)
Base.eltype(::Type{Node{T}}) where T = Node{T}

Base.getindex(n::Node{I,T}, i::Integer) where {I,T} = n.children[i]
Base.setindex!(n::Node{I,T}, x::T, i::Integer) where {I,T} = n.children[i] = x

Base.parent(n::Node) = isdefined(n, :parent) ? n.parent : nothing
Base.parent(root, n::Node) = isdefined(n, :parent) ? n.parent : nothing

Base.delete!(n::Node{I}, i::I) where I =
    deleteat!(n.children, findfirst(c->id(c) == i, children(n)))
Base.delete!(n::Node, m::Node) =
    deleteat!(n.children, findfirst(c->c == m, children(n)))

function Base.push!(n::Node, m::Node)
    isdefined(n, :children) ? push!(n.children, m) : n.children = [m]
    m.parent = n
end

function getroot(n::Node)
    while !isroot(n) n = parent(n) end
    return n
end

Base.show(io::IO, n::Node{I,<:NewickData}) where I = write(io, "$(nwstr(n))")
Base.show(io::IO, n::Node) = write(io,"Node($(id(n)), $(n.data))")

# AbstractTrees interface
AbstractTrees.children(n::Node) = isdefined(n, :children) ?
    n.children : typeof(n)[]
Base.eltype(::Type{<:TreeIterator{Node{I,T}}}) where {I,T} = Node{I,T}
AbstractTrees.nodetype(::Type{Node{I,T}}) where {I,T} = Node{I,T}

# AbstractTrees traits
# Base.IteratorEltype(::Type{<:TreeIterator{Node{I,T}}}) where {I,T} = Base.HasEltype()
# AbstractTrees.parentlinks(::Type{Node{I,T}}) where {I,T} = nothing #AbstractTrees.StoredParents()

# Recursive traversals
"""
    postwalk(n)
"""
function postwalk(t)
    ns = typeof(t)[]
    function walk!(n)
        !isleaf(n) && for c in n.children walk!(c) end
        push!(ns, n)
    end
    walk!(t)
    ns
end

"""
    prewalk(n)
"""
function prewalk(t)
    ns = typeof(t)[]
    function walk!(n)
        push!(ns, n)
        !isleaf(n) && for c in n.children walk!(c) end
    end
    walk!(t)
    ns
end

# Path connecting two nodes
function getpath(n::Node, m::Node)
    p1 = getpath(n)
    p2 = getpath(m)
    mrca = intersect(p1, p2)[1]
    i1 = findfirst(x->x==mrca, p1)
    i2 = findfirst(x->x==mrca, p2)
    p1[1:i1], p2[1:i2]
end

# distance between two nodes
function getdistance(n::Node, m::Node)
    p1, p2 = getpath(n, m)
    sum(distance.(p1[1:end-1])) + sum(distance.(p2[1:end-1]))
end

# path connecting node to the root
function getpath(n::Node)
    path = typeof(n)[]
    while !isnothing(n)
        push!(path, n)
        n = parent(n)
    end
    path
end

height(n::Node) = distance.(getpath(n)) |> x->filter(!isnan, x) |> sum

"""
    extract(n::Node, l::AbstractVector{String})

Extract the tree with leaves in `l` from a given tree, preserving
distances if relevant.
"""
function extract(n::Node{I,T}, l::AbstractVector) where {I,T}
    function walk(n)
        if isleaf(n)
            return name(n) ∈ l ? deepcopy(n) : nothing
        else
            below = Node{I,T}[]
            for c in children(n)
                m = walk(c)
                !isnothing(m) && push!(below, m)
            end
            if length(below) == 0
                return nothing
            elseif length(below) == 1
                setdistance!(below[1].data,
                    distance(below[1]) + distance(n))
                return below[1]
            else
                m = deepcopy(n)
                m.children = below
                for c in below c.parent = m end
                return m
            end
        end
    end
    walk(n)
end

"""
    insertnode!(n, m)

Insert node `m` between `n` and its parent at a distance `distance(n) -
distance(m)` from node `n`.
"""
function insertnode!(n::Node{I,T}, m::Node{I,T}) where {I,T}
    a = parent(n)
    delete!(a, n); push!(a, m); push!(m, n)
    setdistance!(n.data, distance(n) - distance(m))
end

"""
    insertnode!(n, dist=NaN, name="")

Insert a new node with name `name` above `n` at a distance `dist` from `n`.
"""
function insertnode!(n::Node{I,<:NewickData}; dist=NaN, name="") where I
    i = maximum(id.(postwalk(getroot(n)))) + 1
    dist = isnan(dist) ? distance(n) / 2 : dist
    insertnode!(n, Node(I(i), NewickData(d=dist, n=name)))
end

getlca(n::Node, a) = getlca(n, a, a)
function getlca(n::Node, a, b)
    clade = getleaves(n)
    m = clade[findfirst(x->name(x)==a, clade)]
    while !(b ∈ name.(getleaves(m)))
        m = parent(m)
    end
    return m
end

function getleaves(n::N) where N<:Node  # mostly faster than Leaves...
    xs = N[]
    for node in postwalk(n)
        isleaf(node) && push!(xs, node)
    end
    xs
end

# Parsing and I/O
"""
    readnw(s::AbstractString, I::Type)

Read a newick string to a tree. Supports the original Newick standard
(http://evolution.genetics.washington.edu/phylip/newicktree.html). One can
have either support values for internal nodes or a node label, but not both.
"""
readnw(s::AbstractString, I::Type=UInt16) =
    try
        readnw(IOBuffer(s), I)
    catch EOFError
        more = s[end] != ";" ? "(no trailing semicolon?)" : ""
        more = !ispath(s) ? more : "(path-like arg instead of Newick string?)"
        throw("Malformed Newick string '$s' $more")
    end

"""
    nwstr(n::Node{I,N}; internal=false)

Generate a newick tree string for the tree rooted in `n`. To make this
work, `N` (the type of `n.data`) should implement `name()` and `distance()`
methods, and optionally `support()`. If `support` is implemented for the
data type and it is not NaN, it will supersede `name` for the labeling of
internal nodes. See for instance the `NewickData` type.
"""
function nwstr(n; internal=false, dist=true)
    function walk(n)
        d = dist ? stringify(':', distance(n)) : ""
        isleaf(n) && return "$(name(n))$d"
        s = join([walk(c) for c in children(n)], ",")
        sv = hasmethod(support, Tuple{typeof(n)}) ?
            stringify(support(n)) : ""
        sv = sv == "" && internal ? name(n) : sv
        d = dist ? stringify(':', distance(n)) : ""
        return "($s)$sv$d"
    end
    s = walk(n)
    s*";"
end

stringify(x) = isnan(x) ? "" : string(x)
stringify(x, y) = isnan(y) ? "" : string(x, y)

"""
    writenw(io::IO, n)
    writenw(fname::AbstractString, n)

Write a newick representation of the tree rooted in `n`. Note that
the tree data type should allow `nwstr(n)` to work. See the `nwstr`
docstring.
"""
writenw(io::IO, n) = write(io, nwstr(n))
writenw(fname::AbstractString, n) = write(fname, nwstr(n) * "\n")

function readnw(io::IOBuffer, I::Type=UInt16)
    i = I(1)
    c = read(io, Char)
    stack = []
    currdata = NewickData()
    while c != ';'
        if c == '('
            push!(stack, Node(i, NewickData())); i += one(i)
            c = read(io, Char)
        elseif c == ')' || c == ','
            target = pop!(stack)
            source = last(stack)
            push!(source, target)
            target.data = currdata
            if c == ')'
                c = read(io, Char)
                eof(io) || c == ';' ? (break) :
                    (currdata, c) = get_nodedata(io, c)
            else
                c = read(io, Char)
            end
        elseif isspace(c)
            c = read(io, Char)
        else
            push!(stack, Node(i, NewickData())); i += one(i)
            leafname, c = get_leafname(io, c)
            currdata, c = get_nodedata(io, c, leafname)
        end
    end
    last(stack)
end

function get_leafname(io::IOBuffer, c)
    leafname = ""
    while !_isnwdelim(c)
        leafname *= c
        c = read(io, Char)
    end
    String(strip(leafname)), c
end

function get_nodedata(io::IOBuffer, c, name="")
    # get everything up to the next comma or )
    support, c = _readwhile!(io, c)
    distance = ""
    if c == ':'
        c = read(io, Char)
        distance, c = _readwhile!(io, c)
    end
    sv = nanparse(support)
    if typeof(sv) == String
        name = String(strip(sv))
        sv = NaN
    end
    NewickData(nanparse(distance), sv, name), c
end

function _readwhile!(io::IOBuffer, c)
    out = ""
    while !_isnwdelim(c)
        out *= c
        c = read(io, Char)
    end
    out, c
end

_isnwdelim(c::Char) = c == ',' || c == ')' || c == ':' || c == ';'

function nanparse(x)
    y = tryparse(Float64, x)
    isnothing(y) ? (x == "" ? NaN : x) : parse(Float64, x)
end

end # module
