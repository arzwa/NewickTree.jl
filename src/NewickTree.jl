module NewickTree

using AbstractTrees
export Node, isroot, isleaf, postwalk, prewalk, children, print_tree, id, nwstr
export readnw, writenw, distance

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

    function Node(id::I, data::T, p) where {I,T}
        n = new{I,T}(id, data, p)
        isdefined(p, :children) ?
            push!(p.children, n) :
            p.children = [n]
        return n
    end
end

id(n::Node) = n.id
data(n::Node) = n.data
name(n::Node) = name(n.data)
support(n::Node) = support(n.data)
distance(n::Node) = distance(n.data)
isroot(n::Node) = !isdefined(n, :parent)
isleaf(n::Node) = !isdefined(n, :children)
degree(n::Node) = isleaf(n) ? 0 : length(n.children)

Base.parent(n::Node) = n.parent
Base.parent(root, n::Node) = isdefined(n, :parent) ? n.parent : nothing
Base.eltype(::Type{Node{T}}) where T = Node{T}
Base.first(n::Node) = first(children(n))
Base.last(n::Node) = last(children(n))
Base.length(n::Node) = length(n.children)
function Base.delete!(n::Node{I}, i::I) where I
    cs = children(n)
    deleteat!(n, findfirst(c->id(c) == i, cs))
end

function Base.push!(n::Node, m::Node)
    isdefined(n, :children) ? push!(n.children, m) : n.children = [m]
    m.parent = n
end

AbstractTrees.children(n::Node) = isdefined(n, :children) ? n.children : typeof(n)[]
Base.eltype(::Type{<:TreeIterator{Node{I,T}}}) where {I,T} = Node{I,T}

# AbstractTrees traits
# Base.IteratorEltype(::Type{<:TreeIterator{Node{I,T}}}) where {I,T} = Base.HasEltype()
# AbstractTrees.parentlinks(::Type{Node{I,T}}) where {I,T} = nothing #AbstractTrees.StoredParents()
AbstractTrees.nodetype(::Type{Node{I,T}}) where {I,T} = Node{I,T}
Base.show(io::IO, n::Node{I,T}) where {I,T} = write(io, "Node($(id(n)), $(data(n)))")

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

"""
    NewickData{T,S<:AbstractString}

A simple container for the alowed fields in a newick tree. Those are
the `distance` (expected number of substitutions, time, what have you),
`support` (e.g. bootstrap support value, posterior clade probability)
and `name` (leaf names). Note that internal node labels are not allowed.
"""
struct NewickData{T,S}
    distance::T
    support ::T
    name    ::S
end

NewickData(d=NaN, s=NaN, n="") = NewickData(d, s, n)
name(n::NewickData) = n.name
support(n::NewickData) = n.support
distance(n::NewickData) = n.distance

"""
    readnw(s::AbstractString, I::Type)

Read a newick string to a tree.
"""
readnw(s::AbstractString, I::Type=UInt16) = try
        readnw(IOBuffer(s), I)
    catch EOFError
        throw("Malformed Newick string '$s'")
    end

"""
    nwstr(n::Node{I,N})

Generate a newick tree string for the tree rooted in `n`. To make this
work, `N` (the type of `n.data`) should implement `name()`, `distance()`
and `support()` functions. See for instance the `NewickData` type.
"""
function nwstr(n)
    function walk(n)
        d = isnan(distance(n)) ? "" : ":$(distance(n))"
        isleaf(n) && return "$(name(n))$d"
        s = join([walk(c) for c in children(n)], ",")
        sv = isnan(support(n)) ? "" : support(n)
        return "($s)$sv$d"
    end
    s = walk(n)
    s[1:findlast(')', s)]*";"
end

"""
    writenw(io::IO, n)
    writenw(fname::AbstractString, n)

Write a newick representation of the tree rooted in `n`. Note that
the tree data type should allow `nwstr(n)` to work. See the `nwstr`
docstring.
"""
writenw(io::IO, n) = write(io, nwstr(n))
writenw(fname::AbstractString, n) = write(fname, nwstr(n))

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
                eof(io) ? (break) : (currdata, c) = get_nodedata(io, c)
            else
                c = read(io, Char)
            end
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
    leafname, c
end

function get_nodedata(io::IOBuffer, c, name="")
    # get everything up to the next comma or )
    support, c = _readwhile!(io, c)
    distance = ""
    if c == ':'
        c = read(io, Char)
        distance, c = _readwhile!(io, c)
    end
    NewickData(nanparse(distance), nanparse(support), name), c
end

function _readwhile!(io::IOBuffer, c)
    out = ""
    while !_isnwdelim(c)
        out *= c
        c = read(io, Char)
    end
    out, c
end

_isnwdelim(c::Char) = c == ',' || c == ')' || c == ':'

function nanparse(x)
    y = tryparse(Float64, x)
    isnothing(y) ? NaN : parse(Float64, x)
end

end # module
