module NewickTree

import DataStructures: OrderedSet
import AbstractTrees: children, print_tree, Tree
export TreeNode, PhyloTree, isroot, isleaf, postwalk, prewalk, readnw
export children, print_tree, tonw, extract, id

abstract type AbstractNode end

"""
    TreeNode{I,T}

A node of a (phylogenetic) tree. We used OrderedSet as this allows
deterministic iteration over child nodes while allowing constant-time
insertion and deletion.
"""
mutable struct TreeNode{I,T} <: AbstractNode
    id      ::I                        # index
    x       ::T                        # can hold anything, typically a distance
    parent  ::Union{TreeNode{I,T},Nothing}
    children::OrderedSet{TreeNode{I,T}}
end

const Clade{I,T} = OrderedSet{TreeNode{I,T}} where {I,T}
Clade(nodes::Vector{TreeNode}) = OrderedSet(nodes)
newclade(n::TreeNode{I,T}) where {I,T} = Clade{I,T}()

# root initializer
TreeNode{I,T}() where {I,T} = TreeNode(one(I), zero(T), nothing, Clade{I,T}())
TreeNode(i::I, x::T) where {I,T} = TreeNode(i, x, nothing, Clade{I,T}())
TreeNode(i::I, x::T, p) where {I,T} = TreeNode(i, x, p, Clade{I,T}())

id(n::TreeNode) = n.id
isroot(n::TreeNode) = isnothing(n.parent)
isleaf(n::TreeNode) = length(n.children) == 0
distance(n::TreeNode) = n.x
Base.parent(n::TreeNode) = n.parent

"""
    PhyloTree{I,T}

A phylogenetic tree struct, combining nodes, leaf names and
optionally an additional value for internal nodes such as bootstrap
support values or posterior clade probabilities.
"""
struct PhyloTree{I,T}
    nodes ::Vector{TreeNode{I,T}}
    leaves::Dict{I,String}
    values::Vector{T}
end

Base.show(io::IO, t::PhyloTree{I,T}) where {I,T} = write(io, "PhyloTree{$I,$T}")
print_tree(t::PhyloTree) = print_tree(Tree(t.nodes[1]))
print_tree(n::TreeNode) = print_tree(Tree(n))

# Base extensions
Base.getindex(t::PhyloTree, i) = t.nodes[i]
Base.length(t::TreeNode) = length(t.children)
Base.delete!(t::TreeNode, x::TreeNode) = delete!(t.children, x)
Base.push!(t::TreeNode, x::TreeNode) = push!(t.children, x)
Base.first(t::TreeNode) = first(t.children)
Base.show(io::IO, t::TreeNode{I,T}) where {I,T} =
    write(io, "TreeNode{$I,$T}($(t.id); $(length(t)) children)")

children(n::TreeNode) = collect(n.children)
function leafname end

# io
"""
    readnw(nw_str::String [, index_type, value_type])

Read a phylogenetic tree from a newick string. Supports branch
lengths and support-like values (for internal nodes). Example
of a valid newick string:

    str  = "((A:1.2,B:1.4)86:0.2,C:0.6);"
    tree = readnw(str)
"""
function readnw(nw_str::String, I::Type=UInt16, T::Type=Float64)
    l = Dict{I,String}(); v = Dict{I,T}()
    stack = TreeNode{I,T}[]
    i = I(1); x = T(0.); sv = T(0.); j=I(1)
    while i < length(nw_str)
        if nw_str[i] == '('  # add an internal node to stack & tree
            push!(stack, TreeNode(j, T(0.))); j += I(1)
        elseif nw_str[i] == ')'  # add branch + collect data on node
            target = pop!(stack)
            source = stack[end]
            target.parent = source
            target.x = x
            push!(source.children, target)
            sv, x, i = get_node_info(nw_str, i+1)
            v[source.id] = sv
        elseif nw_str[i] == ','  # add branch
            target = pop!(stack)
            source = stack[end]
            target.parent = source
            target.x = x
            v[target.id] = sv
            push!(source.children, target)
        else  # store leaf name and get branch length
            push!(stack, TreeNode(j, T(0.))); j += I(1)
            leaf, i = get_leaf_name(nw_str, i)
            sv, x, i = get_node_info(nw_str, i)
            l[stack[end].id] = leaf
        end
        i += I(1)
    end
    # prewalk(stack[end])
    PhyloTree(prewalk(stack[end]), l, collect(values(sort(v))))
end

function get_leaf_name(nw_str, i)
    j = i
    while (nw_str[j] != ':') & (nw_str[j] != ',') & (nw_str[j] != ')')
        j += 1
    end
    leaf = nw_str[i:j-1]
    return leaf, j
end

function get_node_info(nw_str, i)
    # get everything up to the next comma or semicolon
    substr = nw_str[i:end]
    substr = split(substr, ',')[1]
    substr = split(substr, ')')[1]
    if substr == ";"
        return -1.0, -1.0, i
    end
    substr = split(substr, ';')[1]
    if occursin(":", substr)
        sv, bl = split(substr, ':')
        if sv == ""  # only branch length
            sv = -1.0
            bl = parse(Float64, bl)
        else  # both (B)SV and branch length are there
            sv = parse(Float64, sv)
            bl = parse(Float64, bl)
        end
    else
        # nothing there
        bl = sv = -1.
    end
    return sv, bl, i + length(substr) -1
end

# Recursive traversals
function postwalk(t)
    ns = typeof(t)[]
    function walk!(n)
        if !isleaf(n)
            for c in n.children
                walk!(c)
            end
        end
        push!(ns, n)
    end
    walk!(t)
    ns
end

function prewalk(t::TreeNode{I,T}) where {I,T}
    ns = typeof(t)[]
    function walk!(n)
        push!(ns, n)
        if !isleaf(n)
            for c in n.children
                walk!(c)
            end
        end
    end
    walk!(t)
    ns
end

function extract(t::PhyloTree{I,T}, ℒ::Vector{String}) where {I,T}
    new_leaves = Dict{I,String}()
    function walk(n)
        if isleaf(n)
            new_leaves[n.id] = t.leaves[n.id]
            return t.leaves[n.id] ∈ ℒ ? deepcopy(n) : nothing
        else
            below = TreeNode{I,T}[]
            for c in n.children
                m = walk(c)
                isnothing(m) ? nothing : push!(below, m)
            end
            if length(below) == 0
                return nothing
            elseif length(below) == 1
                below[1].x += n.x  # increment branch length
                return below[1]
            else
                m = TreeNode(n.id, n.x, nothing, OrderedSet(below))
                for c in below
                    c.parent = m
                end
                return m
            end
        end
    end
    nodes = prewalk(walk(t[1]))
    nodes[1].x = 0.
    PhyloTree(nodes, new_leaves, t.values)
end

Base.write(tree::PhyloTree) = write(stdout, tree)
Base.write(io::IO, tree::PhyloTree) = write(io, tonw(tree))

tonw(tree::PhyloTree) = tonw(tree.nodes[1],
    (n)->tree.leaves[n.id],
    blen=(n)->n.x)

function tonw(node, leafname; blen=(n)->1., label=(n)->"")
    function walk(n)
        if isleaf(n)
            return "$(leafname(n)):$(blen(n))"
        else
            nwstr = join([walk(c) for c in children(n)], ",")
            return !isroot(n) ? "($nwstr)$(label(n)):$(blen(n))" : "($(nwstr));"
        end
    end
    walk(node)
end

end # module
