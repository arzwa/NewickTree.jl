# # NewickTree.jl

# Read a [newick tree](http://evolution.genetics.washington.edu/phylip/newicktree.html)
# to a tree data structure. The resulting data structure supports the
# [AbstractTrees.jl](https://github.com/JuliaCollections/AbstractTrees.jl)
# interface.

# ## Reading trees
using NewickTree
t = readnw("((A:1.2,B:1.4)86:0.2,C:0.6);")
print_tree(t)
writenw(stdout, t)

# Use `readnw(readline("your_file.nw"))` to read a newick tree from a file. Use
# `readnw.(readlines("your_file.nw"))` to read a vector of trees from a file
# with a newick tree on each line.

# Note that trees should adhere to the Newick standard, they should end with a
# semicolon and can only contain leaf names, support values and branch lengths
# as node information. Failure to provide a valid Newick string will trigger
# an error:
try
    t = readnw("((A:1.2,B:1.4)86:0.2,C:0.6)")
catch ex
    @show ex
end

# ## Support for writing other tree structured data to Newick strings

# Any data structure that implements the AbstractTrees interface (i.e. defines
# `AbstractTrees.children`) can be written to a Newick structure provided several
# functions re defined. For example:
using AbstractTrees
t = (((1,2),3),(4,5))
print_tree(t)

# The following functions should be defined
NewickTree.name(x::Int) = x
NewickTree.support(x::Union{Int,Tuple}) = NaN
NewickTree.distance(x::Union{Int,Tuple}) = NaN
NewickTree.isleaf(x) = typeof(x) == Int ? true : false

# This enables us to use the `writenw` function
io = IOBuffer()
writenw(io, t)
s = String(take!(io))

# now we can read the Newick string
n = readnw(s)
print_tree(n)
