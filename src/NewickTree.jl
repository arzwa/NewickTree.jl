module NewickTree

using AbstractTrees
using RecipesBase
export Node, NewickData
export isroot, isleaf, postwalk, prewalk, children
export getroot, getlca, getleaves, nodefilter
export insertnode!, print_tree, readnw, writenw
export distance, name, id, nwstr, degree, getheights

include("node.jl")
include("parser.jl")
include("layout.jl")

end # module
