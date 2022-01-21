module NewickTree

using AbstractTrees
using RecipesBase
export Node, NewickData
export isroot, isleaf, postwalk, prewalk, children, sister
export getroot, getlca, getleaves, nodefilter
export insertnode!, print_tree, readnw, writenw, @nw_str
export distance, name, id, nwstr, degree, getheights
export set_outgroup!, set_outgroup

include("node.jl")
include("parser.jl")
include("nexus.jl")
include("layout.jl")

end # module
