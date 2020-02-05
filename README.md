# NewickTree.jl

Read a newick tree to a tree data structure.

```julia
julia> t = readnw("((A:1.2,B:1.4)86:0.2,C:0.6);")
PhyloTree{UInt16,Float64}

julia> t[1]
TreeNode{UInt16,Float64}(1; 2 children)

julia> print_tree(t)
TreeNode{UInt16,Float64}(1; 2 children)
├─ TreeNode{UInt16,Float64}(2; 2 children)
│  ├─ TreeNode{UInt16,Float64}(3; 0 children)
│  └─ TreeNode{UInt16,Float64}(4; 0 children)
└─ TreeNode{UInt16,Float64}(5; 0 children)

julia> for c in NewickTree.children(t[1])
           print_tree(c)
       end
TreeNode{UInt16,Float64}(2; 2 children)
├─ TreeNode{UInt16,Float64}(3; 0 children)
└─ TreeNode{UInt16,Float64}(4; 0 children)

TreeNode{UInt16,Float64}(5; 0 children)

julia> for n in postwalk(t[1])
           println(n)
       end
TreeNode{UInt16,Float64}(3; 0 children)
TreeNode{UInt16,Float64}(4; 0 children)
TreeNode{UInt16,Float64}(2; 2 children)
TreeNode{UInt16,Float64}(5; 0 children)
TreeNode{UInt16,Float64}(1; 2 children)

julia> isleaf(t[5])
true

julia> tonw(t)
"((A:1.2,B:1.4):0.2,C:0.6);"

julia> write("fname.nw", t);

julia> tt = extract(t, ["A","C"])
PhyloTree{UInt16,Float64}

julia> print_tree(tt)
TreeNode{UInt16,Float64}(1; 2 children)
├─ TreeNode{UInt16,Float64}(3; 0 children)
└─ TreeNode{UInt16,Float64}(5; 0 children)
```
