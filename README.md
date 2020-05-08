# NewickTree.jl

Read a [newick tree](http://evolution.genetics.washington.edu/phylip/newicktree.html) to a tree data structure. The resulting data structure supports the [AbstractTrees.jl](https://github.com/JuliaCollections/AbstractTrees.jl) interface.

## Reading trees

```julia
using NewickTree
t = readnw("((A:1.2,B:1.4)86:0.2,C:0.6);")
print_tree(t)
```

```
((A:1.2,B:1.4)86.0:1.4,C:0.6);
├─ (A:1.2,B:1.4);
│  ├─ A:1.2
│  └─ B:1.4
└─ C:0.6
```

Use `readnw(readline("your_file.nw"))` to read a newick tree from a file. Use `readnw.(readlines("your_file.nw"))` to read a vector of trees from a file with a newick tree on each line.

Note that trees should adhere to the Newick standard, they should end with a semicolon and can only contain (1) leaf names, (2) support values *or* internal names and (3) branch lengths as node information. Failure to provide a valid Newick string will trigger an error:

```julia
try
    t = readnw("((A:1.2,B:1.4)86:0.2,C:0.6)")
catch ex
    @show ex
end
```
```
"Malformed Newick string '((A:1.2,B:1.4)86:0.2,C:0.6)' (no trailing semicolon?)"
```

The tree data structure is pretty straightforward, with nodes storing the following fields:

```julia
fieldnames(typeof(t))
```
```
(:id, :data, :parent, :children)
```

Functions from `AbstractTrees` can be used, for instance

```julia
using AbstractTrees
collect(Leaves(t))
```
```
3-element Array{Node{UInt16,NewickData{Float64,String}},1}:
 A:1.2;
 B:1.4;
 C:0.6;
```

or

```julia
collect(PostOrderDFS(t))
```
```
5-element Array{Node{UInt16,NewickData{Float64,String}},1}:
 A:1.2;
 B:1.4;
 (A:1.2,B:1.4)86.0:0.2;
 C:0.6;
 ((A:1.2,B:1.4)86.0:0.2,C:0.6);
```

some simple recursive tree traversals are also implemented

```julia
postwalk(t)
prewalk(t)
```
```
5-element Array{Node{UInt16,NewickData{Float64,String}},1}:
 ((A:1.2,B:1.4)86.0:0.2,C:0.6);
 (A:1.2,B:1.4)86.0:0.2;
 A:1.2;
 B:1.4;
 C:0.6;
```

these tend to be faster (at least for small trees?)

## Writing trees

`nwstr` converts a tree data structure that implements the required functions (see below) to a Newick string:

```julia
nwstr(t)
```
```
"((A:1.2,B:1.4)86.0:0.2,C:0.6);"
```

`writenw` uses this to write to a stream or file.

```julia
io = IOBuffer()
writenw(io, t)
String(take!(io))
```
```
"((A:1.2,B:1.4)86.0:0.2,C:0.6);"
```

## Support for writing other tree structured data to Newick strings

Any data structure that implements the AbstractTrees interface (i.e. defines `AbstractTrees.children`) can be written to a Newick structure provided several functions are defined. For example:

```julia
using AbstractTrees
t = (((1,2),3),(4,5))
print_tree(t)
```

```
(((1, 2), 3), (4, 5))
├─ ((1, 2), 3)
│  ├─ (1, 2)
│  │  ├─ 1
│  │  └─ 2
│  └─ 3
└─ (4, 5)
   ├─ 4
   └─ 5
```

The following functions should be defined

```julia
NewickTree.name(x::Tuple) = ""
NewickTree.name(x::Int) = string(x)
NewickTree.distance(x::Union{Int,Tuple}) = NaN
NewickTree.isleaf(x) = typeof(x) == Int ? true : false
```

This enables us to use the `nwstr` and `writenw` functions

```julia
s = nwstr(t)
```
```
"(((1,2),3),(4,5));"
```

now we can read the Newick string

```julia
n = readnw(s)
print_tree(n)
```

```
(((1,2),3),(4,5));
├─ ((1,2),3);
│  ├─ (1,2);
│  │  ├─ 1
│  │  └─ 2
│  └─ 3
└─ (4,5);
   ├─ 4
   └─ 5
```

```julia
using Literate
Literate.markdown(
    joinpath(@__DIR__, "README.jl"),
    joinpath(@__DIR__, "../"),
    documenter=false, execute=true)
```
```
"/home/arzwa/dev/NewickTree/README.md"
```

using the execute-markdown branch now

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

