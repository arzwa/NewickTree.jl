using Test
using NewickTree, AbstractTrees

@testset "NewickTree" begin
    @testset "Tree reading" begin
        n = readnw(IOBuffer("((A:1,B:2)90:1,C:1);"))
        @test length(collect(getleaves(n))) == 3
        @test distance(n) === NaN
        @test distance(first(n)) == 1.
        @test name(last(n)) == "C"

        n = readnw(readline(joinpath(@__DIR__, "ncov-spike.nw")))
        @test map(x->(name(x), distance(x)), getleaves(n))[3][2] == 0.182372964
        @test name(postwalk(n)[12]) == "Human_coronavirus_OC43__YP_009555241.1"

        t = readnw("((bears : 10 , dogs :1)  carnivores, the penguins here );")
        @test name(t[1][2]) == "dogs"
        @test name(t[2]) == "the penguins here"

        t = readnw("( (a : 10 , b : 1) 97 : 8, hello\nlobster\t )  ;  ")
        @test name(t[2]) == "hello\nlobster"
        @test NewickTree.support(t[1]) == 97
    end

    @testset "AbstractTrees interface" begin
        n = readnw(readline(joinpath(@__DIR__, "ncov-spike.nw")))
        @test Set(collect(Leaves(n))) == Set(getleaves(n))
        # t = Tree(n) # doesn't exist in AbstractTrees@0.4
        # @test t[2][1] == n[2][1]
        @test all(collect(PostOrderDFS(n)) .== postwalk(n))
        @test all(collect(PreOrderDFS(n)) .== prewalk(n))
        # NOTE: the order doesn't have to be exactly the same, but
        # turns out it is.
    end

    @testset "LCA and paths" begin
        t = readnw("(((atr:2.47,(osa:1.82,vvi:1.82):0.65):0.91,(gbi:2.90,"*
            "(gmo:1.77,wmi:1.77):1.13):0.48):0.80,(afi:0.87,scu:0.87):3.31);")
        for n in getleaves(t)
            @test NewickTree.height(n) ≈ 4.18
        end
        @test length(getleaves(getlca(t, "atr", "gbi"))) == 6
        @test NewickTree.height(getlca(t, "atr", "osa")) == 1.71
        l, r = NewickTree.getpath(getlca(t, "atr"), getlca(t, "osa"))
        @test length(l) == 2
        @test length(r) == 3
    end

    @testset "Other (non `Node` based) tree structures" begin
        t = (((1,2),3),(4,5))
        NewickTree.isleaf(::Int) = true
        NewickTree.isleaf(::Tuple) = false
        @test nwstr(t) == "(((1,2),3),(4,5));"
        @test nwstr(t[1]) == "((1,2),3);"
        @test nwstr(t[2][1]) == "4;"
    end

    if isdefined(Base, :get_extension)
        try
            using Clustering

            A = [1 1 5 5; 1 1 5 5; 5 5 1 1; 5 5 1 1]
            hc = hclust(A)

            @testset "Converters to `Node`" begin
                @test hasmethod(convert, (Type{NewickTree.Node},Clustering.Hclust))
                @test nwstr(convert(Node, hclust(A))) == "((1:1.0,2:1.0):4.0,(3:1.0,4:1.0):4.0):0.0;"
            end
        catch
            @warn "Clustering not available, skipping extension tests"
        end
    end
end


# Benchmark on ncov-spike.nw (11/06/2020)
# reading
# NewickTree.jl:
# julia> @btime readnw(s);
#   113.429 μs (3671 allocations: 135.02 KiB)
# Phylo.jl:
# julia> @btime Phylo.parsenewick(s);
#   373.154 μs (3151 allocations: 253.85 KiB)

# postorder
# NewickTree.jl
# julia> @btime postwalk(t);
#   1.102 μs (8 allocations: 1.19 KiB)
# Phylo.jl
# julia> @btime traversal(t, postorder);
#   8.715 μs (174 allocations: 7.97 KiB)
