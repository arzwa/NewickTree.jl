using Test
using NewickTree, AbstractTrees

@testset "Tree reading" begin
    n = readnw(IOBuffer("((A:1,B:2)90:1,C:1);"))
    @test length(collect(Leaves(n))) == 3
    @test distance(n) === NaN
    @test distance(first(n)) == 1.
    @test name(last(n)) == "C"

    n = readnw(readline(joinpath(@__DIR__, "ncov-spike.nw")))
    @test map(x->(name(x), distance(x)), Leaves(n))[3][2] == 0.182372964
    @test name(postwalk(n)[12]) == "Human_coronavirus_OC43__YP_009555241.1"

    t = readnw("((bears : 10 , dogs :1)  carnivores, the penguins are here );")
    T = Tree(t)
    @test name(T[1][2]) == "dogs"
    @test name(T[2]) == "the penguins are here"

    t = readnw("( (a : 10 , b : 1) 97 : 8, hello\nlobster\t )  ;  ")
    @test name(t[2]) == "hello\nlobster"
    @test NewickTree.support(t[1]) == 97
end
