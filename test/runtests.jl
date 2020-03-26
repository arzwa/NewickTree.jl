using Test
using NewickTree
import NewickTree: name

n = readnw(IOBuffer("((A:1,B:2)90:1,C:1);"))
@test length(collect(Leaves(n))) == 3
@test distance(n) === NaN
@test distance(first(n)) == 1.
@test name(last(n)) == "C"

n = readnw(readline(joinpath(@__DIR__, "ncov-spike.nw")))
@test map(x->(name(x), distance(x)), Leaves(n))[3][2] == 0.182372964
@test name(postwalk(n)[12]) == "Human_coronavirus_OC43__YP_009555241.1"
