using Test
using ConstantFunctionMarketMakers

@testset "Test checks" begin
    @test 1 + 1 == 2
end

@testset "ECLP checks" begin
    amm = ECLP(0.9, 1.1, 1.0, 10.0)

    # Regression where this failed b/c I forgot to transpose a vector.
    @test begin
        tc_t(amm, [1.0, 1.0])
        true
    end
end
