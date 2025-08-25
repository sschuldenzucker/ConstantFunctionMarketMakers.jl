using Test
using ConstantFunctionMarketMakers

import ConstantFunctionMarketMakers.ECLPMath

@testset "Test checks" begin
    @test 1 + 1 == 2
end

@testset "ECLP checks" begin
    amm = ECLP(0.9, 1.1, 1.05, 10.0)

    # Regression where this failed b/c I forgot to transpose a vector.
    @test begin
        tc_t(amm, [1.0, 1.0])
        true
    end

    # Regression
    @testset "Comparison old/new ECLPMath impl" begin
        # TODO more property tests like this. We can also have property tests below.
        let t = [50.0, 50.0]
            @test isapprox(
                ECLPMath.residual_eclp(
                    amm.params,
                    t,
                    ECLPMath.calc_invariant(amm.params, t),
                ),
                0.0;
                atol = 1e-12,
            )
        end

        tc = tc_v_p(amm, 100.0, 1.0)

        let y = 50.0
            @test isapprox(
                ECLPMath.mk_x_y(amm.params, tc.common.l)(y),
                ECLPMath.x_y(amm.params, tc.common.l, y),
            )
        end

        let x = 25.0
            @test isapprox(
                ECLPMath.mk_y_x(amm.params, tc.common.l)(x),
                ECLPMath.y_x(amm.params, tc.common.l, x),
            )
        end
    end
end
