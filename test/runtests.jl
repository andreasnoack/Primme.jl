using Base.Test
using Primme

@testset "svds" begin
    @testset "svds m=$m, n=$n, k=$k" for n in [200, 400],
                                         m in [200, 400],
                                         k = [10,20],
                                         method = [Primme.svds_op_AAt, Primme.svds_op_AtA, Primme.svds_op_augmented]
        A = randn(m, n)
        svdPrimme = Primme.svds(A, k, method = method)
        svdLAPACK = svd(A)
        @test svdLAPACK[2][1:k] ≈ svdPrimme[2]
        @test abs.(svdLAPACK[1][:, 1:k]'svdPrimme[1]) ≈ eye(k)
    end
end