using Base.Test
using Primme

@testset "svds" begin
    @testset "svds m=$m, n=$n, k=$k" for n in [200, 400], m in [200, 400], k = [10,20]
        A = randn(m, n)
        svdPrimme = Primme.svds(A, k)
        svdLAPACK = svd(A)
        @test svdLAPACK[2][1:k] ≈ svdPrimme[2]
        @test abs.(svdLAPACK[1][:, 1:k]'svdPrimme[1]) ≈ eye(k)
    end
end