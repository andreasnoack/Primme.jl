# Use this version of Elemental
ver = "56b6d0838b09718206ce1f33e163ce6486bcb086" # v2.1

if is_windows()
    error("Primme.jl only works on Unix Platforms")
end

depdir = dirname(@__FILE__)

srcdir = joinpath(depdir, "primme")

if !isdir(srcdir)
    LibGit2.clone("https://github.com/andreasnoack/primme.git", "$srcdir") # Use my fork until Makefile has been fixed
    # LibGit2.clone("https://github.com/primme/primme.git", "$srcdir")
end
cd(srcdir) do
    LibGit2.checkout!(LibGit2.GitRepo("."), "$ver")
end

BLAS.check()
blas = BLAS.vendor()
mathlib = Libdl.dlpath(BLAS.libblas)
blas64 = LinAlg.USE_BLAS64 ? "ON" : "OFF"
blas_suffix = blas === :openblas64 ? "_64_" : "_"

cd(srcdir) do
    run(`make -j solib`)
end
