# Use this version of Elemental
ver = "9b2bd87efa69287eaa21d5fd4685b57c6d1f5278" # v2.1

if is_windows()
    error("Primme.jl only works on Unix Platforms")
end

depdir = dirname(@__FILE__)

srcdir = joinpath(depdir, "primme")

if !isdir(srcdir)
    LibGit2.clone("https://github.com/primme/primme.git", "$srcdir") # Use my fork until Makefile has been fixed
    # LibGit2.clone("https://github.com/primme/primme.git", "$srcdir")
else
    LibGit2.fetch(LibGit2.GitRepo("$srcdir"))
end
cd(srcdir) do
    LibGit2.checkout!(LibGit2.GitRepo("."), "$ver")
end

cflags = "-O3 -fPIC"
if LinAlg.USE_BLAS64
    cflags *= " -DPRIMME_BLASINT_SIZE=64"
    if BLAS.vendor() == :openblas64
        cflags *= " -DPRIMME_BLAS_SUFFIX=_64"
    end
end

libpath = dirname(Libdl.dlpath(Base.liblapack_name))

cd(srcdir) do
    run(`make solib CFLAGS=$cflags LDFLAGS="-l$(Base.liblapack_name[4:end]) -l$(Base.libblas_name[4:end]) -L$libpath"`)
end
