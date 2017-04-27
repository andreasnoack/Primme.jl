depdir = dirname(@__FILE__)
srcdir = joinpath(depdir, "primme")

if isdir(srcdir)
    rm(srcdir, recursive=true)
end