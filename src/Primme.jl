module Primme

const libprimme = joinpath(dirname(@__FILE__()), "../deps/primme/lib/libprimme")

const PRIMME_INT = Int # might be wrong. Should be detected.

@enum(Target,
    smallest,        # leftmost eigenvalues */
    largest,         # rightmost eigenvalues */
    closest_geq,     # leftmost but greater than the target shift */
    closest_leq,     # rightmost but less than the target shift */
    closest_abs,     # the closest to the target shift */
    largest_abs      # the farthest to the target shift */
)

@enum(Init,         # Initially fill up the search subspace with: */
    init_default,
    init_krylov, # a) Krylov with the last vector provided by the user or random */
    init_random, # b) just random vectors */
    init_user    # c) provided vectors or a single random vector */
)

@enum(Projection,
    proj_default,
    proj_RR,          # Rayleigh-Ritz */
    proj_harmonic,    # Harmonic Rayleigh-Ritz */
    proj_refined      # refined with fixed target */
)
const C_projection_params = Projection

@enum(Restartscheme,
    thick,
    dtr
)

abstract type PrimmeCStruct end

struct C_restarting_params <: PrimmeCStruct
    scheme::Restartscheme
    maxPrevRetain::Cint
end

struct JD_projectors
    LeftQ::Cint
    LeftX::Cint
    RightQ::Cint
    RightX::Cint
    SkewQ::Cint
    SkewX::Cint
end

@enum(Convergencetest,
    full_LTolerance,
    decreasing_LTolerance,
    adaptive_ETolerance,
    adaptive
)

struct C_correction_params <: PrimmeCStruct
    precondition::Cint
    robustShifts::Cint
    maxInnerIterations::Cint
    projectors::JD_projectors
    convTest::Convergencetest
    relTolBase::Cdouble
end

struct C_stats <: PrimmeCStruct
    numOuterIterations::PRIMME_INT
    numRestarts::PRIMME_INT
    numMatvecs::PRIMME_INT
    numPreconds::PRIMME_INT
    numGlobalSum::PRIMME_INT         # times called globalSumReal
    volumeGlobalSum::PRIMME_INT      # number of SCALARs reduced by globalSumReal
    numOrthoInnerProds::Cdouble      # number of inner prods done by Ortho
    elapsedTime::Cdouble
    timeMatvec::Cdouble              # time expend by matrixMatvec
    timePrecond::Cdouble             # time expend by applyPreconditioner
    timeOrtho::Cdouble               # time expend by ortho
    timeGlobalSum::Cdouble           # time expend by globalSumReal
    estimateMinEVal::Cdouble         # the leftmost Ritz value seen
    estimateMaxEVal::Cdouble         # the rightmost Ritz value seen
    estimateLargestSVal::Cdouble     # absolute value of the farthest to zero Ritz value seen
    maxConvTol::Cdouble              # largest norm residual of a locked eigenpair
    estimateResidualError::Cdouble   # accumulated error in V and W
end

struct C_params <: PrimmeCStruct

    # The user must input at least the following two arguments
    n::PRIMME_INT
    matrixMatvec::Ptr{Void}
    # void (*matrixMatvec)
       # ( void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize,
         # struct primme_params *primme, int *ierr);

    # Preconditioner applied on block of vectors (if available)
    applyPreconditioner::Ptr{Void}
    # void (*applyPreconditioner)
       # ( void *x, PRIMME_INT *ldx,  void *y, PRIMME_INT *ldy, int *blockSize,
         # struct primme_params *primme, int *ierr);

    # Matrix times a multivector for mass matrix B for generalized Ax = xBl
    massMatrixMatvec::Ptr{Void}
    # void (*massMatrixMatvec)
       # ( void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize,
         # struct primme_params *primme, int *ierr);

    # input for the following is only required for parallel programs */
    numProcs::Cint
    procID::Cint
    nLocal::PRIMME_INT
    commInfo::Ptr{Void}
    globalSumReal::Ptr{Void}
    # void (*globalSumReal)
       # (void *sendBuf, void *recvBuf, int *count, struct primme_params *primme,
        # int *ierr );

    # Though Initialize will assign defaults, most users will set these
    numEvals::Cint
    target::Target
    numTargetShifts::Cint             # For targeting interior epairs,
    targetShifts::Ptr{Cdouble}        # at least one shift must also be set

    # the following will be given default values depending on the method
    dynamicMethodSwitch::Cint
    locking::Cint
    initSize::Cint
    numOrthoConst::Cint
    maxBasisSize::Cint
    minRestartSize::Cint
    maxBlockSize::Cint
    maxMatvecs::PRIMME_INT
    maxOuterIterations::PRIMME_INT
    intWorkSize::Cint
    realWorkSize::Csize_t
    iseed::NTuple{4,PRIMME_INT}
    intWork::Ptr{Cint}
    realWork::Ptr{Void}
    aNorm::Cdouble
    eps::Cdouble

    printLevel::Cint
    outputFile::Ptr{Void}

    matrix::Ptr{Void}
    preconditioner::Ptr{Void}
    ShiftsForPreconditioner::Ptr{Cdouble}
    initBasisMode::Init
    ldevecs::PRIMME_INT
    ldOPs::PRIMME_INT

    projectionParams::C_projection_params
    restartingParams::C_restarting_params
    correctionParams::C_correction_params
    stats::C_stats

    convTestFun::Ptr{Void}
    # void (*convTestFun)(double *eval, void *evec, double *rNorm, int *isconv, 
          # struct primme_params *primme, int *ierr);
    convtest::Ptr{Void}
    monitorFun::Ptr{Void}
    # void (*monitorFun)(void *basisEvals, int *basisSize, int *basisFlags,
       # int *iblock, int *blockSize, void *basisNorms, int *numConverged,
       # void *lockedEvals, int *numLocked, int *lockedFlags, void *lockedNorms,
       # int *inner_its, void *LSRes, primme_event *event,
       # struct primme_params *primme, int *err);
    monitor::Ptr{Void}
end

@enum(Svds_target,
    svds_largest,
    svds_smallest,
    svds_closest_abs
)

@enum(Svds_operator,
    svds_op_none,
    svds_op_AtA,
    svds_op_AAt,
    svds_op_augmented
)

struct C_svds_stats <: PrimmeCStruct
    numOuterIterations::PRIMME_INT
    numRestarts::PRIMME_INT
    numMatvecs::PRIMME_INT
    numPreconds::PRIMME_INT
    numGlobalSum::PRIMME_INT         # times called globalSumR
    volumeGlobalSum::PRIMME_INT      # number of SCALARs reduced by globalSumReal
    numOrthoInnerProds::Cdouble      # number of inner prods done by Ortho
    elapsedTime::Cdouble
    timeMatvec::Cdouble              # time expend by matrixMatvec
    timePrecond::Cdouble             # time expend by applyPreconditioner
    timeOrtho::Cdouble               # time expend by ortho
    timeGlobalSum::Cdouble           # time expend by globalSumReal
end

struct C_svds_params <: PrimmeCStruct
    # Low interface: configuration for the eigensolver
    primme::C_params # Keep it as first field to access primme_svds_params from
                          # primme_params
    primmeStage2::C_params # other primme_params, used by hybrid

    # Specify the size of the rectangular matrix A
    m::PRIMME_INT # number of rows
    n::PRIMME_INT # number of columns

    # High interface: these values are transferred to primme and primmeStage2 properly
    matrixMatvec::Ptr{Void}
    # void (*matrixMatvec)
    #    (void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize,
    #     int *transpose, struct primme_svds_params *primme_svds, int *ierr);
    applyPreconditioner::Ptr{Void}
    # void (*applyPreconditioner)
       # (void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize,
        # int *transpose, struct primme_svds_params *primme_svds, int *ierr);

    # Input for the following is only required for parallel programs
    numProcs::Cint
    procID::Cint
    mLocal::PRIMME_INT
    nLocal::PRIMME_INT
    commInfo::Ptr{Void}
    globalSumReal::Ptr{Void}
    # void (*globalSumReal)
       # (void *sendBuf, void *recvBuf, int *count,
        # struct primme_svds_params *primme_svds, int *ierr);

    # Though primme_svds_initialize will assign defaults, most users will set these
    numSvals::Cint
    target::Svds_target
    numTargetShifts::Cint  # For primme_svds_augmented method, user has to
    targetShifts::Ptr{Cdouble} # make sure  at least one shift must also be set
    method::Svds_operator # one of primme_svds_AtA, primme_svds_AAt or primme_svds_augmented
    methodStage2::Svds_operator # hybrid second stage method; accepts the same values as method */

    # These pointers are not for users but for d/zprimme_svds function
    intWorkSize::Cint
    realWorkSize::Csize_t
    intWork::Ptr{Cint}
    realWork::Ptr{Void}

    # These pointers may be used for users to provide matrix/preconditioner
    matrix::Ptr{Void}
    preconditioner::Ptr{Void}

    # The following will be given default values depending on the method
    locking::Cint
    numOrthoConst::Cint
    aNorm::Cdouble
    eps::Cdouble

    precondition::Cint
    initSize::Cint
    maxBasisSize::Cint
    maxBlockSize::Cint
    maxMatvecs::PRIMME_INT
    iseed::NTuple{4,PRIMME_INT}
    printLevel::Cint
    outputFile::Ptr{Void}
    stats::C_svds_stats

    monitorFun::Ptr{Void}
    # void (*monitorFun)(void *basisSvals, int *basisSize, int *basisFlags,
       # int *iblock, int *blockSize, void *basisNorms, int *numConverged,
       # void *lockedSvals, int *numLocked, int *lockedFlags, void *lockedNorms,
       # int *inner_its, void *LSRes, primme_event *event, int *stage,
       # struct primme_svds_params *primme_svds, int *err);
    monitor::Ptr{Void}
end


# Julia API

free(r::Ref{C_svds_params}) = ccall((:primme_svds_free, libprimme), Void, (Ptr{C_svds_params},), r)

function svds_initialize()
    r = Ref{C_svds_params}()
    ccall((:primme_svds_initialize, libprimme), Void, (Ptr{C_svds_params},), r)
    finalizer(r, free)
    return r
end

# matrix-vector product, y = a * x (or y = a^t * x), where
# (void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize,
       # int *transpose, struct primme_svds_params *primme_svds, int *ierr);
const _A_ = Base.Ref{Any}()
function matrixMatvec(xp, ldxp, yp, ldyp, blockSizep, trp, parp, ierrp)
    ldx, ldy           = unsafe_load(ldxp), unsafe_load(ldyp)
    blockSize, tr, par = Int(unsafe_load(blockSizep)), unsafe_load(trp), unsafe_load(parp)
    x, y = unsafe_wrap(Array, xp, (ldx, blockSize)), unsafe_wrap(Array, yp, (ldy, blockSize))

    if tr == 0
        A_mul_B!( view(y, 1:par.m, :), _A_[], view(x, 1:par.n, :))
    else
        Ac_mul_B!(view(y, 1:par.n, :), _A_[], view(x, 1:par.m, :))
    end
    unsafe_store!(ierrp, 0)
    return nothing
end
_fp_ = cfunction(matrixMatvec, Void,
        (Ptr{Float64}, Ptr{Int}, Ptr{Float64}, Ptr{Int}, Ptr{Cint}, Ptr{Cint},
         Ptr{C_svds_params}, Ptr{Cint}))

_print(r::Ref{C_svds_params}) = ccall((:primme_svds_display_params, libprimme), Void, (C_svds_params,), r[])

function Base.setindex!(r::Ref{T}, x, sym::Symbol) where T<:PrimmeCStruct
    p  = Base.unsafe_convert(Ptr{T}, r)
    pp = convert(Ptr{UInt8}, p)
    i  = findfirst(t -> t == sym, fieldnames(r[]))
    o  = fieldoffset(T, i)
    if o == 0
        throw(ArgumentError("no such field"))
    end
    S = fieldtype(T, i)
    unsafe_store!(convert(Ptr{S}, pp + o), x)
    return x
end

function _svds(r::Ref{C_svds_params})
    m, n, k = r[].m, r[].n, r[].numSvals
    svals  = Vector{Float64}(k)
    svecs  = rand(Float64, m + n, k)
    rnorms = Vector{Float64}(k)

    err = ccall((:dprimme_svds, libprimme), Cint,
        (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{C_svds_params}),
         svals, svecs, rnorms, r)
    if err != 0
        warn("err = $err")
    end

    nConv = Int(r[].initSize)

    return reshape(svecs[r[].numOrthoConst*m + (1:(m*nConv))], m, nConv),
        svals,
        reshape(svecs[(r[].numOrthoConst + nConv)*m + r[].numOrthoConst*n + (1:(n*nConv))], n, nConv)
end

function svds(A::AbstractMatrix, k = 5; tol = 1e-12, maxBlockSize = 2k, debuglevel::Int = 0, method::Svds_operator = svds_op_AtA)
    r = svds_initialize()
    _A_[]            = A
    r[:m]            = size(A, 1)
    r[:n]            = size(A, 2)
    r[:matrixMatvec] = _fp_
    r[:numSvals]     = k
    r[:printLevel]   = debuglevel
    r[:eps]          = tol
    r[:maxBlockSize] = maxBlockSize
    r[:method]       = method
    if debuglevel > 0
        _print(r)
    end
    out = _svds(r)
    if debuglevel > 0
        _print(r)
    end
    return out
end

end # module