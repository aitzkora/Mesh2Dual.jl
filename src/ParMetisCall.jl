"""
parmetis\\_mesh\\_to\\_dual(;elmdist, eptr, eind, baseval, ncommon, comm)

call the PARMETIS_Mesh2Dual routine which computes the dual graph of a mesh
using ncommon points to define adjacency relationship between elements.
elmdist is a common integer vector to all process such that
`[elmdist[rank+1],elmdist[rank+2]-1]` is the range of elements of 
the process of rank rank
"""
function parmetis_mesh_to_dual(;elmdist::Array{T,1}, 
                               eptr::Array{T,1}, 
                               eind::Array{T,1}, 
                               baseval::T, 
                               ncommon::T, 
                               comm::MPI.Comm) where {T}
    if "METIS_LIB" in keys(ENV)
        metis_str = ENV["METIS_LIB"]
    else
        metis_str = "/usr/lib/libmetis.so"
    end
    lib_metis = dlopen(metis_str; throw_error=false)
     if "PARMETIS_LIB" in keys(ENV)
        parmetis_str = ENV["PARMETIS_LIB"]
    else
        parmetis_str = "/usr/lib/libparmetis.so"
    end
    lib_parmetis = dlopen(parmetis_str; throw_error=false)
    @debug lib_parmetis
    @assert lib_parmetis != nothing
    @debug lib_metis
    @assert lib_metis != nothing
    mesh_to_dual_ptr = dlsym(lib_parmetis, :ParMETIS_V3_Mesh2Dual)
    @debug "PARMETIS_MeshToDual Pointer", mesh_to_dual_ptr
    r_xadj = Ref{Ptr{T}}()
    r_adjncy = Ref{Ptr{T}}()
    t₀ = time()
    ccall(mesh_to_dual_ptr, Cvoid, 
          (Ptr{T}, Ptr{T}, Ptr{T}, Ref{T}, Ref{T}, Ref{Ptr{T}}, Ref{Ptr{T}}, Ref{MPI.Comm}),
          elmdist,
          eptr,
          eind,
          Ref{T}(baseval),
          Ref{T}(ncommon),
          r_xadj,
          r_adjncy,
          Ref{MPI.Comm}(comm)
         )
    @printf "total time = %.3e\n"  (time() - t₀)
    rank = MPI.Comm_rank(comm)
    ne_local = elmdist[rank+2]-elmdist[rank+1]
    x_adj = GC.@preserve r_xadj [unsafe_load(r_xadj[] ,i) for i=1:ne_local+1]
    x_adjncy = GC.@preserve r_adjncy [unsafe_load(r_adjncy[],i) for i=1:x_adj[end] ]
    return x_adj, x_adjncy
end  
# function parmetis_mesh_to_dual

"""
ptscotchparmetis_mesh_to_dual\\_mesh\\_to\\_dual(;elmdist, eptr, eind, baseval, ncommon, comm)

call the PARMETIS_Mesh2Dual routine from libscotchmetis which computes the dual graph of a mesh
"""
function ptscotchparmetis_mesh_to_dual(;elmdist::Array{T,1}, 
                               eptr::Array{T,1}, 
                               eind::Array{T,1}, 
                               baseval::T, 
                               ncommon::T, 
                               comm::MPI.Comm) where {T}
    if "SCOTCHERR_LIB" in keys(ENV)
        scotcherr_str = ENV["SCOTCHERR_LIB"]
    end
    #lib_scotcherr = dlopen(scotcherr_str, RTLD_LAZY|RTLD_DEEPBIND|RTLD_GLOBAL ; throw_error=false)
    lib_scotcherr = dlopen(scotcherr_str, RTLD_LAZY|RTLD_GLOBAL ; throw_error=false)
    @assert lib_scotcherr != nothing
    if "PARSCOTCHMETIS_LIB" in keys(ENV)
        parmetis_str = ENV["PARSCOTCHMETIS_LIB"]
    end
    lib_ptscotchparmetis = dlopen(parmetis_str; throw_error=false)
    @debug lib_ptscotchparmetis
    @assert lib_ptscotchparmetis != nothing
    mesh_to_dual_ptr = dlsym(lib_ptscotchparmetis, :SCOTCH_ParMETIS_V3_Mesh2Dual)
    @debug "PTSCOTCHPARMETIS_MeshToDual Pointer", mesh_to_dual_ptr
    r_xadj = Ref{Ptr{T}}()
    r_adjncy = Ref{Ptr{T}}()
    t₀ = time()
    ccall(mesh_to_dual_ptr, Cvoid, 
          (Ptr{T}, Ptr{T}, Ptr{T}, Ref{T}, Ref{T}, Ref{Ptr{T}}, Ref{Ptr{T}}, Ref{MPI.Comm}),
          elmdist,
          eptr,
          eind,
          Ref{T}(baseval),
          Ref{T}(ncommon),
          r_xadj,
          r_adjncy,
          Ref{MPI.Comm}(comm)
         )
    @printf "total time = %.3e\n"  (time() - t₀)
    rank = MPI.Comm_rank(comm)
    ne_local = elmdist[rank+2]-elmdist[rank+1]
    @info "ne_local =", ne_local
    x_adj = GC.@preserve r_xadj [unsafe_load(r_xadj[] ,i) for i=1:ne_local+1]
    x_adjncy = GC.@preserve r_adjncy [unsafe_load(r_adjncy[],i) for i=1:x_adj[end] ]
    return x_adj, x_adjncy
end  # function ptscotchparmetis_mesh_to_dual
