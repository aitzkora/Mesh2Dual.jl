using Libdl: dlopen, dlsym

"""
metis\\_graph\\_dual(m::Mesh, n_common::Int)
computes the graph dual (i.e. elements graph using Metis)
"""

function metis_graph_dual(m::Mesh{T}, n_common::Int) where {T}
    n_common = T(n_common)
    if "METIS_LIB" in keys(ENV)
        metis_str = ENV["METIS_LIB"]
    else
        metis_str = "/usr/lib/libmetis.so"
    end
    lib_metis = dlopen(metis_str; throw_error=false)
    @debug lib_metis
    @assert lib_metis != nothing
    grf_dual_ptr = dlsym(lib_metis, :libmetis__CreateGraphDual)
    @debug "CreateGraphDual Pointer", grf_dual_ptr
    eptr, eind, mini_node = mesh_to_metis_fmt(m)
    r_xadj = Ref{Ptr{T}}()
    r_adjncy = Ref{Ptr{T}}()
    ccall(grf_dual_ptr, Cvoid, 
          (Cint, Cint, Ptr{T}, Ptr{T}, Cint, Ref{Ptr{T}}, Ref{Ptr{T}}),
          size(m.elements, 1),
          size(m.nodes, 1),
          eptr,
          eind,
          n_common,
          r_xadj,
          r_adjncy
         )
    x_adj = [unsafe_load(r_xadj[] ,i) for i=1:length(m.elements)+1]
    x_adjncy = [unsafe_load(r_adjncy[],i) for i=1:x_adj[end] ]
    return Graph(metis_fmt_to_vector(x_adj, x_adjncy, mini_node))
end 

"""
metis\\_mesh\\_to\\_dual(;ne::Int64, nn::Int64, eptr::Array{Int64,1}, eind::Array{Int64,1}, baseval::Int64, ncommon::Int64)
call the METIS_MeshToDual function
"""

function metis_mesh_to_dual(;ne::T, nn::T , eptr::Vector{T}, eind::Vector{T}, ncommon::T, baseval::T) where {T}
    if "METIS_LIB" in keys(ENV)
        metis_str = ENV["METIS_LIB"]
    else
        metis_str = "/usr/lib/libmetis.so"
    end
    lib_metis = dlopen(metis_str; throw_error=false)
    @debug lib_metis
    @assert lib_metis != nothing
    mesh_to_dual_ptr = dlsym(lib_metis, :METIS_MeshToDual)
    @debug "METIS_MeshToDual Pointer", mesh_to_dual_ptr
    r_xadj = Ref{Ptr{T}}()
    r_adjncy = Ref{Ptr{T}}()
    ccall(mesh_to_dual_ptr, Cvoid, 
          (Ref{T}, Ref{T}, Ptr{T}, Ptr{T}, Ref{T}, Ref{T}, Ref{Ptr{T}}, Ref{Ptr{T}}),
          Ref{T}(ne),
          Ref{T}(nn),
          eptr,
          eind,
          Ref{T}(ncommon),
          Ref{T}(baseval),
          r_xadj,
          r_adjncy
         )
    x_adj = GC.@preserve r_xadj [unsafe_load(r_xadj[] ,i) for i=1:ne+1]
    x_adjncy = GC.@preserve r_adjncy [unsafe_load(r_adjncy[],i) for i=1:x_adj[end] ]
    return x_adj, x_adjncy
end 


