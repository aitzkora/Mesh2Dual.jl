module Mesh2Dual

using MPI
using Revise

export Mesh, Graph, graph_dual, mesh_to_metis_fmt, metis_graph_dual, metis_fmt_to_graph, 
       mesh_to_scotch_fmt, graph_dual_new, metis_mesh_to_dual, SimplexMesh, parmetis_mesh_to_dual,
       dgraph_dual

"""
`Mesh` implements a very basic topological structure for meshes
"""
struct Mesh
    elements::Array{Array{Int64,1},1}
    nodes::Array{Int64,1}
    Mesh(elem::Array{Array{Int64,1},1})=
    begin
        new(elem, sort(reduce(union,elem)))
    end
end

"""
Mesh structure with simplicial elements
"""
struct SimplexMesh{D}
  dim::Int64
  nodes::Array{Float64,2}
  elements::Array{Int64,2}
  function SimplexMesh{D}(nodes::Array{Float64,2}, elements::Array{Int64,2}) where {D}
    @assert Val(D) isa Union{map(x->Val{x},1:3)...}
    @assert D == size(nodes, 2)
    @assert D+1 == size(elements, 2)
    new{D}(D, nodes, elements)
  end
end


"""
very basic Graph structure
"""
struct Graph
    adj::Array{Array{Int64,1},1}
end

"""
Overloaded equality operator for Graphs 
"""
Base.:(==)(a::Graph, b::Graph) = Base.:(==)(a.adj, b.adj)


using Test

"""
gen_parts(l::Array{Int64,1}, k::Int)

generates an Array contaning the (#l choose k) sets of k elements in the set l
"""
function gen_parts(l::Array{Int64,1}, k::Int64)
    n = size(l, 1)
    gen =[ (1:n) .*i for i in filter(x->sum(x)==k,digits.(0:2^n-1,base=2, pad=n))] 
    return map(x->Set(l[x]), setdiff.(unique.(gen),[0]))
end 

@testset "gen_parts_test" begin
   l = [1, 2, 4]
   @test Set(gen_parts(l, 2)) == Set(Set.([(1,2),(1,4),(2,4)]))
   @test Set(gen_parts(l, 1)) == Set(Set.([(1),(2),(4)])) 
   @test Set(gen_parts(l, 3)) == Set(Set.([(1,2,4)])) 
end

""" 
graph\\_dual(m::Mesh, n\\_common::Int64)

Compute the elements graph of the Mesh m respect to the following adjacency
relationship
e₁ ~ e₂ ↔ e₁,e₂ have n\\_common shared vertices 
"""
function graph_dual(m::Mesh, n_common::Int64)
    T = Dict()
    # build the hash table for each n_common selection
    for (i, e) in enumerate(m.elements)
         for r in gen_parts(e, n_common)
              if (! haskey(T, r) )
                  T[r] = Set([i])
              else
                  push!(T[r],i)
              end
         end
    end
    for ke in T
        @debug Tuple(ke.first), " → " ,typeof(ke.second)
    end
    #now build adjacency
    adj = [ Set{Int64}() for _ in 1:length(m.elements)]
    for (i, e) in enumerate(m.elements)
        for r in gen_parts(e, n_common)
            union!(adj[i], T[r])
        end
    end
    # suppress reflexive adjacency
    for (i, e) in enumerate(m.elements)
        setdiff!(adj[i], i)
    end
    return Graph(map(x->sort(collect(keys(x.dict))), adj))
end 

"""
mesh\\_to\\_metis\\_fmt(m::Mesh)
converts a Mesh struct to a adjacency list defining a Metis Mesh
i.e : nodes indexes must start to zero for metis
"""
function mesh_to_metis_fmt(m::Mesh)
    eptr = Cint[0]
    eind = Cint[] 
    mini_node = Int32(minimum(m.nodes))
    for (i, e) in enumerate(m.elements)
        append!(eptr, eptr[i]+length(e))
        append!(eind, e .- mini_node)
    end
    return (eptr, eind, mini_node)
end 

"""
mesh\\_to\\_scotch\\_fmt(m::Mesh)

converts a mesh to the SCOTCH format (i.e. a bipartite graph)
"""
function mesh_to_scotch_fmt(m::Mesh)
    eptr = Cint[0]
    eind = Cint[]
    baseval = minimum( m.nodes )
    T = Dict()
    # we construct an adjacency for nodes
    for (i, e) in enumerate(m.elements)
        append!(eptr, eptr[i]+length(e))
        append!(eind, e .- baseval)
        for n in e
            if (! haskey(T, n))
                T[n] = Cint[i]
            else
                union!(T[n], Cint[i])
            end
        end
    end
    ne = length(m.elements)
    for (i, n) in enumerate(m.nodes)
        append!(eptr, eptr[ne+i] + length(T[n]))
        append!(eind, T[n] .- baseval)
    end
    eind[1:eptr[ne+1]] .+= Cint[ne]
    return (eptr, eind , ne)
end 

"""
graph\\_dual\\_new(m::Mesh, n_common::Int = 1)
computes a dual graph(i.e. element graph) using a alternative method
prototype for the futur scotch function
"""
function graph_dual_new(m::Mesh, n_common::Int = 1)
  ptr, ind , ne = mesh_to_scotch_fmt(m)
  adj = [Array{Int64,1}() for _ in 1:ne]

  # first pass : compute adjacency for n_common = 1
  for i=1:ne # ∀ e ∈ elems(m)
    for n ∈ ind[1+ptr[i]:ptr[i+1]] # ∀ n ∈ e
      for e₂ ∈ ind[ptr[n+1]+1:ptr[n+2]] # ∀ e₂ ⊃ { n }
        if e₂ ≠ (i-1)
          union!(adj[i],e₂)
        end 
      end 
    end
  end

  if (n_common == 1)
    return adj
  else
    adjp = [Array{Int64,1}() for _ in 1:ne]
    for e₁=1:ne 
      for e₂ ∈ adj[e₁] 
        accu = 0
        for n₁ ∈ ind[1+ptr[e₁]:ptr[e₁+1]] 
          for n₂ ∈ ind[1+ptr[e₂+1]:ptr[e₂+2]] 
            if (n₁ == n₂ )
                accu+=1
            end 
          end
        end
        if (accu >= n_common)
            union!(adjp[e₁], e₂)
        end
      end 
    end
    return adjp
  end 

end

"""
metis\\_fmt\\_to\\_graph(eptr::Array{Cint,1}, eind::Array{Cint,1}, min_node::Cint)

converts a metis adjncy list to a graph
"""
function metis_fmt_to_graph(eptr::Array{Cint,1}, eind::Array{Cint,1}, min_node::Cint = Cint(1))
    elems = fill(Int64[],size(eptr,1)-1)
    nodes = Int64[]
    for i=1:length(eptr)-1
        elems[i] = eind[eptr[i]+1:eptr[i+1]] .+ min_node
    end
    return Graph(elems)
end

using Libdl: dlopen, dlsym

"""
metis\\_graph\\_dual(m::Mesh, n_common::Int)
computes the graph dual (i.e. elements graph using Metis)
"""

function metis_graph_dual(m::Mesh, n_common::Int)
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
    r_xadj = Ref{Ptr{Cint}}()
    r_adjncy = Ref{Ptr{Cint}}()
    ccall(grf_dual_ptr, Cvoid, 
          (Cint, Cint, Ptr{Cint}, Ptr{Cint}, Cint, Ref{Ptr{Cint}}, Ref{Ptr{Cint}}),
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
    return metis_fmt_to_graph(x_adj, x_adjncy, mini_node)
end 

"""
metis\\_mesh\\_to\\_dual(;ne::Int64, nn::Int64, eptr::Array{Int64,1}, eind::Array{Int64,1}, baseval::Int64, ncommon::Int64)
call the METIS_MeshToDual function
"""

function metis_mesh_to_dual(;ne::Int64, nn::Int64 , eptr::Array{T,1}, eind::Array{T,1}, ncommon::Int64, baseval::Int64) where {T}
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
    x_adj = [unsafe_load(r_xadj[] ,i) for i=1:ne+1]
    x_adjncy = [unsafe_load(r_adjncy[],i) for i=1:x_adj[end] ]
    return x_adj, x_adjncy
end 

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
    rank = MPI.Comm_rank(comm)
    ne_local = elmdist[rank+2]-elmdist[rank+1]
    println("je suis $rank, j'ai $ne_local elements")
    x_adj = GC.@preserve r_xadj [unsafe_load(r_xadj[] ,i) for i=1:ne_local+1]
    x_adjncy = GC.@preserve r_adjncy [unsafe_load(r_adjncy[],i) for i=1:x_adj[end] ]
    return x_adj, x_adjncy
end  # function parmetis_mesh_to_dual

"""
dgraph\\_dual\\(;elmdist, eptr, eind, baseval, ncommon, comm)
computes a distributed dual graph
"""
function dgraph_dual(;elmdist::Vector{T}, eptr::Vector{T}, eind::Vector{T}, baseval::T, ncommon::T, comm::MPI.Comm) where T
  r = MPI.Comm_rank(comm)
  p = MPI.Comm_size(comm)
  ne = elmdist[p+1]
  neLoc = elmdist[r+2]-elmdist[r+1]

  nnLoc = maximum(eind) - baseval + 1
  nMax = Ref{T}(nnLoc)
  n2e = [Vector{T}() for _ in 1:nnLoc]
  MPI.Allreduce!(nMax, MPI.MAX, comm)
  nn = nMax[]
  if (r == 0) 
    println("ne = $ne, nn = $nn")
  end
  # first pass : compute local n2e 
  nLocs = Set{T}()
  for i=1:neLoc
    for n ∈ eind[1+eptr[i]:eptr[i+1]] # ∀ n ∈ e
          union!(n2e[n+1],i)
          push!(nLocs, n)
    end 
  end 
  n2p = [ Vector{T}() for _ in 1:nn]
  for n=baseval:nn-1+baseval
    nInNlocs = (n in nLocs)
    countLoc = [(nInNlocs) ? T(1) : T(0)]
    counts = MPI.Allgather(countLoc, comm)
    nBuffLoc = (nInNlocs) ? [T(r)] : T[]
    n2p[n-baseval+1] = fill(zero(T), sum(counts))
    MPI.Allgatherv!(nBuffLoc, VBuffer(n2p[n-baseval+1], counts), comm)
  end
  if (r == 0)
    println("n2p = $n2p")
  end
end
end # module
