module Mesh2Dual

using MPI
using Revise

export Mesh, Graph, graph_dual, mesh_to_metis_fmt, metis_graph_dual, metis_fmt_to_vector, 
       mesh_to_scotch_fmt, graph_dual_new, metis_mesh_to_dual, SimplexMesh, parmetis_mesh_to_dual,
       dgraph_dual, gen_parts, read_par_mesh

"""
`Mesh` implements a very basic topological structure for meshes
"""
struct Mesh{T}
    elements::Vector{Vector{T}}
    nodes::Vector{T}
    function Mesh{T}(elem::Vector{Vector{T}}) where {T}
      new{T}(elem, sort(reduce(union,elem)))
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
struct Graph{T}
    adj::Vector{Vector{T}}
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
function gen_parts(l::Vector{T}, k::T) where {T}
    n = size(l, 1)
    gen =[ (T(1):T(n)) .*i for i in filter(x->sum(x)==k,digits.(T(0):T(2^n-1),base=T(2), pad=T(n)))] 
    return map(x->Set(l[x]), setdiff.(unique.(gen),[T(0)]))
end 

""" 
graph\\_dual(m::Mesh{T}, n\\_common::T)

Compute the elements graph of the Mesh m respect to the following adjacency
relationship
e₁ ~ e₂ ↔ e₁,e₂ have n\\_common shared vertices 
"""
function graph_dual(m::Mesh{T}, n_common::Int) where {T}
  n_common = T(n_common)
  table = Dict()
  # build the hash table for each n_common selection
  for (i, e) in enumerate(m.elements)
    for r in gen_parts(e, n_common)
      if (! haskey(table, r) )
        table[r] = Set([T(i)])
      else
        push!(table[r], T(i))
      end
    end
  end
  for ke in table
      @debug Tuple(ke.first), " → " ,typeof(ke.second)
  end
  #now build adjacency
  adj = [ Set{T}() for _ in 1:length(m.elements)]
  for (i, e) in enumerate(m.elements)
    for r in gen_parts(e, n_common)
      union!(adj[i], table[r])
    end
  end
  # suppress reflexive adjacency
  for (i, e) in enumerate(m.elements)
    setdiff!(adj[i], T(i))
  end
  return Graph(map(x->sort(collect(keys(x.dict))), adj))
end 

"""
mesh\\_to\\_metis\\_fmt(m::Mesh)
converts a Mesh struct to a adjacency list defining a Metis Mesh
i.e : nodes indexes must start to zero for metis
"""
function mesh_to_metis_fmt(m::Mesh{T}) where {T}
    eptr = T[0]
    eind = T[] 
    mini_node = T(minimum(m.nodes))
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
function mesh_to_scotch_fmt(m::Mesh{T}) where {T}
    eptr = T[0]
    eind = T[]
    baseval = minimum( m.nodes )
    table = Dict()
    # we construct an adjacency for nodes
    for (i, e) in enumerate(m.elements)
        append!(eptr, eptr[i]+length(e))
        append!(eind, e .- baseval)
        for n in e
            if (! haskey(table, n))
                table[n] = T[i]
            else
                union!(table[n], T[i])
            end
        end
    end
    ne = length(m.elements)
    for (i, n) in enumerate(m.nodes)
        append!(eptr, eptr[ne+i] + length(table[n]))
        append!(eind, table[n] .- baseval)
    end
    eind[1:eptr[ne+1]] .+= T[ne]
    return (eptr, eind , ne)
end 

"""
graph\\_dual\\_new(m::Mesh, n_common::Int = 1)
computes a dual graph(i.e. element graph) using a alternative method
prototype for the futur scotch function
"""
function graph_dual_new(m::Mesh{T}, n_common::Int = Int(1)) where {T}
  n_common = T(n_common)
  ptr, ind , ne = mesh_to_scotch_fmt(m)
  adj = [Vector{T}() for _ in 1:ne]

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
    adjp = [Vector{T}() for _ in 1:ne]
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
metis\\_fmt\\_to\\_graph(eptr::Vector{T}, eind::Vector{T}, min_node::T)

converts a metis adjacency list to a Vector of Vector
"""
function metis_fmt_to_vector(eptr::Vector{T}, eind::Vector{T}, min_node::T = T(1)) where {T}
    elems = fill(T[],size(eptr,1)-1)
    nodes = T[]
    for i=1:length(eptr)-1
        elems[i] = eind[eptr[i]+1:eptr[i+1]] .+ min_node
    end
    return elems
end

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
    x_adj = GC.@preserve r_xadj [unsafe_load(r_xadj[] ,i) for i=1:ne_local+1]
    x_adjncy = GC.@preserve r_adjncy [unsafe_load(r_adjncy[],i) for i=1:x_adj[end] ]
    return x_adj, x_adjncy
end  # function parmetis_mesh_to_dual

"""
dgraph\\_dual\\(;elmdist, eptr, eind, baseval, ncommon, comm)
computes a distributed dual graph
"""
function dgraph_dual(;elmdist::Vector{T}, eptr::Vector{T}, eind::Vector{T}, baseval::T, ncommon::T, comm::MPI.Comm) where {T}
  r = MPI.Comm_rank(comm)
  p = MPI.Comm_size(comm)
  ne = elmdist[p+1]
  neLoc = elmdist[r+2]-elmdist[r+1]

  nnLoc = maximum(eind) - baseval + 1
  nMax = Ref{T}(nnLoc)
  n2e = [Vector{T}() for _ in 1:nnLoc]
  MPI.Allreduce!(nMax, MPI.MAX, comm)
  nn = nMax[]

  # first pass : compute local n2e 
  nLocs = Set{T}()
  for i=1:neLoc # ∀ e local element
    for n ∈ eind[1+eptr[i]:eptr[i+1]]  #∀ n ∈ e
      union!(n2e[n-baseval+1],i-1+elmdist[r+1])
      push!(nLocs, n)
    end 
  end 
  n2p = [ Vector{T}() for _ in 1:nn]
  for n=baseval:nn-1+baseval # ∀ n 
    nInNlocs = (n in nLocs)
    countLoc = [(nInNlocs) ? T(1) : T(0)]
    counts = MPI.Allgather(countLoc, comm)
    nBuffLoc = (nInNlocs) ? [T(r)] : T[]
    n2p[n-baseval+1] = fill(zero(T), sum(counts))
    MPI.Allgatherv!(nBuffLoc, VBuffer(n2p[n-baseval+1], counts), comm)
  end
  # compute the sizes of exchanges
  max_size = maximum(map(length, n2p))
  size_swap = Dict( i => Dict{T,Vector{T}}() for i in nLocs if length(n2p[i-baseval+1]) > 1 )
  for n=nLocs
    if length(n2p[n-baseval+1])> 1
      for r2=n2p[n-baseval+1]
        if (r2 != r) 
          tag1 = n * p * p + r * p + r2 
          tag2 = n * p * p + r2 * p + r 
          size_swap[n][r] = [T(length(n2e[n-baseval+1]))]
          size_swap[n][r2] = [T(0)]
          MPI.Sendrecv!(@view(size_swap[n][r][1]), r2, tag2, @view(size_swap[n][r2][1]), r2, tag1, comm)
        end
      end
    end
  end

  # do the exchange
  tab_swap = Dict( i => Dict{T,Vector{T}}() for i in nLocs if length(n2p[i-baseval+1]) > 1)
  for n=nLocs
    if length(n2p[n-baseval+1])> 1
      for r2=n2p[n-baseval+1]
        if (r2 != r) 
          tab_swap[n][r] = n2e[n-baseval+1]
          tab_swap[n][r2] = fill(T(0), size_swap[n][r2][1])
          tag1 = n * p * p + r * p + r2 
          tag2 = n * p * p + r2 * p + r 
          MPI.Sendrecv!(@view(tab_swap[n][r][:]), r2, tag2, @view(tab_swap[n][r2][:]), r2, tag1, comm)
        end
      end
    end
  end

  # merge local node to elements with exchanged data
  for n=nLocs
    if length(n2p[n-baseval+1])> 1
      for r2=n2p[n-baseval+1]
        if (r2 != r) 
          union!(n2e[n-baseval+1], tab_swap[n][r2])
        end
      end
    end
  end

  # Now we have all the information to build 1-dual graph
  adj = [Array{T,1}() for _ in 1:neLoc]
  for i=1:neLoc # ∀ e local element
    for n ∈ eind[1+eptr[i]:eptr[i+1]]  #∀ n ∈ e
      for e₂ ∈ n2e[n-baseval+1] #∀ e₂ ⊃ { n }
        if e₂ ≠ ((i-1) + elmdist[r+1])
          union!(adj[i],e₂) #! beware for baseval of i
       end
     end
    end 
  end 
  if (ncommon > 1) 
    adjp = [Vector{T}() for _ in 1:neLoc]
    for i=1:neLoc
      accu = fill(T(0), length(adj[i]))
      for n ∈ eind[1+eptr[i]:eptr[i+1]] 
        for (j,e₂) ∈ enumerate(adj[i])
          if ( (e₂ ∈ n2e[n-baseval+1]) & (e₂ ≠ ((i-1) + elmdist[r+1] )))
              accu[j] += 1
          end 
        end
      end 
      for (j,e₂) ∈ enumerate(adj[i])
        if (accu[j] >= ncommon)
          union!(adjp[i], e₂)
        end
      end
    end 
  else
    adjp = adj
  end 
  return adjp
end # function

function tile(size_obj, rank) # TODO : parametrize it! 
  # tiling indexes 
  p = MPI.Comm_size(MPI.COMM_WORLD)
  size_chunk = (size_obj - 1) ÷ p + 1
  mod_chunk = size_obj % p 
  if (mod_chunk != 0 & (rank + 1 > mod_chunk))
    size_chunk = size_chunk - 1
    startIndex = rank * size_chunk + mod_chunk + 1
  else
    startIndex = rank * size_chunk + 1
  end
  endIndex = startIndex + size_chunk -1
  return startIndex, endIndex, size_chunk
end


function read_par_mesh(filename, ::Val{D}, baseval) where {D}
  IndexType = typeof(baseval)
  @assert Val(D) isa Union{map(x->Val{x},1:3)...}
  comm = MPI.COMM_WORLD
  p = MPI.Comm_size(comm)
  r = MPI.Comm_rank(comm)
  isMaster = (r == 0)
  if (isMaster)
    lines = readlines(filename)
    lines = lines[map(x-> length(x) > 0 && (x[1] != '#'), lines)] # beware to the '
    f_s(x) = findall(map(y->occursin(x,y),lines))[1]
    dim = parse(IndexType,lines[f_s("Dimension")+1])
    @debug "dim =", dim
    @assert  D == dim
    of_nodes = f_s("Vertices")
    nb_nodes = parse(IndexType, lines[of_nodes + 1])
    nodes = zeros(Float64, nb_nodes, dim)
    for i=1:nb_nodes
      nodes[i, : ] = map(x->parse(Float64,x), split(lines[of_nodes + 1 + i])[1:end-1])
    end
    @debug "nodes=", nodes
    if (dim == 2)
      of_elem = f_s("Triangles")
    else
      of_elem = f_s("Tetrahedra")
    end
    nb_elem = parse(IndexType, lines[of_elem + 1])
  else
    nb_elem = IndexType(0)
    dim = IndexType(0)
  end
  nb_elem = MPI.bcast(nb_elem, 0, comm)
  dim = MPI.bcast(dim, 0, comm)
  elmdist  = zeros(IndexType, p+1)

  if (isMaster) 
    elements = zeros(IndexType, nb_elem, dim + 1)
    for i=1:nb_elem
      elements[i, : ] = map(x->parse(IndexType,x), split(lines[of_elem + 1 + i])[1:dim+1])
    end
    @debug "$elements"
    startIndex, endIndex, _ = tile(nb_elem, 0)
    eind = elements[startIndex:endIndex,:]
    elmdist[1] = startIndex - 1
    elmdist[2] = endIndex
    for i=2:p
      startIndex, endIndex, _ = tile(nb_elem, i-1)
      elmdist[i+1] = elmdist[i] + endIndex-startIndex + 1
      eind_send = elements[startIndex:endIndex,:]
      MPI.send(eind_send, i - 1, i , comm)
    end 
  end
  elmdist = MPI.bcast(elmdist,0, comm)
  if (!isMaster)
    (eind, _) = MPI.recv(0, r + 1, comm)
  end  

  size_chunk = (nb_elem - 1) ÷ p + 1
  eptr = [IndexType(0):IndexType(dim+1):IndexType(size_chunk * (dim + 1));]
  eind = eind'[:]
  @debug "$r -> $eptr, $eind, $elmdist"
  return (eptr, eind, elmdist)
end
end # module

