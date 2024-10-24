module Mesh2Dual

using MPI
using Revise
using Printf
using Libdl
export Mesh, Graph, graph_dual, mesh_to_metis_fmt, metis_graph_dual, metis_fmt_to_vector, Dgraph_header,
       mesh_to_scotch_fmt, graph_dual_new, metis_mesh_to_dual, SimplexMesh, parmetis_mesh_to_dual,
       dgraph_dual, gen_parts, read_par_mesh, toProc, tile, send_lists, ptscotchparmetis_mesh_to_dual,
       write_par_dmesh, write_par_dgraph, shift_par_msh!, list_to_csr, csr_to_list, read_dgraph_header, read_adj, parse_file_name,
       write_dmesh, shift_msh!, read_dmesh, read_scotch_mesh


include("misc.jl")
include("MPI_tools.jl")
include("graphIO.jl")
include("meshIO.jl")

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

# to do some comparisons with metis
include("metisCalls.jl")
include("ParMetisCall.jl")



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

"""
dgraph\\_dual\\(;elmdist, eptr, eind, baseval, ncommon, comm)
computes a distributed dual graph
"""
function dgraph_dual(;elmdist::Vector{T}, eptr::Vector{T}, eind::Vector{T}, baseval::T, ncommon::T, comm::MPI.Comm) where {T}
  t₀ = time()
  r = MPI.Comm_rank(comm)
  p = MPI.Comm_size(comm)
  ne = elmdist[p+1]
  neLoc = elmdist[r+2]-elmdist[r+1]

  nMax = Ref{T}(maximum(eind))
  MPI.Allreduce!(nMax, MPI.MAX, comm)
  nMin = Ref{T}(minimum(eind))
  MPI.Allreduce!(nMin, MPI.MIN, comm)
  nn = nMax[] - nMin[] + T(1)
  nMin = nMin[]
  nMax = nMax[]
  toSend = [ Vector{T}() for _ in 1:p]
  for i=1:neLoc # ∀ e local element
    for n ∈ eind[1+eptr[i]:eptr[i+1]]  # ∀ n ∈ e
      append!(toSend[toProc(nMax, n, nMin) + 1], i-1+elmdist[r+1], n)
    end 
  end 
  t0 = time()
  verttab_s, edgetab_s = list_to_csr(toSend)
  t0 = time()
  verttab, edgetab = send_lists(verttab_s, edgetab_s)
  #@printf "tps first All2All %.3e\n" time() - t0
  t0 = time()
  #@info "apres envoi", verttab, edgetab
  startIndex, endIndex , sizeChunk = tile(nMax, T(r), nMin)
  # n2p[idx] will contain the list of procs containing the node idx
  n2p = [ Vector{T}() for _ in startIndex:endIndex ]
  println("nodes range : [$startIndex, $endIndex], $sizeChunk")
  frsttab = fill(zero(T), 2 * sizeChunk)
  frsttab[2:2:2*sizeChunk] .= -one(T)
  # frsttab [2i] will have the value :
  #   - -1 if node i does not appear
  #   - the adress j in edgetab where to start the linked list of elements containing i
  # frstab [2i+1] will store  how many elements contain i
  # egdetab and frsttab will implement a linked list following the scheme
  # frsttab : .... | n_j ,  j |  
  #                | 2k+1,  2k|
  #                         |
  #                  ┌──────┘ 
  #                  | 
  #                  v
  # edgetab : ...  | e_j,  z |       | e_z, -1 | # -1 indicates the end of the linking process
  #                | j  , j+1|       | z  , z+1| 
  #                        |           ∧  
  #                        └───────────┘
  for proc=1:p
    @inbounds for j=1+verttab[proc]:2:verttab[proc+1]
       node = edgetab[j+1] # node global number included in proc
       nodex = node-startIndex+1 # local numerotation
       union!(n2p[nodex], proc-1) # acts like an hashtable
       frsttab[2*(nodex-1)+1] += one(T) 
       if (frsttab[2*nodex] == -one(T))
         edgetab[j+1] = -one(T)
       else
         edgetab[j+1] = frsttab[2*nodex]
       end
       frsttab[2*nodex] = j # pointer on edgetab
    end 
  end
  totSize = sum(frsttab[1:2:2*sizeChunk]) # sum the effectives of each node
  n2eptr = T[zero(T)]
  n2etab = fill(zero(T), totSize)
  # we will build a csr list (n2eptr,n2etab) which contains list of elements containing nodes
  curr = 0 
  for j=1:sizeChunk
    w = frsttab[2*j]
    while w != -one(T)
      n2etab[curr+1] = edgetab[w]
      curr += 1
      w = edgetab[w+1]
    end
    push!(n2eptr, curr) # neptr contains the the index of the csr list in n2etab
  end
  toSend = [ Vector{T}() for _ in 1:p]
  map(x->sizehint!(x,totSize + 2*sizeChunk), toSend) # not give a significant reduction of CPU time
  for i=1:sizeChunk
    for proc in n2p[i]
      push!(toSend[proc+1],i+startIndex-1, n2eptr[i+1]-n2eptr[i], n2etab[1+n2eptr[i]:n2eptr[i+1]]...)
    end
  end
  verttab, edgetab = list_to_csr(toSend)
  verttab_t, edgetab_t = send_lists(verttab, edgetab)
  # build nn2
  # nn2e is an hashmap containing the list of edges which contains n
  nn2e = Dict{T,Vector{T}}()
  for proc=1:p
    curr = 1+verttab_t[proc]
    while curr < verttab_t[proc+1]
      currNode = edgetab_t[curr]
      nNodes = edgetab_t[curr+1]
      if !haskey(nn2e, currNode) 
        nn2e[currNode] = edgetab_t[curr+2:curr+1+nNodes]
      end 
      curr += 2 + nNodes 
    end
  end
  # build 1-adjacency
  adj = [Set{T}() for _ in 1:neLoc]
  for i=1:neLoc 
    e = i-1+elmdist[r+1]
    for n ∈ eind[1+eptr[i]:eptr[i+1]]  
      for e₂ ∈ nn2e[n]
        if (e₂ ≠ e)
          push!(adj[i], e₂)
        end  
      end 
    end 
  end
  adj = map(x->[ i for i in x], adj)
  if (ncommon > 1) 
    t0 = time()
    adjp = [Set{T}() for _ in 1:neLoc]
    for i=1:neLoc
      accu = fill(T(0), length(adj[i]))
      for n ∈ eind[1+eptr[i]:eptr[i+1]] 
        for (j,e₂) ∈ enumerate(adj[i])
          if ((e₂ ∈ nn2e[n]) & (e₂ ≠ ((i-1) + elmdist[r+1] )))
            accu[j] += 1
          end 
        end
      end 
      for (j,e₂) ∈ enumerate(adj[i])
        if (accu[j] >= ncommon)
          push!(adjp[i], e₂)
        end
      end
    end 
    adjp = map(x->[i for i in x], adjp)
  else
    adjp = adj
  end 
  return adjp
end 
end # module
