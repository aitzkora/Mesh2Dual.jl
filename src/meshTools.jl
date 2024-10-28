"""
shift_par_msh!(;elmdist::Vector{T}, eptr::Vector{T}, eind::Vector{T}, shift::T) 
shifts parallely mesh in place by a constant 
"""
function shift_par_msh!(;elmdist::Vector{T}, eptr::Vector{T}, eind::Vector{T}, shift::T) where {T}
  elmdist .+= shift
  eptr .+= shift
  eind .+= shift
end

"""
shift_msh!(;elmdist, eptr, eind, shift)
shifts sequentially a parallel mesh in place by a constant 
"""
function shift_msh!(; eptr::Vector{Vector{T}}, eind::Vector{Vector{T}}, elmdist::Vector{T}, shift::T) where {T}
  size = length(elmdist) - 1
  elmdist .+= shift
  for rank=0:size-1
    eptr[rank+1] .+= shift
    eind[rank+1] .+= shift
  end
end

function split_msh(eptr::Vector{T}, eind::Vector{T}, p::Int64) where {T}
  IndexT = typeof(length(T[0]))
  elmdist = Vector{T}(undef, p+1)
  eptr_p = [T[] for _=1:p]
  eind_p = [T[] for _=1:p] 
  elmmax = length(eptr) - 1

  startIndex = datascan(elmmax, IndexT(0), IndexT(p))
  endIndex = startIndex + datasize(elmmax, IndexT(0), IndexT(p)) - 1

  elmdist[1] = startIndex
  elmdist[2] = endIndex
  eptr_p[1] = eptr[startIndex+1:endIndex+1]
  eind_p[1] = eind[eptr[startIndex+1]+1:eptr[endIndex+2]]
  for i=2:p
    startIndex = datascan(elmmax, IndexT(i-1) , IndexT(p))
    endIndex = startIndex + datasize(elmmax, IndexT(i-1), IndexT(p)) -1
    elmdist[i+1] = elmdist[i] + endIndex-startIndex + 1
    eptr_p[i] = eptr[startIndex+1:endIndex]
    eind_p[i] = eind[eptr[startIndex+1]+1:eptr[endIndex+2]]
  end
  return eptr_p, eind_p, elmdist
end

function split_msh2(eptr::Vector{T}, eind::Vector{T}, p::Int64) where {T}
  adj = csr_to_list(eptr, eind)
  IndexT = typeof(length(T[0]))
  elmdist =Vector{T}(undef, p+1)
  eptr_p = [T[] for _=1:p]
  eind_p = [T[] for _=1:p] 
  n = length(adj)

  startIndex = datascan(n, IndexT(0), IndexT(p))
  nChunk = datasize(n, IndexT(0), IndexT(p))
  adj_loc = adj[startIndex+1:startIndex + nChunk]
  eptr_p[1], eind_p[1] = list_to_csr(adj_loc)
  elmdist[1] = T(0)
  elmdist[2] = T(0) + nChunk
  for i=2:p
    startIndex = datascan(n,  IndexT(i-1) , IndexT(p))
    nChunk = datasize(n, IndexT(i-1), IndexT(p))  
    elmdist[i+1] = elmdist[i] + nChunk
    adj_loc = adj[startIndex+1:startIndex + nChunk]
    eptr_p[i], eind_p[i] = list_to_csr(adj_loc)
  end
  return eptr_p, eind_p, elmdist
end
