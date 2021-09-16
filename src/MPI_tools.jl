"""
listToCsr(vec::Vector{Vector{T}}) where {T:<Int}

converts a list of list to CSR format
"""

function listToCsr(vec::Vector{Vector{T}}) where {T}
  eptr = T[0]
  eind = T[] 
  for (i, e) in enumerate(vec)
      append!(eptr, eptr[i]+length(e))
      append!(eind, e)
  end
  return (eptr, eind)
end 

"""
function send_lists(verttab::Vector{T}, edgetab::Vector{T}) where T

do a MPI_AlltoAllv with a CSR graph

"""
function send_lists(verttab::Vector{T}, edgetab::Vector{T}) where T
  comm = MPI.COMM_WORLD
  p = MPI.Comm_size(comm)
  send_counts = diff(verttab)
  @assert size(send_counts, 1) == p
  recv_counts = fill(T(0), p)
  MPI.Alltoall!(UBuffer(send_counts, 1), UBuffer(recv_counts, 1), comm)
  recv = fill(T(0), sum(recv_counts))
  MPI.Alltoallv!(VBuffer(edgetab, send_counts), VBuffer(recv, recv_counts), comm)
  verttab_t = pushfirst!(accumulate(+, recv_counts), T(0))
  edgetab_t = recv
  return verttab_t, edgetab_t
end

"""
tile(nMini::T, rank::T; nMini)

computes a partition of \\{nMini...nMax\\}, into p MPI processes
"""

function tile(nMax::T, rank::T, nMini::T = T(1)) where{T<:Integer}
  # tiling indexes 
  p = MPI.Comm_size(MPI.COMM_WORLD)
  n = nMax - nMini + 1
  m = convert(T, ceil(n / p))
  r = convert(T, n % p)
  if r == T(0) 
    startIndex = rank * m +  nMini
    return startIndex, startIndex + m - T(1), m
  end
  if rank + T(1) > r
    m -= 1 
    startIndex = rank * m + r + nMini
  else
    startIndex = rank * m + nMini
  end
  endIndex = startIndex + m - T(1)
  return startIndex, endIndex, m
end

"""
toProc(n::T, k::T) where {T<:Integer}

maps a k index to a proc number belonging to [0, MPI.Comm_size]
according to the formula taken from "Parallel programming with Coarrays", p. 12

""" 
function toProc(nMax::T, k::T, nMini::T = T(1)) where {T<:Integer}
  p = MPI.Comm_size(MPI.COMM_WORLD)
  n = nMax - nMini + 1
  m = convert(T, ceil(n / p))
  r = n % convert(T, p)
  if r == 0
    return convert(T, floor((k - nMini) / m))
  end
  if k ≤ m * r
    return convert(T, floor((k - nMini) / m))
  else
    return convert(T, floor((k - r - nMini) / (m-1)))
  end
end
