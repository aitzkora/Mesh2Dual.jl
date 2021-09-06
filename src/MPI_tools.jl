function send_lists(vec::Vector{Vector{T}}) where T
  send_counts = T[length(i) for i in vec]
  comm = MPI.COMM_WORLD
  recv_counts = fill(T(0), MPI.Comm_size(comm))
  MPI.Alltoall!(UBuffer(send_counts, 1), UBuffer(recv_counts, 1), comm)
  recv = fill(T(0), sum(recv_counts))
  send_vec = collect(Iterators.flatten([v[:] for v in vec]))
  MPI.Alltoallv!(VBuffer(send_vec, send_counts), VBuffer(recv, recv_counts), comm)
  displs = cumsum(recv_counts)
  recv_vec = pushfirst!([recv[range(i...)] for i in zip(displs[1:end-1].+1,displs[2:end])], recv[1:recv_counts[1]])
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
