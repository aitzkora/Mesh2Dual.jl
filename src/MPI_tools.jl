"""
function get_com_size_rank()
returns the triplet comm, size and rank of MPI : syntastic sugar
"""
function get_com_size_rank()
  comm = MPI.COMM_WORLD
  size = MPI.Comm_size(comm)
  rank = MPI.Comm_rank(comm)
  return comm, size, rank
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
  return verttab_t, recv
end

"""
tile(nMini::T, rank::T; nMini)

computes a partition of \\{nMini...nMax\\}, into p MPI processes
rank is in [0, nbproc-1]
"""

function tile(nMax::T, rank::T, nMini::T = T(0)) where{T<:Integer}
  # tiling indexes 
  p = convert(T, MPI.Comm_size(MPI.COMM_WORLD))
  n = nMax - nMini + T(1)
  m = ceil(T, n / p) # could be replace by (n-1)/p + 1
  r = n % p
  if (r > T(0) && rank >= r)
    m = m - 1
    startIndex = rank * m + r + nMini
  else
    startIndex = rank * m + nMini
  end
  endIndex = startIndex + m - T(1)
  return startIndex, endIndex, m
end

"""
toProc(n::T, k::T, nMini::T = T(0)) where {T<:Integer}

maps a k index in [nMini,nMax] to a proc number belonging to [0, MPI.Comm_size]
according to the formula taken from "Parallel programming with Coarrays", p. 12

""" 
function toProc(nMax::T, k::T, nMini::T = T(0)) where {T<:Integer}
  p = MPI.Comm_size(MPI.COMM_WORLD)
  n = nMax - nMini + 1
  m = convert(T, ceil(n / p))
  r = n % convert(T, p)
  if (r == 0  || (k-nMini +1) <= m * r)
    return floor(T, (k - nMini) / m)
  else
    return floor(T, (k - r - nMini) / (m-1))
  end 
end
