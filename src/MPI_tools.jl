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
print_master(str)
print only on master process
"""

function print_master(str)
   _, _, rank = get_com_size_rank()
   if (rank == 0)
       println(str)
   end
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
tile(nMax, r, nMini, p)

computes bounds of a r-th part of the global partition of set \\{nMini...nMax\\} in p parts
"""

function tile(nn::T, rank::T, baseval::T, p::T) where {T<:Integer}
  startIndex = datascan(nn, rank, p) + baseval
  chunkSize = datasize(nn, rank, p)
  endIndex = startIndex + chunkSize - 1
  return startIndex, endIndex, chunkSize
end

"""
toProc(nMax, k, nMini, p) 

maps a k index in [nMini,nMax] to a proc number belonging to [0, p]
according to the formula taken from "Parallel programming with Coarrays", p. 12

""" 
function toProc(vertglbnbr::T, k::T, baseval::T, p::T) where {T<:Integer}
  k -= baseval
  blokrmnval = 1 + (vertglbnbr - 1) % p
  bloksizmax = ((vertglbnbr - 1) + p ) ÷ p
  bloksizmin = bloksizmax - 1
  proclocnum = k ÷ bloksizmax
  if ((proclocnum >= blokrmnval) && (bloksizmin > 0))
     proclocnum = blokrmnval + (k - blokrmnval * bloksizmax) ÷ bloksizmin
  end 
  return proclocnum
end

"""
taken from scotch
give
"""
function datascan(n, i, p) 
  return i * (n ÷ p) + ((i > (n % p)) ? (n % p) : i)
end

function datasize(n, i, p)
  return  (n + (p - 1 - i)) ÷ p
end
