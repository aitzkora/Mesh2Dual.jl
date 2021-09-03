using MPI
using Test

MPI.Init()
comm = MPI.COMM_WORLD
size = MPI.Comm_size(comm)
rank = MPI.Comm_rank(comm)

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
@assert size == 3
if rank == 0
  send = [Int32[] , Int32[1, 2], Int32[ 3, 4, 5]]
elseif rank == 1
  send = [Int32[6], Int32[], Int32[7, 8, 9]]
elseif rank == 2
  send = [Int32[10], Int32[11, 12], Int32[] ]
end
recv = send_lists(send)
println("I'm $rank, recv = $recv")
MPI.Finalize()
@test MPI.Finalized()
