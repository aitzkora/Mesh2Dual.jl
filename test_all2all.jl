using MPI
using Test

MPI.Init()
comm = MPI.COMM_WORLD
size = MPI.Comm_size(comm)
rank = MPI.Comm_rank(comm)
root = 0

@assert size == 3
if rank == 0
  send = Int32[1, 2, 3, 4, 5]
  recv = fill(Int32(0), 3)
  send_counts = Int32[0, 2, 3]
  recv_counts = Int32[0, 1, 2]
elseif rank == 1
  send = Int32[6, 7, 8, 9]
  recv = fill(Int32(0), 4)
  send_counts = Int32[1, 0, 3]
  recv_counts = Int32[2, 0, 2]
elseif rank == 2
  send = Int32[10, 11, 12]
  recv = fill(Int32(0), 6)
  send_counts = Int32[1, 2, 0]
  recv_counts = Int32[3, 3, 0]
end
MPI.Alltoallv!(VBuffer(send,send_counts), VBuffer(recv,recv_counts), comm)
println("I'm $rank, recv = $recv")
MPI.Barrier(comm)
#end
MPI.Finalize()
@test MPI.Finalized()
