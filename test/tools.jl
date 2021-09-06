using MPI
using Mesh2Dual
using Test

MPI.Init()
comm = MPI.COMM_WORLD
size = MPI.Comm_size(comm)
rank = MPI.Comm_rank(comm)

# tile
n = 8
mCheck = floor(n/size)
s,e,m = tile(n, rank, 1)
if rank == 0
  sCheck, mCheck = 1, 3
elseif rank == 1
  sCheck, mCheck = 4, 3
else rank == 2
  sCheck, mCheck = 7, 2
end
eCheck = sCheck + mCheck - 1
@test s == sCheck
@test m == mCheck
@test e == eCheck

# toProc
nMax = 8
toProcCheck=[0,0,0,1,1,1,2,2]
for k=1:8
  #println("$k -> $(toProc(nMax,k))")
  @test toProc(nMax,k) == toProcCheck[k]
end
MPI.Finalize()
@test MPI.Finalized()
