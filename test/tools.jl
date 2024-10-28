using MPI
using Mesh2Dual
using Test

MPI.Init()
comm, size, rank = get_com_size_rank()

# tile
n = 8
mCheck = floor(n/size)
s,e,m = tile(n, rank, 1, Int64(size))
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
for k=0:7
  @test toProc(nMax, k, 0, size) == toProcCheck[k+1]
end

# send_lists

@assert size == 3
if rank == 0
  verttab = [0, 0, 2, 5]
  edgetab = [1, 2, 3, 4, 5]
elseif rank == 1
  verttab = [0, 1, 1, 4] 
  edgetab = [6, 7, 8, 9]
elseif rank == 2
  verttab = [0, 1, 3, 3]
  edgetab = [10, 11, 12]
end

# do the transpose
verttab_t, edgetab_t = send_lists(verttab, edgetab)

if rank == 0
  verttab_c = [0, 0, 1, 2]
  edgetab_c = [6, 10]
elseif rank == 1
  verttab_c = [0, 2, 2, 4] 
  edgetab_c = [1, 2, 11,12]
elseif rank == 2
  verttab_c = [0, 3, 6, 6]
  edgetab_c = [3, 4, 5, 7, 8, 9]
end
@test verttab_c == verttab_t
@test edgetab_c == edgetab_t

MPI.Finalize()
@test MPI.Finalized()
MPI.Finalize()
@test MPI.Finalized()
