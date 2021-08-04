using Mesh2Dual
using MPI


MPI.Init()
comm = MPI.COMM_WORLD
size = MPI.Comm_size(comm)
rank = MPI.Comm_rank(comm)
root = 0

@assert size == 3
eptr = Int32[]
eind = Int32[]
if rank == 0
  eptr = Int32[ 0, 4, 8 ]
  eind = Int32[ 0, 1, 5, 6, 1, 2, 6, 7 ]
elseif rank == 1
  eptr = Int32[ 0, 4, 8 ]
  eind = Int32[ 2, 3, 7, 8, 3, 4, 8, 9 ]
elseif rank == 2
  eptr = Int32[0, 4, 8, 12, 16]
  eind = Int32[5, 6, 10, 11, 6, 7, 11, 12, 7, 8, 12, 13, 8, 9, 13, 14]
end
elmdist = Int32[0, 2, 4, 8]
xadj, adjcy = parmetis_mesh_to_dual(;elmdist=elmdist, eptr=eptr, eind=eind, baseval=Int32(0), ncommon=Int32(2), 
                                    comm = comm)
MPI.Barrier(comm)
println("I'm $rank, xadj = $xadj, adjcy = $adjcy")
MPI.Finalize()
