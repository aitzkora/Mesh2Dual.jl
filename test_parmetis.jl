using Mesh2Dual
using MPI
using Test

MPI.Init()
comm = MPI.COMM_WORLD
size = MPI.Comm_size(comm)
rank = MPI.Comm_rank(comm)

@assert size == 3
eptr = Int32[]
ind = Int32[]
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


# our algorithm
xadj_m, adjcy_m = parmetis_mesh_to_dual(;elmdist=elmdist, eptr=eptr, eind=eind, baseval=Int32(0), ncommon=Int32(2), comm = comm)
adj_check_m = metis_fmt_to_vector(xadj_m, adjcy_m, Int32(0))
MPI.Barrier(comm)
adj = dgraph_dual(;elmdist=elmdist, eptr=eptr, eind=eind, baseval=Int32(0), ncommon=Int32(2), comm = comm) 
MPI.Barrier(comm)
# call parmetisj
xadj, adjcy = ptscotchparmetis_mesh_to_dual(;elmdist=elmdist, eptr=eptr, eind=eind, baseval=Int32(0), ncommon=Int32(2), 
                                    comm = comm)
adj_check = metis_fmt_to_vector(xadj, adjcy, Int32(0))
MPI.Barrier(comm)
@test length(adj) == length(adj_check)
for i=1:length(adj)
  @test sort(adj[i]) == sort(adj_check[i])
end
for i=1:length(adj)
  @test sort(adj[i]) == sort(adj_check_m[i])
end
@info "je suis", rank, xadj, adjcy
sleep(100)
MPI.Finalize()
@test MPI.Finalized()
