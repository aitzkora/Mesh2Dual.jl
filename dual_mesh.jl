using Mesh2Dual
using MPI
using Test

MPI.Init()
comm = MPI.COMM_WORLD
size = MPI.Comm_size(comm)
rank = MPI.Comm_rank(comm)
root = 0
println("loading mesh")
#eptr, eind, elmdist  = read_par_mesh("mesh_oliv.mesh", Val(3), Int32(0)) 
eptr, eind, elmdist  = read_par_mesh("cube.1.mesh", Val(3), Int32(0)) 
# call parmetis
println("computing dual mesh by parmetis")
xadj, adjcy = parmetis_mesh_to_dual(;elmdist=elmdist, eptr=eptr, eind=eind, baseval=Int32(0), ncommon=Int32(2), 
                                    comm = comm)
adj_check = metis_fmt_to_vector(xadj, adjcy, Int32(0))
MPI.Barrier(comm)
# our algorithm
println("computing dual mesh by our algorithm")
adj = dgraph_dual(;elmdist=elmdist, eptr=eptr, eind=eind, baseval=Int32(0), ncommon=Int32(2), comm = comm) 
MPI.Barrier(comm)
@test length(adj) == length(adj_check)
for i=1:length(adj)
  @test sort(adj[i]) == sort(adj_check[i])
end
MPI.Finalize()
@test MPI.Finalized()
