using Mesh2Dual
using MPI
using Test

MPI.Init()
comm = MPI.COMM_WORLD
size = MPI.Comm_size(comm)
rank = MPI.Comm_rank(comm)
root = 0
# taken from SO (https://stackoverflow.com/questions/27929766/how-to-print-text-with-a-specified-rgb-value-in-julia)
function print_rgb(r, g, b, t)
  print("\e[1m\e[38;2;$r;$g;$b;249m",t)
end

if (rank == 0)
  print_rgb(200, 0, 200, "loading mesh\n")
end
if (length(ARGS) >= 2)
  filename = ARGS[1]
  dim = parse(Int64, ARGS[2])
else
  filename ="cube.1.mesh" 
  dim = 3
end
if (dim == 2)
  eptr, eind, elmdist  = read_par_mesh(filename, Val(2), Int32(0)) 
else
  eptr, eind, elmdist  = read_par_mesh(filename, Val(3), Int32(0)) 
end
println("I'm $rank, nodes ∈ [$(minimum(eind)), $(maximum(eind))]")
println("I'm $rank, elements ∈ [$(elmdist[rank+1]), $(elmdist[rank+2])]")
# call parmetis
if (rank == 0)
  print_rgb(200, 0, 200, "computing dual mesh by parmetis\n")
end
xadj, adjcy = parmetis_mesh_to_dual(;elmdist=elmdist, eptr=eptr, eind=eind, baseval=Int32(0), ncommon=Int32(2), 
                                    comm = comm)
adj_check = metis_fmt_to_vector(xadj, adjcy, Int32(0))
# our algorithm
if (rank == 0)
  print_rgb(200,0, 200, "computing dual mesh by our algorithm\n")
end
adj = dgraph_dual(;elmdist=elmdist, eptr=eptr, eind=eind, baseval=Int32(0), ncommon=Int32(2), comm = comm) 
if (rank == 0)
  print_rgb(200, 0, 200, "checking")
end
@test length(adj) == length(adj_check)
for i=1:length(adj)
  @test sort(adj[i]) == sort(adj_check[i])
end
MPI.Finalize()
@test MPI.Finalized()
