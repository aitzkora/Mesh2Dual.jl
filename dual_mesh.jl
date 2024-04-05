using Mesh2Dual
using MPI
using Test
using Printf
using DelimitedFiles

MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
root = 0
if (rank == 0)
  println("loading mesh")
end
if (length(ARGS) >= 2)
  filename = ARGS[1]
else
  filename ="cube.1.mesh" 
end
dim, eptr, eind, elmdist  = read_par_mesh(filename, Int32(0)) 
#map([(eptr, "eptr.$rank"), (eind, "eind.$rank"), (elmdist, "elmdist.$rank")])  do x
#  io = open(x[2], "w")
#  writedlm(io,x[1])
#end
io = open("eptr.$rank", "w")
writedlm(io, eptr)
@info "eptr_size=", size(eptr,1)
io = open("eind.$rank", "w")
writedlm(io, eind)
io = open("elmdist.$rank", "w")
writedlm(io, elmdist)
println("I'm $rank, nodes ∈ [$(minimum(eind)), $(maximum(eind))]")
println("I'm $rank, elements ∈ [$(elmdist[rank+1]), $(elmdist[rank+2])]")
# call parmetis
if (rank == 0)
  println("computing dual mesh by parmetis")
end
t0 = time()
xadj, adjcy = parmetis_mesh_to_dual(;elmdist=elmdist, eptr=eptr, eind=eind, baseval=Int32(0), ncommon=Int32(2), comm = comm)
@printf "tps parmetis %.3e\n" time() - t0
adj_check = metis_fmt_to_vector(xadj, adjcy, Int32(0))
# call ptscotchparmetis
# our algorithm in pure julia
if (rank == 0)
  println("computing dual mesh by scotchparmetis")
end
t0 = time()
xadj2, adjcy2 = ptscotchparmetis_mesh_to_dual(;elmdist=elmdist, eptr=eptr, eind=eind, baseval=Int32(0), ncommon=Int32(2), comm = comm)
@printf "tps ptscotchparmetis %.3e\n" time() - t0
adj2_check = metis_fmt_to_vector(xadj2, adjcy2, Int32(0))

if (rank == 0)
  println("checking")
end
@test length(adj_check) == length(adj_check)
for i=1:length(adj_check)
  @test sort(adj_chek[i]) == sort(adj_check)
end
MPI.Finalize()
@test MPI.Finalized()
