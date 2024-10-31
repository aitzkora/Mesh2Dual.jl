using Mesh2Dual
using MPI
using Test
using Printf
using DelimitedFiles
using Comonicon

Comonicon.@main function compute_dual(;filename::String="cube.1.mesh")
  MPI.Init()
  comm, size, rank = get_com_size_rank()
  print_master("loading mesh")
  dim, eptr, eind, elmdist  = read_par_mesh(filename, Int32(0)) 
  t0 = time()
  print_master("computing dual mesh...")
  adj = dgraph_dual(;elmdist=elmdist, eptr=eptr, eind=eind, baseval=Int32(0), ncommon=Int32(3), comm = comm) 
  tps_dual = time() - t0
  print_master("tps pure julia $(@sprintf("%.2fs", tps_dual))")
  xadj, adjncy = list_to_csr(adj)
  basename = split(filename,".")[1]*""
  print_master("writing mesh...")
  write_par_dgraph(;xadj=xadj, adjncy=adjncy, elmdist=elmdist, baseval=Int32(0), filename=basename)
  MPI.Finalize()
  return
end
