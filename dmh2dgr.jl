using Mesh2Dual
using MPI
using Test
using Printf
using DelimitedFiles
using Comonicon

Comonicon.@main function compute_dual(;filename::String, ncommon::Int64)
  MPI.Init()
  comm, size, rank = get_com_size_rank()
  eptr, eind, elmdist, baseval, basename  = read_par_dmh(filename, Int32) 
  print_master("baseval = $baseval\n")
  print_master("computing dual mesh...")
  t0 = time()
  adj = dgraph_dual(;elmdist=elmdist, eptr=eptr, eind=eind, baseval=Int32(baseval), ncommon=Int32(ncommon), comm = comm) 
  tps_dual = time() - t0
  print_master("tps pure julia $(@sprintf("%.2fs", tps_dual))")
  xadj, adjncy = list_to_csr(adj, baseval)
  print_master("writing file...")
  write_par_dgraph(;xadj=xadj, adjncy=adjncy, elmdist=elmdist, baseval=Int32(baseval), filename=basename)
  MPI.Finalize()
  return
end
