using Mesh2Dual
using MPI
using Test
using Printf
using DelimitedFiles
using Comonicon

Comonicon.@main function compute_dual(;filename::String="cube.1.mesh", shift::Int32=Int32(0))
  MPI.Init()
  comm, size, rank = get_com_size_rank()

  #print_master("loading mesh")

  dim, eptr, eind, elmdist  = read_par_mesh(filename, Int32(0)) 
  println("$rank -> elements ∈ [$(elmdist[rank+1]), $(elmdist[rank+2])]")
  t0 = time()
  #print_master("computing dual mesh...")
  adj = dgraph_dual(;elmdist=elmdist, eptr=eptr, eind=eind, baseval=Int32(0), ncommon=Int32(3), comm = comm) 
  tps_dual = time() - t0
  #print_master("tps pure julia $(@sprintf("%.2fs", tps_dual))")
  xadj, adjncy = list_to_csr(adj)
  basename = split(filename,".")[1]*""
  write_par_dgraph(;xadj=xadj, adjncy=adjncy, elmdist=elmdist, baseval=Int32(0), filename=basename)
  # write mesh
  if shift != Int32(0)
    shift_par_msh!(;elmdist=elmdist, eind=eind, eptr=eptr, shift)
  end
  write_par_dmh(;elmdist=elmdist, eptr=eptr, eind=eind, baseval=shift, filename=basename*"_"*string(shift))
  MPI.Finalize()
  return
end
