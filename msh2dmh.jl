using Mesh2Dual
using MPI
using Test
using Printf
using DelimitedFiles
using Comonicon

Comonicon.@main function msh2dmh(;filename::String="data/ship001", 
                                  parts::Int64 = 3)
  eptr, eind = read_msh(filename, Int32) 
  eptr_p, eind_p, elmdist = split_msh2(eptr, eind, parts)
  write_dmh(;elmdist=elmdist, eptr=eptr_p, eind=eind_p, baseval=Int32(0), filename=filename)
  return 0
end
