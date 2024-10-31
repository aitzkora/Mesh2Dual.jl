using Mesh2Dual
using MPI
using Test
using Printf
using DelimitedFiles
using Comonicon

Comonicon.@main function msh2dmh(;filename::String, parts::Int64 = 3)
  name_we, ext = splitext(filename)
  @assert ext == ".msh"
  eptr, eind = read_msh(filename, Int32) 
  eptr_p, eind_p, elmdist = split_msh2(eptr, eind, parts)
  write_dmh(;elmdist=elmdist, eptr=eptr_p, eind=eind_p, baseval=Int32(0), filename=name_we)
  return 0
end
