using Mesh2Dual
using MPI
using Test
using Printf
using DelimitedFiles
using Comonicon

Comonicon.@main function shift_dmh(;filename::String, shift::Int32=Int32(0))

  eptr, eind, elmdist, baseval, basename  = read_dmh(filename, Int32) 
  if shift != Int32(0)
    shift_msh!(;elmdist=elmdist, eind=eind, eptr=eptr, shift)
  end
  write_dmh(;elmdist=elmdist, eptr=eptr, eind=eind, baseval=baseval+shift, filename=basename*"_"*string(shift))
  return 0
end
