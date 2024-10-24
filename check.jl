#!/usr/bin/env julia
using Mesh2Dual

function check(::Type{T}, f1::String, f2::String, baseval::T = zero(T)) where {T<:Integer}
  println("comparing $f1 and $f2")
  io1 = open(f1)
  io2 = open(f2)
  h1 = read_dgraph_header(io1)
  h2 = read_dgraph_header(io2)
  adj1 = read_adj(io1)
  adj2 = read_adj(io2)
  close(io1)
  close(io2)
  for i=1:length(adj1)
    if (sort(adj1[i]) != sort(adj2[i]) .+ baseval)
      println("$f1 and $f2 are differents")
      exit(-1)
    end
  end
end

if (length(ARGS) < 2)
  println("usage : julia check.jl fileIndexedOne fileIndexedTwo [baseval]")
  exit(-1) 
end

f1 = parse_file_name(ARGS[1])
f2 = parse_file_name(ARGS[2])
if (length(ARGS) >= 3)
  baseval = parse(Int64, ARGS[3])
end
if (isempty(f1))
  println(ARGS[1] * " is not in the current dir !"); exit(-1)
end 
if (isempty(f2))
  println(ARGS[2] * " is not in the current dir !"); exit(-1)
end 

if (length(f1) != length(f2))
  println("two generic filenames do not have the same number of processes"); exit(-1)
else
  for i=1:length(f1)
    check(Int32, f1[i], f2[i], Int32(baseval))
  end
  println("****** SUCCESS *****");
end
