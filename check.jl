#!/usr/bin/env julia
struct header
  version :: Int32
  procglbnum :: Int32
  proclocnum :: Int32
  vertglbnum :: Int32
  edgeglbnum :: Int32
  vertlocnum :: Int32
  edgelocnum :: Int32
  baseval    :: Int32
  chaco      :: Int32
end

function extract!(s::IO, n::Int64)
  tab = split(readline(s))
  return tuple((x->parse(Int32,x)).(tab[1:n])...)
end

function read_header(s::IO)
  version = extract!(s, 1)[1]
  procglbnum, proclocnum = extract!(s, 2)
  vertglbnum, edgeglbnum = extract!(s, 2)
  vertlocnum, edgelocnum = extract!(s, 2)
  baseval, chaco = extract!(s, 2)   
  return header(version, procglbnum, proclocnum, vertglbnum, edgeglbnum, vertlocnum, edgelocnum, baseval, chaco)
end

function read_adj(s::IO)
  adj = Vector{Int32}[]
  i = 5
  while !eof(s)
    linea = split(readline(s), "\t")
    tab = (x->parse(Int32,x)).(linea)
    if tab[1] != length(tab[2:end])
        println("mauvaise taille", " ", i, " ", linea)
    else
      push!(adj, tab[2:end])
    end
    i+=1
  end
  return adj
end

function check(f1, f2)
  println("comparing $f1 and $f2")
  io1 = open(f1)
  io2 = open(f2)
  h1 = read_header(io1)
  h2 = read_header(io2)
  adj1 = read_adj(io1)
  adj2 = read_adj(io2)
  close(io1)
  close(io2)
  for i=1:length(adj1)
    if (sort(adj1[i]) != sort(adj2[i]))
      println("$f1 and $f2 are differents")
      exit(-1)
    end
  end
end

function parse_file_name(str)
  pattern_rng = findfirst("%r", str)
  if isnothing(pattern_rng)
    @warn "file $str does not contains %r"
    return
  end
  fst_part = str[1:pattern_rng[1]-1]
  snd_part = str[pattern_rng[2]+1:end]
  p = 0
  files = String[] 
  for s in readdir()
      rx = Regex(fst_part*"(\\d)+"*snd_part)
      mx = match(rx,s)
      if !isnothing(mx)
          p+=1
          push!(files, fst_part*mx.captures[1]*snd_part)
      end
  end   
  return files
end

if (length(ARGS) < 2)
    println("usage : julia check.jl fileIndexedOne fileIndexedTwo")
    exit(-1) 
end

f1 = parse_file_name(ARGS[1])
f2 = parse_file_name(ARGS[2])
if (isempty(f1))
    println(ARGS[1]* " is not in the current dir !"); exit(-1)
end 
if (isempty(f2))
    println(ARGS[2]* " is not in the current dir !"); exit(-1)
end 

if (length(f1) != length(f2))
    println("two generic filenames do not have the same number of processes"); exit(-1)
else
    for i=1:length(f1)
        check(f1[i], f2[i])
    end
    println("****** SUCCESS *****");
end
