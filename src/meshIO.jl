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
  for s in readdir() # FIXME : works only in current dir
      rx = Regex(fst_part*"(\\d)+"*snd_part)
      mx = match(rx,s)
      if !isnothing(mx)
          p+=1
          push!(files, fst_part*mx.captures[1]*snd_part)
      end
  end   
  return files
end

function read_par_mesh(filename, baseval)
  IndexType = typeof(baseval)
  comm = MPI.COMM_WORLD
  p = MPI.Comm_size(comm)
  r = MPI.Comm_rank(comm)
  lines = readlines(filename)
  lines = lines[map(x-> length(x) > 0 && (x[1] != '#'), lines)]
  f_s(x) = findall(map(y->occursin(x,y),lines))[1]
  dim = parse(IndexType,lines[f_s("Dimension")+1])
  if (r == 0)
    @info "dim =", dim
  end
  of_nodes = f_s("Vertices")
  nb_nodes = parse(IndexType, lines[of_nodes + 1])
  nodes = zeros(Float64, nb_nodes, dim)
  for i=1:nb_nodes
    nodes[i,:] = map(x->parse(Float64,x), split(lines[of_nodes + 1 + i])[1:end-1])
  end
  if (r == 0)
    @info "nb_nodes=", nb_nodes
  end
  if (dim == 2)
    of_elem = f_s("Triangles")
  else
    of_elem = f_s("Tetrahedra")
  end
  nb_elem = parse(IndexType, lines[of_elem + 1])
  nb_elem = convert(typeof(length(IndexType[0])), nb_elem)
  if r == 0
    @info "nb_elem=", nb_elem
  end
  elmdist  = zeros(IndexType, p+1)
  elements = zeros(IndexType, nb_elem, dim + 1)
  for i=1:nb_elem
    elements[i, : ] = map(x->parse(IndexType,x), split(lines[of_elem + 1 + i])[1:dim+1])
  end
  nMin = convert(typeof(length(IndexType[0])), minimum(elements)) 
  elements .-= nMin;
  nMin = convert(typeof(length(IndexType[0])), minimum(elements)) 
  nMax = convert(typeof(length(IndexType[0])), maximum(elements)) 
  #@assert nMin == 1
  # build elmdist
  startIndex, endIndex, _ = tile(nMax, 0, nMin)
  elmdist[1] = startIndex
  elmdist[2] = endIndex
  for i=2:p
    startIndex, endIndex, _ = tile(nMax, i-1, nMin)
    elmdist[i+1] = elmdist[i] + endIndex-startIndex + 1
  end 
  startIndex, endIndex, sizeChunk = tile(nMax, r, nMin)
  eind = zeros(IndexType, endIndex - startIndex + 1 , dim + 1)
  for i=startIndex:endIndex 
      eind[i-startIndex+1,:] = map(x->parse(IndexType,x)-1, split(lines[of_elem + i + 2])[1:dim+1])
  end
  eptr = [IndexType(0):IndexType(dim+1):IndexType(sizeChunk * (dim + 1));]
  eind = eind'[:]
  return dim, eptr, eind, elmdist
end


"""
shift_par_msh!(;elmdist::Vector{T}, eptr::Vector{T}, eind::Vector{T}, shift::T) where {T}
shifts a mesh in place by a constant 
"""
function shift_par_msh!(;elmdist::Vector{T}, eptr::Vector{T}, eind::Vector{T}, shift::T) where {T}
  comm = MPI.COMM_WORLD
  size = MPI.Comm_size(comm)
  rank = MPI.Comm_rank(comm)
  for i in eachindex(elmdist)
    elmdist[i] += shift
  end
  for i in eachindex(eptr)
    eptr[i] += shift
  end
  for i in eachindex(eind)
    eind[i] += shift
  end
end

"""
shift_msh!(;elmdist::Vector{T}, eptr::Vector{T}, eind::Vector{T}, shift::T) where {T}
shifts a mesh in place by a constant 
"""
function shift_msh!(; eptr::Vector{Vector{T}}, eind::Vector{Vector{T}}, elmdist::Vector{T}, shift::T) where {T}
  size = length(elmdist) - 1
  for i in eachindex(elmdist)
    elmdist[i] += shift
  end
  for rank=0:size-1
    for i in eachindex(eptr[rank+1])
        eptr[rank+1][i] += shift
    end
    for i in eachindex(eind[rank+1])
        eind[rank+1][i] += shift
    end
  end
end

function write_par_dmesh(;eptr::Vector{T}, eind::Vector{T}, elmdist::Vector{T}, baseval::T, filename::String) where {T}
  comm = MPI.COMM_WORLD
  size = MPI.Comm_size(comm)
  rank = MPI.Comm_rank(comm)
  io = open(filename * "-$rank" * ".dmh", "w")
  format_version = 0
  elmlocnbr = elmdist[rank+2] - elmdist[rank+1]
  println(io, format_version) # write version
  println(io, size," ", rank) # write procglbnum and procglbnim
  println(io, elmdist[end]) # write total number of elements 
  println(io, elmlocnbr, " ", eptr[end]) # write local number of elements and size of vertlocnbr
  println(io, baseval, " 000") # write baseval and chaco code
  print(io, length(elmdist), " ")
  for elm=elmdist[1:end-1]
      print(io, elm," ")
  end
  println(io, elmdist[end])
  for i=1:length(eptr)-1
      print(io, eptr[i+1] - eptr[i], " ")
      for idx=eptr[i]-baseval+1:eptr[i+1]-baseval-1
          print(io, eind[idx], " ")
      end
      println(io,eind[eptr[i+1]-baseval])
  end
end

function read_dmesh(x::Val{T}, filename::String) where {T}
  files = parse_file_name(filename)
  nb_files = length(files)
  println("nbprocs = $nb_files")
  elmlocnum = zeros(T, nb_files)
  vertlocnum = zeros(T, nb_files)
  eptr = [T[] for _ in 1:nb_files]
  eind = [T[] for _ in 1:nb_files]
  elmdist = []
  baseval = 0
  for (i,f) in enumerate(files)
    fIO = open(f)
    version = extract!(fIO, 1)[1]
    procglbnum, proclocnum = extract!(fIO, 2)
    elmglbnum = extract!(fIO, 1)[1]
    elmlocnum[i], vertlocnum[i] = extract!(fIO, 2)
    baseval, chaco = extract!(fIO, 2)   
    linea = split(readline(fIO))
    tab= (x->parse(Int32,x)).(linea)
    elmdist = tab[2:end]
    elmlocnbr = elmdist[i] - elmdist[i]
    j = 1 
    eptr[i] = T[T(baseval)]
    eind[i] = T[] 
    while !eof(fIO)
      linea = split(readline(fIO))
      tab=(x->parse(Int32,x)).(linea)
      @assert tab[1] == length(tab[2:end])
      e = tab[2:end]
      append!(eptr[i], eptr[i][j]+length(e))
      append!(eind[i], e)
      j += 1
    end
    close(fIO)
  end
  return eptr, eind, elmdist, baseval
end

"""
write_dmesh(;eptr::Vector{T}, eind::Vector{T}, elmdist::Vector{T}, baseval::T, filename::String) where {T}

write sequentially a distributed mesh in a file sequence
"""
function write_dmesh(;eptr::Vector{Vector{T}}, eind::Vector{Vector{T}}, elmdist::Vector{T}, baseval::T, filename::String) where {T}
  size = length(elmdist) - 1
  for rank=0:size-1
    io = open(filename * "-$rank" * ".dmh", "w")
    format_version = 0
    elmlocnbr = elmdist[rank+2] - elmdist[rank+1]
    println(io, format_version) # write version
    println(io, size," ", rank) # write procglbnum and procglbnim
    println(io, elmdist[end]) # write total number of elements 
    println(io, elmlocnbr, " ", eptr[rank+1][end]) # write local number of elements and size of vertlocnbr
    println(io, baseval, " 000") # write baseval and chaco code
    print(io, length(elmdist), " ")
    for elm=elmdist[1:end-1]
        print(io, elm," ")
    end
    println(io, elmdist[end])
    for i=1:length(eptr[rank+1])-1
        print(io, eptr[rank+1][i+1] - eptr[rank+1][i], " ")
        for idx=eptr[rank+1][i]-baseval+1:eptr[rank+1][i+1]-baseval-1
            print(io, eind[rank+1][idx], " ")
        end
        println(io,eind[rank+1][eptr[rank+1][i+1]-baseval])
    end
    close(io)
  end
end

function write_par_dgraph(;xadj::Vector{T}, adjncy::Vector{T}, elmdist::Vector{T}, baseval::T, filename::String) where {T}
  comm = MPI.COMM_WORLD
  size = MPI.Comm_size(comm)
  rank = MPI.Comm_rank(comm)
  io = open(filename * "-$rank" * ".dgr", "w")
  format_version = 2
  elmlocnbr = elmdist[rank+2] - elmdist[rank+1]
  vertglbnbr = Ref{T}(xadj[elmlocnbr])
  edgelocnbr = 0
  for i=1:(length(xadj)-1)
    for n ∈ adjncy[1+xadj[i]:xadj[i+1]]
        edgelocnbr = edgelocnbr + 1
    end 
  end 
  edgeglbnbr = Ref{T}(edgelocnbr)
  MPI.Allreduce!(edgeglbnbr, MPI.SUM, comm)
  println(io, format_version) # write version
  println(io, size,"\t", rank) # write procglbnum and procglbnim
  println(io, elmdist[end], "\t",  edgeglbnbr[]) # write global vertices and edges number
  println(io, elmlocnbr, "\t", edgelocnbr) # write local number of elements and size of vertlocnbr
  println(io, baseval, "\t000") # write baseval and chaco code
  for i=1:length(xadj)-1
    sizeLoc=xadj[i+1]-xadj[i]
    if (sizeLoc > 0)
      print(io, xadj[i+1] - xadj[i], "\t")
      for idx=xadj[i]+1:xadj[i+1]-1
        print(io, adjncy[idx], "\t")
      end
      println(io, adjncy[xadj[i+1]])
    else
      println(io,0)
    end
  end   
  close(io) 
end     
        
