function read_par_mesh(filename, ::Val{D}, baseval) where {D}
  IndexType = typeof(baseval)
  @assert Val(D) isa Union{map(x->Val{x},1:3)...}
  comm = MPI.COMM_WORLD
  p = MPI.Comm_size(comm)
  r = MPI.Comm_rank(comm)
  
  lines = readlines(filename)
  lines = lines[map(x-> length(x) > 0 && (x[1] != '#'), lines)] # beware to the '
  f_s(x) = findall(map(y->occursin(x,y),lines))[1]
  dim = parse(IndexType,lines[f_s("Dimension")+1])
  if (r == 0)
    @info "dim =", dim
  end
  @assert  D == dim
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
  return (eptr, eind, elmdist)
end


"""
shift a mesh by a constant 
"""
function shift_msh(;elmdist::Vector{T}, eptr::Vector{T}, eind::Vector{T}, shift::T) where {T}
  comm = MPI.COMM_WORLD
  size = MPI.Comm_size(comm)
  rank = MPI.Comm_rank(comm)
  elmlocnbr = elmdist[rank+2] - elmdist[rank+1]
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
    for idx=eptr[i]+1:eptr[i+1]-1
      print(io, eind[idx], " ")
    end
    println(io,eind[eptr[i+1]])
  end
end

function write_dgraph(;xadj::Vector{T}, adjncy::Vector{T}, elmdist::Vector{T}, baseval::T, comm::MPI.Comm, filename::String) where {T}
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
end     
        
