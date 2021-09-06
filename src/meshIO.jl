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
  if (r == 0)
    @info "nb_elem=", nb_elem
  end
  elmdist  = zeros(IndexType, p+1)

  elements = zeros(IndexType, nb_elem, dim + 1)
  for i=1:nb_elem
    elements[i, : ] = map(x->parse(IndexType,x), split(lines[of_elem + 1 + i])[1:dim+1])
  end
  nMin = convert(typeof(length(IndexType[0])), minimum(elements)) 
  @assert nMin == 1
  # build elmdist
  startIndex, endIndex, _ = tile(nb_elem, 0)
  elmdist[1] = startIndex
  elmdist[2] = endIndex
  for i=2:p
    startIndex, endIndex, _ = tile(nb_elem, i-1)
    elmdist[i+1] = elmdist[i] + endIndex-startIndex + 1
  end 
  startIndex, endIndex, sizeChunk = tile(nb_elem, r)
  eind = zeros(IndexType, endIndex - startIndex + 1 , dim + 1)
  @info "I'm" , r,  startIndex, endIndex
  for i=startIndex:endIndex 
    eind[i-startIndex+1,:] =  map(x->parse(IndexType,x), split(lines[of_elem + i + 1])[1:dim+1])
  end
  eptr = [IndexType(0):IndexType(dim+1):IndexType(sizeChunk * (dim + 1));]
  eind = eind'[:]
  return (eptr, eind, elmdist)
end
