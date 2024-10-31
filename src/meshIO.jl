# Reader functions
"""
function read_par_mesh(filename::String, baseval::T) where {T<:Integer}
reads a medit mesh file in parallel
"""
function read_par_mesh(filename::String, baseval::T) where {T<:Integer}
  IndexT = typeof(length(T[0]))
  comm, p, r = get_com_size_rank()
  lines = readlines(filename)
  lines = lines[map(x-> length(x) > 0 && (x[1] != '#'), lines)]
  f_s(x) = findall(map(y->occursin(x,y),lines))[1]
  dim = parse(T,lines[f_s("Dimension")+1])
  if (r == 0)
    @info "dim =", dim
  end
  of_nodes = f_s("Vertices")
  nb_nodes = parse(T, lines[of_nodes + 1])
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
  nb_elem = parse(T, lines[of_elem + 1])
  print_master("nb_elem=" * nb_elem)

  elements = zeros(T, nb_elem, dim + 1)
  for i=1:nb_elem
    elements[i, : ] = map(x->parse(T,x), split(lines[of_elem + 1 + i])[1:dim+1])
  end
  nMin = convert(IndexT, minimum(elements)) 
  elements .-= nMin;
  nMin = convert(IndexT, minimum(elements)) # FIXME : why ????
  nMax = convert(IndexT, maximum(elements)) 

  # build elmdist
  elmdist  = zeros(T, p+1)
  startIndex, endIndex, _ = tile(nMax, 0, nMin, IndexT(p))
  elmdist[1] = startIndex
  elmdist[2] = endIndex
  for i=2:p
      startIndex, endIndex, _ = tile(nMax, i-1, nMin, IndexT(p))
    elmdist[i+1] = elmdist[i] + endIndex-startIndex + 1
  end 
  startIndex, endIndex, sizeChunk = tile(nMax, r, nMin, IndexT(p))
  eind = zeros(T, endIndex - startIndex + 1 , dim + 1)
  for i=startIndex:endIndex 
      eind[i-startIndex+1,:] = map(x->parse(T,x)-1, split(lines[of_elem + i + 2])[1:dim+1])
  end
  eptr = [T(0):T(dim+1):T(sizeChunk * (dim + 1));]
  eind = eind'[:]
  return dim, eptr, eind, elmdist
end

function read_par_dmh(filename::String, T::DataType = Int32)
  files, basename = parse_file_name(filename)
  com, size, rank = get_com_size_rank()
  if (length(files) != size)
    @error "bad number of MPI processes, must be equal to", length(files)
  end
  nb_files = length(files)
  elmdist = Vector{T}(undef, size+1)
  baseval = 0
  fIO = open(files[rank+1])
  version = extract!(fIO, T)[1]
  procglbnum, proclocnum = extract_tuple!(fIO, 2)
  elmglbnum = extract!(fIO, T)[1]
  elmlocnum, vertlocnum = extract_tuple!(fIO, 2)
  baseval, chaco = extract_tuple!(fIO, 2)   
  elmdist = extract!(fIO, T)[2:end]
  j = 1 
  eptr = T[baseval]
  eind = T[] 
  while !eof(fIO)
    tab = extract!(fIO, T)
    @assert tab[1] == length(tab[2:end])
    e = tab[2:end]
    append!(eptr, eptr[j]+length(e))
    append!(eind, e)
    j += 1
  end
  close(fIO)
  return eptr, eind, elmdist, baseval, basename
end

"""
function read_dmh(filename, T)
reads sequentially a distributed mesh scotch file
"""
function read_dmh(filename::String, T::DataType = Int32)
  files, basename = parse_file_name(filename)
  nb_files = length(files)
  elmlocnum = zeros(T, nb_files)
  vertlocnum = zeros(T, nb_files)
  eptr = Vector{Vector{T}}(undef, nb_files)
  eind = Vector{Vector{T}}(undef, nb_files)
  elmdist = T[]
  baseval = 0
  for (i,f) in enumerate(files)
    fIO = open(f)
    version = extract!(fIO, T)[1]
    procglbnum, proclocnum = extract_tuple!(fIO, 2)
    elmglbnum = extract!(fIO, T)[1]
    elmlocnum[i], vertlocnum[i] = extract_tuple!(fIO, 2)
    baseval, chaco = extract_tuple!(fIO, 2)   
    elmdist = extract!(fIO, T)[2:end]
    j = 1 
    eptr[i] = T[T(baseval)]
    eind[i] = T[] 
    while !eof(fIO)
      tab = extract!(fIO, T)
      @assert tab[1] == length(tab[2:end])
      e = tab[2:end]
      append!(eptr[i], eptr[i][j]+length(e))
      append!(eind[i], e)
      j += 1
    end
    close(fIO)
  end
  return eptr, eind, elmdist, baseval, basename
end

"""
read\\_msh(filename, T) 
read sequentially a scotch mesh
"""
function read_msh(filename, T::DataType = Int32)
  io = open(filename, "r")
  format_version = extract!(io, T)[1]
  @assert format_version == 1
  velmnbr, vnodnbr, edgenbr = extract_tuple!(io, 3, T)
  velmbas, vnodbas, chaco_flags = extract_tuple!(io, 3, T)
  if (velmbas > vnodbas)
    for i=1:vnodnbr
      readline(io)
    end 
  end
  eptr = T[vnodbas]
  eind = T[]
  j = 1
  for i = velmbas : velmbas + velmnbr -1 
    tab = extract!(io, T)
    e = tab[2:end]
    @assert tab[1] == length(e)
    append!(eptr, eptr[j]+length(e))
    append!(eind, e)
    j+=1
  end
  close(io)
  eptr .-= vnodbas
  eind .-= vnodbas
  return eptr, eind
end

#Writers

function write_par_dmh(;eptr::Vector{T}, eind::Vector{T}, elmdist::Vector{T}, baseval::T, filename::String) where {T}
  comm, size, rank = get_com_size_rank()
  io = open(filename * "-$rank" * ".dmh", "w")
  format_version = 0
  elmlocnbr = elmdist[rank+2] - elmdist[rank+1]
  nMax = Ref{T}(maximum(eind))
  MPI.Allreduce!(nMax, MPI.MAX, comm)
  nMin = Ref{T}(minimum(eind))
  MPI.Allreduce!(nMin, MPI.MIN, comm)
  @info baseval, nMin
  @assert baseval == nMin[]
  
  vertglbnbr = nMax[] - baseval + 1
  println(io, format_version) # write version
  println(io, size, " ", rank) # write procglbnum and procglbnim
  println(io, elmdist[end], " ", vertglbnbr) # write total number of elements 
  println(io, elmlocnbr, " ", eptr[end]) # write local number of elements and size of vertlocnbr
  println(io, baseval, " 000") # write baseval and chaco code
  print(io, length(elmdist), " ")
  for elm=elmdist[1:end-1]
    print(io, elm, " ")
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


"""
write_dmh(;eptr, eind, elmdist, baseval, filename) 
write sequentially a distributed mesh in a file sequence
"""
function write_dmh(;eptr::Vector{Vector{T}}, eind::Vector{Vector{T}}, elmdist::Vector{T}, baseval::T, filename::String) where {T}
  size = length(elmdist) - 1
  nMax = eind[1][1]
  nMin = nMax
  for rank=0:size-1
      nMax = max(nMax, maximum(eind[rank+1])) 
      nMin = min(nMin, minimum(eind[rank+1])) 
  end
  @assert nMin == baseval
  vertglbnbr = nMax - baseval + 1
  for rank=0:size-1
    @assert baseval == eptr[rank+1][1]
    io = open(filename * "-$rank" * ".dmh", "w")
    format_version = 0
    elmlocnbr = elmdist[rank+2] - elmdist[rank+1]
    println(io, format_version) # write version
    println(io, size, " ", rank) # write procglbnum and procglbnim
    println(io, elmdist[end] - baseval, " ", vertglbnbr) # write total number of elements 
    println(io, elmlocnbr, " ", eptr[rank+1][end] - eptr[rank+1][1]) # write local number of elements and size of vertlocnbr
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
