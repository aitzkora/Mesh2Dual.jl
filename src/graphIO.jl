struct dgraph_header
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

function read_dgraph_header(s::IO)
  version = extract!(s, Int32)[1]
  procglbnum, proclocnum = extract_tuple!(s, 2)
  vertglbnum, edgeglbnum = extract_tuple!(s, 2)
  vertlocnum, edgelocnum = extract_tuple!(s, 2)
  baseval, chaco = extract_tuple!(s, 2)   
  return dgraph_header(version, procglbnum, proclocnum, vertglbnum, edgeglbnum, vertlocnum, edgelocnum, baseval, chaco)
end

function read_adj(s::IO)
  adj = Vector{Int32}[]
  i = 5
  while !eof(s)
    tab = extract!(s, Int32, "\t")
    @assert tab[1] == length(tab[2:end])
    push!(adj, tab[2:end])
  end
  return adj
end

function write_par_dgraph(;xadj::Vector{T}, adjncy::Vector{T}, elmdist::Vector{T}, baseval::T, filename::String) where {T}
  comm, size, rank = get_com_size_rank()
  io = open(filename * "$rank" * ".dgr", "w")
  format_version = 2
  elmlocnbr = elmdist[rank+2] - elmdist[rank+1]
  vertglbnbr = Ref{T}(xadj[elmlocnbr])
  edgelocnbr = 0
  for i=1:(length(xadj)-1)
    for n ∈ adjncy[1+xadj[i]-baseval:xadj[i+1]-baseval]
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
        print(io, adjncy[idx-baseval], "\t")
      end
      println(io, adjncy[xadj[i+1]-baseval])
    else
      println(io,0)
    end
  end   
  close(io) 
end
