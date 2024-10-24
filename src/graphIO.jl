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
  version = extract!(s)[1]
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
