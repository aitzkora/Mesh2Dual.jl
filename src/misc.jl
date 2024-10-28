"""
parse_file_name(str::String)
find all files matching with the regular expression where %r stands for rank number
"""
function parse_file_name(str::String)
  pattern_rng = findfirst("%r", str)
  if isnothing(pattern_rng)
    @warn "file $str does not contains %r"
    return
  end
  fst_part = str[1:pattern_rng[1]-1]
  snd_part = str[pattern_rng[2]+1:end]
  p = 0
  files = String[] 
  dir = dirname(str)
  for s in readdir(dir, join=true)
      rx = Regex(fst_part*"(\\d)+"*snd_part)
      mx = match(rx,s)
      if !isnothing(mx)
          p+=1
          push!(files, fst_part*mx.captures[1]*snd_part)
      end
  end   
  return files, fst_part
end

"""
list_to_csr(adj::Vector{Vector{T}}) 
converts a adjacency list to a csr pair based at 0
"""
function list_to_csr(adj::Vector{Vector{T}}) where {T}
    xadj = T[0]
    adjncy = T[] 
    for (i, e) in enumerate(adj)
        append!(xadj, xadj[i]+length(e))
        append!(adjncy, e)
    end
    return (xadj, adjncy)
end 


"""
csr_to_list(eptr, eind)
converts csr pair based at 0 (eptr, eind) to an adjacency list
"""
function csr_to_list(eptr::Vector{T}, eind::Vector{T}) where {T}
  n = length(eptr)
  adj = Vector{Vector{T}}(undef, n  - 1)
  for i=1:n-1
    adj[i] = eind[eptr[i]+1:eptr[i+1]]
  end
  return adj
end 

"""
extract!(s::IO, T::DataType)
extract a Vector{T} in a text line in s 
"""

function extract!(s::IO, T::DataType = Int32)
  tab = split(readline(s))
  return (x->parse(T,x)).(tab[1:end])
end

"""
extract_tuple!(s::IO, n::Int64, T::DataType)
extract a text line in s and convert it to a n-tuple of base type T

example
========

we want to read "100 12 2" in a line of a file `s`
z = extract_tuple!(s, 3, Int32) -> z = (100, 12, 2)
"""
extract_tuple!(s::IO, n::Int64, T::DataType = Int32) = tuple(extract!(s, T)[1:n]...)

