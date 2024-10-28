# Mesh2Dual.jl
This julia package aims to make experiments with some algorithms on meshes to compute dual graphs, i.e. also called graphs of elements, to implement after the corresponding algorithm in the Scotch software. Some functions manipulating meshes and graphs could be used for others purposes

## Compile the package

```julia
using Pkg
Pkg.activate("path/to/Mesh2Dual.jl")
using Mesh2Dual
```

## Testing
```julia
Pkg.test("Mesh2Dual")
```

## compute a dual graph of a mesh (medit file) in parallel

```bash
mpirun -np 3 julia mesh2dgr.jl --filename mesh_oliv.mesh --shift 1
```

## compute a dual graph of a dmh (distributed mesh scotch file) to a scotch distributed graph

```bash
mpirun -np 3 julia dmh2dgr.jl --filename data/ship001-%r.dmh
```

## How to check if two distributed graphs are equal

```bash
check.jl mesh_oliv-%r.dgr mesh_0-%r.dgr
```
Beware, the two dgr files must be in the working directory (or calling)

# Dependencies
- to run some tests you need of shared versions of `libmetis` and `libparmetis` exported with their paths into  `METIS_LIB` and `PARMETIS_LIB` environments variables
