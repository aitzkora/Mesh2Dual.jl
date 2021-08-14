using Test
using Mesh2Dual

@testset "gen_parts_test" begin
   l = [1, 2, 4]
   @test Set(gen_parts(l, 2)) == Set(Set.([(1,2),(1,4),(2,4)]))
   @test Set(gen_parts(l, 1)) == Set(Set.([(1),(2),(4)])) 
   @test Set(gen_parts(l, 3)) == Set(Set.([(1,2,4)])) 
end

@testset "mesh_to_metis" begin

  m_dat = reshape(Int32[1,2,3, 1,2,6, 2,6,5, 2,5,7, 2,7,4, 2,4,3], 3, 6)
  m = Mesh{Int32}([m_dat[:, j] for j=1:size(m_dat,2)])
  @test mesh_to_metis_fmt(m) == (Int32[0,3,6,9,12,15,18],Int32[0,1,2,0,1,5,1,5,4,1,4,6,1,6,3,1,3,2], Int32(1))
  @test metis_fmt_to_vector(mesh_to_metis_fmt(m)...) == m.elements

end

@testset "graph_dual" begin

    m_dat = reshape([1,2,3, 1,2,6, 2,6,5, 2,5,7, 2,7,4, 2,4,3], 3, 6)
    m = Mesh{Int64}([m_dat[:, j] for j=1:size(m_dat,2)])

    @test graph_dual(m, 1) == Graph([setdiff([1:6;],[i]) for i=1:6]) 
    @test graph_dual(m, 3) == Graph([ Int64[] for _ in 1:6])
    @test graph_dual(m, 2) == Graph([[2,6], [1,3], [2,4], [3,5], [4,6], [1,5]])
end

@testset "graph_dual_new" begin
  if !haskey(ENV, "METIS_LIB") 
    @error "you must define METIS_LIB in your variable environment"
  end
 
  m_dat = reshape(Int32[1,2,3, 1,2,6, 2,6,5, 2,5,7, 2,7,4, 2,4,3], 3, 6)
  m = Mesh{Int32}([m_dat[:, j] for j=1:size(m_dat,2)])
 
  ptr, ind, ne  = mesh_to_scotch_fmt(m)
  
  g1 = graph_dual_new(m, 1)
  g2 = graph_dual_new(m, 2)
  g3 = graph_dual_new(m, 3)
   
  m1 = metis_graph_dual(m, 1)
  m2 = metis_graph_dual(m, 2)
  m3 = metis_graph_dual(m, 3)
  
  m1p = graph_dual(m, 1);
  m2p = graph_dual(m, 2);
  m3p = graph_dual(m, 3);
  
  @assert map(x->x .- 1, m1p.adj) == g1
  @assert map(x->x .- 1, m2p.adj) == g2
  @assert map(x->x .- 1, m3p.adj) == g3
  
  
  @test m1 == m1p
  @test m2 == m2p
  @test m3 == m2p # yes it is m2, due to the special heuristic of metis
end

@testset "metis_dual" begin
  if !haskey(ENV, "METIS_LIB") 
    @error "you must define METIS_LIB in your variable environment"
  end

  eptr = Int32[0, 3, 6, 9, 12, 15, 18]
  eind = Int32[0, 1, 2, 0, 1, 5, 1, 5, 4, 1, 4, 6, 1, 6, 3, 1, 3, 2]
  xadj, adjncy = metis_mesh_to_dual(ne=Int32(6), nn= Int32(7), eptr=eptr, eind=eind, ncommon = Int32(1), baseval=Int32(0))
  
  #eptr = Int32[ 1, 4, 7, 10, 13, 16, 19]
  #eind = Int32[ 0, 1, 2, 3, 1, 2, 6, 2, 6, 5, 2, 5, 7, 2, 7, 4, 2, 4, 3]
  #xadj1, adjncy1 = metis_mesh_to_dual(ne=Int32(6), nn= Int32(7), eptr=eptr, eind=eind, ncommon = Int32(1), baseval=Int32(1))

  @test true
end


@testset "dgraph_dual" begin
  using MPI
  if !haskey(ENV, "METIS_LIB") | !haskey(ENV, "PARMETIS_LIB")
    @error "you must define PARMETIS_LIB and METIS_LIB in your variable environment"
  end
  mpiexec() do cmd
    testdir = @__DIR__
    run(`$cmd -n 3 $(Base.julia_cmd()) $(joinpath(testdir, "test_parmetis.jl"))`)
    @test true
  end
end
