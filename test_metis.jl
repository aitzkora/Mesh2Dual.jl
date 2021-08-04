using Mesh2Dual

eptr = Int32[0, 4, 8, 12, 16, 20, 24, 28, 32]
eind = Int32[0, 1, 5, 6, 1, 2, 6, 7, 2, 3, 7, 8, 3, 4, 8, 9, 5, 6, 
             10, 11, 6, 7, 11, 12, 7, 8, 12, 13, 8, 9, 13, 14   ]
xadj, adjcy = metis_mesh_to_dual(ne = 8, nn = 15,  eptr=eptr, eind=eind, baseval=0, ncommon=2)
println("xadj = $xadj \nadjcy = $adjcy")
