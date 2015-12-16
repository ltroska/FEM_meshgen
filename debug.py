from dolfin import *

mesh = Mesh("uniform.xml")
facets = MeshFunction('size_t',mesh,'uniform_facets.xml')
plot(mesh)
plot(facets)
interactive()

"""mesh = Mesh("/home/lukas/Downloads/lshape_fine.xml")
facets = MeshFunction('size_t',mesh,'/home/lukas/Downloads/lshape_fine_facet_region.xml')
plot(mesh)
plot(facets)
interactive()"""