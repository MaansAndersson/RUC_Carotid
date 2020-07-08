from dolfin import *
import math
import time
import os
import sys
import os

import numpy as np

vtk = File('u.pvd')
#set_log_level(1)
lumen = Mesh('XMLs/19-03-22--12-44Lumen.xml')
wall = Mesh('XMLs/19-03-22--12-47dilatedWall.xml')

lumen.scale(1/165)
wall.scale(1/165)
lumen.translate(Point(0.25, 0.25, 0))
wall.translate(Point(0.25, 0.25, 0))


domain = UnitCubeMesh(7, 7, 2)
#domain.scale(3.0)

for i in range(0,3):
    cell_markers = MeshFunction("bool", domain, domain.topology().dim())
    for cell in cells(domain):
      cell_markers[cell] = False
      p = cell.midpoint()
      if (abs(p.y()-0.5) < 0.2 and abs(p.x()-0.5) < 0.2):
          cell_markers[cell] = True
    domain = refine(domain, cell_markers)
print(len(domain.cells()))



V = FunctionSpace(domain, "Lagrange", 1)
Z = FunctionSpace(domain, "DG", 0)


lumenMeshFunction = MeshFunction('size_t',lumen, lumen.topology().dim())
wallMeshFunction = MeshFunction('size_t',wall, wall.topology().dim())
"""
nu = Function(V)

v = TestFunction(V)
u = TrialFunction(V)

ds1 = Measure('ds', domain=domain, subdomain_data=lumenMeshFunction)
ds2 = Measure('ds', domain=domain, subdomain_data=wallMeshFunction)

n1 = FacetNormal(lumen)
n2 = FacetNormal(wall)
#x = Expression(('x[0]','x[1]','x[2]'),degree = 1)

#a = u*v*dx(mesh)
"""
"""b = dot(x,n2)*ds1(mesh) # - dot(x,n1)*ds2
assemble(b)
solve(a,nu,b)
"""
vtk << (lumenMeshFunction)
vtk << (wallMeshFunction)
#vtk << (nu)

#nL = FacetNormal(lumen)
#nW = FacetNormal(wall)

lumenbound = lumen.bounding_box_tree() #BoundingBoxTree() #A.build(lumen)
wallbound = wall.bounding_box_tree() #BoundingBoxTree() #B.build(wall)

class InitialCondition(UserExpression):
    #def __init__(self,lumenbound,wallbound):
    #    self.lumenbound = lumenbound
    #    self.wallbound = wallbound
    #print(lumen.coordinates()[a_intersect[0],:])
    def eval(self, values, x):
        p = Point(x[0], x[1], x[2])
        a_intersect_entity = (lumenbound.compute_closest_entity(p))
        b_intersect_entity = (wallbound.compute_closest_entity(p))
        
        """
        print(lumen.num_cells())
        print(a_intersect_entity)
        print(wall.num_cells())
        print(b_intersect_entity)
      
        print(len(lumen.coordinates()[:,0]))
        print(len(wall.coordinates()[:,0]))
        
        
        print("hej: ", len(lumen.cells()))
        
        print(lumen.cells()[a_intersect_entity[0]])
        
        print(lumen.coordinates()[lumen.cells()[a_intersect_entity[0]][0]])
        print(lumen.coordinates()[lumen.cells()[a_intersect_entity[0]][1]])
        print(lumen.coordinates()[lumen.cells()[a_intersect_entity[0]][2]])
        
        
        print(wall.cells()[b_intersect_entity[0]])
        
        print(wall.coordinates()[wall.cells()[b_intersect_entity[0]][0]])
        print(wall.coordinates()[wall.cells()[b_intersect_entity[0]][1]])
        print(wall.coordinates()[wall.cells()[b_intersect_entity[0]][2]])
        """
        a_intersect_point = lumen.coordinates()[lumen.cells()[a_intersect_entity[0]][0]]
        b_intersect_point = wall.coordinates()[wall.cells()[b_intersect_entity[0]][0]]
        
        if True: # or (b_intersect_entity[0] < wall.num_cells() and a_intersect_entity[0] < lumen.num_cells()):

            
            a_intersect_dir = a_intersect_point - x
            b_intersect_dir = b_intersect_point - x
            if np.dot(a_intersect_dir,b_intersect_dir) < 0. :
                #values[0] = 1
            #if a_intersect_entity[1] < 0.001: #Temoprary def around Lumen
                values[0] = 1;
            else:
                values[0] = 0;
            
            #print(type(lumen.coordinates()[b_intersect_entity[0],:]))
            #print(a_intersect_entity)
            #print(b_intersect_entity)
            #values[0] = 1.
        else :
            values[0] = 0
        return values
    
    #def value_shape(self):
     #   return(1,)


        
II = InitialCondition()

vtk << project(II,Z)
        
"""
cell_markers = MeshFunction("bool", mesh, mesh.topology().dim())
for c in cells(mesh):
    cell_markers[c] = False
    p = c.midpoint[c]
   if(i < 3):
        if(c.midpoint().distance(p0) < 0.75 and c.circumradius() > 0.05):
            cell_markers[c] = True

    if((c.midpoint().distance(p1) < 0.2 or c.midpoint().distance(p2) < 0.2) and c.circumradius() > 0.03):
        cell_markers[c] = True
"""
#nu.vector()[:] = 0



