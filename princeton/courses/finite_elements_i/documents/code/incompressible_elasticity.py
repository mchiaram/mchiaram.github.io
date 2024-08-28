
from dolfin import *
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":

	#Order of the polynomial interpolation
	order = 1

	# Space dimension in which we are solving the problem
	space_dim = 2

	# Create the mesh object
	lenght = 1. 
	height = lenght/10.
	cell_size = 0.025
	ndivx = int(lenght/cell_size)
	ndivy = int(height/cell_size)

	mesh = RectangleMesh(Point(0,0), Point(lenght, height), ndivx, ndivy, "crossed")

	# Define the function space
	poly_order_u = 2
	poly_order_p = 1

	# Define a vector element for displacement
	Ve_u = VectorElement("CG", mesh.ufl_cell(), poly_order_u ,dim=space_dim)

	# Define a scalar element for pressure
	Ve_p = FiniteElement("CG", mesh.ufl_cell(), poly_order_p)

	# Create the mixed function space
	V = FunctionSpace(mesh, Ve_u*Ve_p)

	# Specify the test and trial function
	(u,p) = TrialFunctions(V)
	(v,q) = TestFunctions(V)

	# Define the boundary condition

	# First describe the dirichlet boundary
	# The left edge is clammed
	class dirichlet_boundary(SubDomain):
		def inside(self, x, on_boundary):
				return abs(x[0]) < DOLFIN_EPS and on_boundary 
			
	# Describe the Neumann boundary
	class neumann_boundary(SubDomain):
		def inside(self, x, on_boundary):
				return abs(x[0]-lenght) < DOLFIN_EPS and on_boundary 

	# Next we mark the boundaries for neumann bc
	boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1 )
	dirichlet = dirichlet_boundary()
	boundaries.set_all(3)
	neumann = neumann_boundary()
	dirichlet.mark(boundaries, 0)
	neumann.mark(boundaries, 1)
	ds = Measure("ds", subdomain_data=boundaries)

	file_mesh = File('mesh.pvd')
	file_mesh << boundaries

	# We specify the dirichlet boundary on the displacement
	bc = DirichletBC(V.sub(0), Constant((0.0,0.0)), dirichlet_boundary())

	# Specify the test and trial function
	(u,p) = TrialFunctions(V)
	(v,q) = TestFunctions(V)

	# Define the traction term
	t = Constant((0.0,1.))
		
	# Define the body force term
	f = Constant((0.0,0.0))
	
	# Save the results
	file_u = File("displacement_incompressible_%i%i.pvd"%(poly_order_u,poly_order_p))
	file_p = File("pressure_incompressible_%i%i.pvd"%(poly_order_u,poly_order_p))

	# Define stress
	
	# The material parameters
	E = 1.0 # Young's modulus
	
	# Compute the lame` parameters
	nu = 0.5
	mu = Constant(E/(2.0*(1.0+nu)) )
	
	# Define the strain
	eps = sym(grad(u))
	
	# The weighted residual
	a = inner(2*mu*sym(grad(u)),grad(v))*dx - p*(div(v))*dx
	a += div(u)*q*dx
	
	# Define the forcing functional
	F = dot(t,v)*ds(1) -dot(f,v)*dx 

	# Compute the soultion
	uph = Function(V)
	
	solve(a == F, uph, bcs=bc)
	(uh, ph) = uph.split()

	# Rename the functions
	uh.rename('displacement', 'displacement')
	ph.rename('pressure', 'pressure')
	
	# Save the file
	file_u << uh
	file_p << ph
		
