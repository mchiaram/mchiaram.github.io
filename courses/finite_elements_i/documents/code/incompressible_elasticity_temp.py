## @author CEE 513 @ Princeton University
#  @date 12/1/2017
#  @brief The following solves the cook's
# 	membrane problem to test the incompressible
#	eleasticity. A mixed formulation is implemeneted.

from __future__ import print_function
from dolfin import *
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt


if __name__ == '__main__':

	
	#---------------------------------------#
	# Define the problem parameters

	#Order of the polynomial interpolation
	order = 1

	# Space dimension in which we are solving the problem
	space_dim = 2

	# Create the mesh object
	mesh = Mesh()

	# We import the mesh
	XDMFFile("mesh/cook_mesh.xdmf").read(mesh)
	# The mesh uses the typical cook membrane geometry
	# ref. https://dealii.org/developer/doxygen/deal.II/code_gallery_Quasi_static_Finite_strain_Compressible_Elasticity.html

	# The material parameters
	nu = 0.49999999# poisson's ratio
	E = 1.0 # Young's modulus

	# Compute the lame` parameters
	mu = E/(2.0*(1.0+nu)) 
	lmbda = mu*nu/(0.5-nu)

	# We solve the problem using mixed formulation.
	# The strong and weak form of the problem could be found in 
	# Sec. 4.3 of Hughes book.

	# Define the function space
	# Define a vector element for displacement
	Ve_u = VectorElement("CG", mesh.ufl_cell(), order ,dim=space_dim)

	# Define a scalar element for pressure
	Ve_p = FiniteElement("CG", mesh.ufl_cell(), order)

	# Create the mixed function space
	V = FunctionSpace(mesh, Ve_u*Ve_p)

	# Specify the test and trial function
	(u,p) = TrialFunctions(V)
	(v,q) = TestFunctions(V)

	# Define the constitutive relation
	sigma = -p*Identity(space_dim)+ 2.0*mu*sym(grad(u)) 

	# Define the boundary condition

	# First describe the dirichlet boundary
	# The left edge is clammed
	class dirichlet_boundary(SubDomain): 
		def inside(self, x, on_boundary):
				return abs(x[0]) < DOLFIN_EPS and on_boundary #<--- Fill here

	# Describe the Neumann boundary
	class neumann_boundary(SubDomain):
		def inside(self, x, on_boundary):
				return abs(x[0]-48.0) < DOLFIN_EPS and on_boundary #<--- Fill here

	# Next we mark the boundaries for neumann bc
	boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
	dirichlet = dirichlet_boundary()
	neumann = neumann_boundary()
	dirichlet.mark(boundaries, 0)
	neumann.mark(boundaries, 1)
	ds = Measure("ds", subdomain_data=boundaries)

	# We specify the dirichlet boundary on the displacement
	# V.sub(0) applies bc on the first function space, displacement
	# in our case.
	bc = DirichletBC(V.sub(0), Constant((0.0,0.0)), dirichlet_boundary())

	# Define the variational form
	a = inner(sigma,grad(v))*dx - q*(p/lmbda + div(u))*dx #<--- Fill here

	# Define the forcing term
	t = Constant((0.0,1.0/16.))

	# Define the forcing functional
	f = Constant((0,0))
	F = dot(t,v)*ds(1) - dot(f , v)*dx

	# Compute the soultion
	sol = Function(V)
	solve(a==F, sol, bc)
	(uh, ph) = sol.split()

	# Add code to compute the maximum displacement
	V_u = FunctionSpace( mesh, Ve_u ) 
	uhu = interpolate(uh, V_u )
	max_displ = max(uhu.vector().array()[1::2] )
	print('Maximum vertical displacement', max_displ   )

	# Rename the functions
	uh.rename('displacement', 'displacement')
	ph.rename('pressure', 'pressure')

	file_u = XDMFFile("incompressible_elasticity/displacement.xdmf")
	file_u.write(uh,0)

	file_p = XDMFFile("incompressible_elasticity/pressure.xdmf")
	file_p.write(ph,0)
