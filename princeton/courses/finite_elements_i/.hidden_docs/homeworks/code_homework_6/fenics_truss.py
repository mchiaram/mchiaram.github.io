## @author CEE 513 @ Princeton University
#  @date 10/30/2018
#  @brief a simple example of how to use FEniCS 
#  by solving the Poisson problem in 1-D

from __future__ import print_function
from fenics import *
import sympy as sp
from sympy.plotting import plot
import numpy as np
import matplotlib.pyplot as plt


if __name__ == "__main__":

	# Compute the boundary terms and source 
	# terms for a given solution 
	x = sp.symbols('x[0]')
	solution = 1./2*sp.cos( 2.*sp.pi*x ) + sp.exp(x) + x**3
	ue_code = sp.ccode(solution).replace('M_PI','pi')
	t_code = sp.ccode( sp.diff(solution,x,1)).replace('M_PI','pi')
	f_code = sp.ccode( sp.diff(solution,x,2)).replace('M_PI','pi')
	solution = sp.lambdify(x,solution)
	
	# Start by creating a unit interval mesh
	# subdivided in ndiv elements
	ndiv = 10
	mesh = UnitIntervalMesh(ndiv)

	# We are going to shift the above mesh 
	# such that it dicretizes the interval [-1,1]
	mesh.coordinates()[:] *= 2
	mesh.coordinates()[:] -= 1
	
	# Create the function space of 
	poly_order = 1
	V = FunctionSpace(mesh, 'Lagrange', poly_order)
	
	# Define boundary condition
	
	# Define the dirichlet boundary 
	class dirichlet_boundary(SubDomain):
		def inside(self, x, on_boundary):
			return on_boundary and abs( x + 1. ) < DOLFIN_EPS

	# Define the neumann boundary 
	class neumann_boundary(SubDomain):
		def inside(self, x, on_boundary):
			return on_boundary and abs( x - 1. ) < DOLFIN_EPS

	# Define the boundary measure (effectively a way to 
	# distinguish between Neumann and Dirichlet)
	boundaries = FacetFunction("size_t", mesh)
	dirichlet_boundary().mark(boundaries,1)
	neumann_boundary().mark(boundaries,2)
	ds = ds(subdomain_data=boundaries)
	
	# Define the analaytical solution
	ue = Expression(ue_code, pi=np.pi, degree=5)

	# Define the dirichlet boundary conditions
	bc = DirichletBC(V, ue, dirichlet_boundary())
	
	# Define the trial and test function
	u = TrialFunction(V)
	v = TestFunction(V)
	
	# Define the forcing function
	f = Expression(f_code,pi=np.pi, degree=5)
	t = Expression(t_code,pi=np.pi, degree=5)
	
	# Define the bilinear form
	a = dot(u.dx(0), v.dx(0))*dx

	# Define the forcing 
	F = -f*v*dx + t*v*ds(2)
	
	# Compute solution
	uh = Function(V)
	solve(a == F, uh, bc)

	# Compute error in L2 norm
	error_L2 = errornorm(ue, uh, 'L2')
	
	# Compute maximum error at vertices
	vertex_values_ue = ue.compute_vertex_values(mesh)
	vertex_values_uh = uh.compute_vertex_values(mesh)
	
	error_max = np.max(np.abs(vertex_values_ue - vertex_values_uh))

	# Print errors
	print('error_L2  =', error_L2)
	print('error_max =', error_max)
	
	# Plot the values of the exact vs the computed 
	x = np.linspace(-1,1,100)
	plt.plot(x,solution(x),'r-')
	plt.plot(mesh.coordinates()[:],vertex_values_uh,'bo-')
	plt.show()
	

