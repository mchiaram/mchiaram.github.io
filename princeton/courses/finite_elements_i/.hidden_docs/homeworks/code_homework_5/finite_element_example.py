## @author Maurizio M. Chiaramonte CEE 513 @ Princeton University
#  @date 10/30/2018
#  @brief constructing the element 

import matplotlib 
import matplotlib.pyplot as plt
import numpy as np
from numpy import polynomial
from lagrange_polynomials import *
import sympy as sp
from sympy.plotting import plot


# @brief Compute the element stiffness matrix and forcing vector
# @param[in] element the element number assuming a numbering starting at zero
# @param[in] coordinate the array of coordinates 
# @param[in] connectivity the connectivity array which contains for element
#						the global nodal index (starting at zero) of the left and right node 
# @param[in] poly_order the degree of the polynomial interpolant
# @param[in,optional] f the function of the source term
# @return ke, fe the tuple of the local element stiffness and vector source
def element_stiffness(element, coordinates, connectivity, poly_order, f = 0 ):

	# The number of degrees of freedom
	num_dofs = poly_order + 1 

	# Create the local stiffness matrix
	ke = np.zeros( (num_dofs, num_dofs ) )

	# Create the local force vector
	fe = np.zeros( num_dofs )

	# Get the quadrature rule for the interval [-1,1]
	quadrature_order = int(np.ceil( (poly_order+1)/2 )+1)
	gauss_points, gauss_weights = np.polynomial.legendre.leggauss( quadrature_order  )

	# Get the coordinates of the nodes in the parametric domain
	# effectively this are equally spaced nodes in the interval [-1,1]
	nodes_crds_param = np.arange(-1., 1. + 1e-9,  2./(poly_order) )

	# Here we are simply interpolating the physical cordinates of the interior nodes
	# this could be more general but for the purpose of the homework it will suffice
	xi, xj = coordinates[ connectivity[ element ] ]
	nodes_crds_phys = xi + (nodes_crds_param + 1. )/2.*(xj - xi )

	# Fill in the stiffness matrix
	for q in range(len(gauss_points)):

		# Compute the map and jacobian of the map 
		x = 0
		jacobian = 0

		# Here we are using an isoparametric map
		for i in range(num_dofs):
			# Get the Gauss point position in physical space 
			x += lagrange_basis(nodes_crds_param,i,poly_order,gauss_points[q])*\
									nodes_crds_phys[i]
			# Get the jacobian of the mapping at the gauss point
			jacobian += d_lagrange_basis(nodes_crds_param,i,poly_order,gauss_points[q])*\
									nodes_crds_phys[i]

		# Sum over all degree of freedoms
		for i in range(num_dofs):
			
			for j in range(num_dofs):
				ke[i,j] += gauss_weights[q]*\
									d_lagrange_basis(nodes_crds_param,i,poly_order,gauss_points[q])*\
									d_lagrange_basis(nodes_crds_param,j,poly_order,gauss_points[q])*\
									jacobian**(-1.)

			# If we specified a source term
			if f:
				fe[i] -= gauss_weights[q]*\
								lagrange_basis(nodes_crds_param,i,poly_order,gauss_points[q])*\
								f( x )*jacobian

	return ke, fe

# @brief local to global map for arbitrary order finite elements. The numbering is
# done first for the exterior dofs and then for the internal dofs. The degrees of
# freedom are numbered starting from zero
# @param[in] poly_order the degree of the polynomial interpolant
# @param[in] element the element number assuming a numbering starting at zero
# @param[in] connectivity_array the connectivity array which contains for element
#						the global nodal index (starting at zero) of the left and right node 
# @ param[in] local_dof the local degree of freedom starting from zero going
#						in order from left to right
def local_to_global_map(poly_order, element, connectivity_array, local_dof):

	if not local_dof%(poly_order): 
		global_dof = connectivity_array[element][local_dof/poly_order]
		return global_dof
	
	# The number of elements in the mesh
	num_elements = len( connectivity_array )

	# The nodes per element (this are the external nodes)
	nodes_per_element = 2

	# Number of internal nodes is simply the total number of 
	# nodes minus the external nodes
	interal_nodes = poly_order + 1 - nodes_per_element

	# Get the global number of dof
	global_dof = num_elements + element*interal_nodes + local_dof

	return global_dof	

# @brief assembles the local arrays into the global arrays
# @param[in] poly_order the degree of the polynomial interpolant
# @param[in] element the element number assuming a numbering starting at zero
# @param[in] connectivity_array the connectivity array which contains for element
#						the global nodal index (starting at zero) of the left and right node 
# @param[in] ke the element stiffness matrix
# @param[in] fe the element stiffness matrix
# @param[out] K the reference to the global stiffness matrix
def assemble_local(poly_order, element, connectivity_array,ke, fe, K, F ):
	
	# The number of degrees of freedom
	num_dofs = poly_order + 1 

	# Loop over all degrees of freedom 
	for i in range(num_dofs):

		# Get the global index of the local i dof
		i_global = local_to_global_map(poly_order, element, connectivity_array, i)
		
		# Add the contribution to the source vector
		F[i_global] += fe[i]

		for j in range(num_dofs):

			# Get the global index of the local j dof
			j_global = local_to_global_map(poly_order, element, connectivity_array, j)

			# Add the contribution of the local stiffness to the global stiffness matrix		
			K[i_global, j_global] += ke[ i, j]

	return


# @brief assembles the local arrays into the global arrays
# @param[in] poly_order the degree of the polynomial interpolant
# @param[in] element the element number assuming a numbering starting at zero
# @param[in] connectivity_array the connectivity array which contains for element
#						the global nodal index (starting at zero) of the left and right node 
# @param[in] ke the element stiffness matrix
# @param[in] fe the element stiffness matrix
# @param[out] K the reference to the global stiffness matrix
# @param[out] F the reference to the global source vector
# @param[in,optional] f the function of the source term
def assemble_global(coordinates_array, connectivity_array, poly_order, K, F, f = 0 ):
	
	# The number of elements
	num_elements = len(connectivity_array)

	# Loop over all elements to assemble the global stiffness matrix
	for e in range(num_elements):

		# The element stiffness matrix and source vector
		ke, fe = element_stiffness(e, coordinates_array, connectivity_array, poly_order, f )
		
		# Assemble them 
		assemble_local(poly_order, e, connectivity_array, ke, fe, K, F )

	return


# @brief applies boundary conditions where appropriate
# @param[in] dofs the degree of freedom indeces of the boundary conditions
# @param[in] values thevalues of the boundary conditions
# @param[out] K the global stiffness matrix 
# @param[out] F the global forcing vector
def apply_bc(dofs, values, K, F ):

	# Loop over all degrees of freedom
	for i in range(len(dofs)):
		
		# Zero out the corresponding row 
		K[dofs[i],:] = 0 

		# Place one on the diagonal 
		K[dofs[i],dofs[i]] = 1 

		# Place the value of the boundary in the forcing vector
		F[dofs[i]] = values[i]
	return


# @brief compute the l2 norm of the error
# @param[in] coordinates the array of coordinates 
# @param[in] connectivity the connectivity array which contains for element
#						the global nodal index (starting at zero) of the left and right node 
# @param[in] poly_order the degree of the polynomial interpolant
# @param[in] uh the vector of degrees of freedom values
# @param[in] ue a function corresponding to the exact solution that can be called
# @return the l2 norm of the error
def compute_l2_norm_error(coordinates, connectivity, poly_order, uh, ue ):

	# The number of elements 
	num_elements = len(connectivity_array)

	# The number of degrees of freedom
	num_dofs = poly_order + 1 

	# Get the coordinates of the nodes in the parametric domain
	# effectively this are equally spaced nodes in the interval [-1,1]
	nodes_crds_param = np.arange(-1., 1. + 1e-9,  2./(poly_order) )

	# Get the quadrature rule for the interval [-1,1]
	quadrature_order = max( int(np.ceil( (poly_order+1)/2 )+1), 10 )
	gauss_points, gauss_weights = np.polynomial.legendre.leggauss( quadrature_order  )

	# Initialize the l2 norm of the error
	l2_err = 0

	# Perform the integration by summing the integrals ove all elements
	for e in range( num_elements ):
		
		# Here we are simply interpolating the physical cordinates of the interior nodes
		# this could be more general but for the purpose of the homework it will suffice
		xi, xj = coordinates[ connectivity[ e ] ]
		nodes_crds_phys = xi + (nodes_crds_param + 1. )/2.*(xj - xi )
	
		# Fill in the stiffness matrix
		for q in range(len(gauss_points)):
	
			# Compute the map and jacobian of the map 
			x_g = 0
			jacobian = 0
	
			# Here we are using an isoparametric map
			for i in range(num_dofs):
				# Get the Gauss point position in physical space 
				x_g += lagrange_basis(nodes_crds_param,i,poly_order,gauss_points[q])*\
										nodes_crds_phys[i]
				# Get the jacobian of the mapping at the gauss point
				jacobian += d_lagrange_basis(nodes_crds_param,i,poly_order,gauss_points[q])*\
										nodes_crds_phys[i]
				
			# Initialize the value of the approximate solution at the gauss point
			uh_g = 0 

			# Reconstruct the function at the gauss point
			for i in range(num_dofs):

				# Get the global index of the local i dof
				i_global = local_to_global_map(poly_order, e, connectivity_array, i)

				# Add the contribution of the i th dof to the value of the function
				uh_g += lagrange_basis(nodes_crds_param,i,poly_order,gauss_points[q])*uh[ i_global ] 
			
			# Add contribution of intergral at gauss point
			l2_err += pow( ue( x_g ) - uh_g , 2 )*jacobian*gauss_weights[q]

	return np.sqrt(l2_err)

if __name__ == "__main__":

	# Compute the boundary terms and source 
	# terms for a given solution 
	x = sp.symbols('x')
	ue = 1./2*sp.cos( 2.*sp.pi*x ) + sp.exp(x) + x**3
	f = sp.lambdify(x,sp.diff(ue,x,2))
	ue = sp.lambdify(x,ue)
	
	# The total number of elements 
	num_elements = pow(2,4) 

	# The polynomial degree 
	poly_order = 1

	# This are the nodes (only the external nodes) per element
	nodes_per_element = 2 

	# Fill in the connectivity array 
	connectivity_array = np.zeros((num_elements,nodes_per_element),dtype=int)
	for e in range(num_elements):
		connectivity_array[e,:] = np.arange(e,e+2)
	
	# The array of coordinates
	w = 1.
	coordinates_array = np.linspace(-w,w,num_elements+1)
	
	# The total number of degrees of freedom
	exterior_dofs_per_element = 2
	interior_dofs = num_elements*( poly_order + 1 - exterior_dofs_per_element )
	exterior_dofs = num_elements + 1
	total_dofs = interior_dofs + exterior_dofs

	# allocate memory for the global stiffness matrix
	K = np.zeros( (total_dofs, total_dofs) )

	# allocate memory for the global sourcevector
	F = np.zeros( total_dofs )

	# Assemble global stiffnes and source vector 
	assemble_global(coordinates_array, connectivity_array, poly_order, K, F, f  )

	# Apply boundary dirichlet bc at the left and right
	bc_vals = [ ue(-w), ue(w) ]
	bc_dofs = [ local_to_global_map(poly_order, 0, connectivity_array, 0),\
							local_to_global_map(poly_order, num_elements - 1, connectivity_array, poly_order ) ]

	# Apply the Dirichlet BC
	apply_bc(bc_dofs, bc_vals, K, F ) 

	# Solve finally
	u = np.linalg.solve( K, F ) 
 
	# Plot all the solution 
	x = np.linspace(-w,w,100)
	plt.plot( x , map( ue, x ) , 'b--')
	plt.plot( coordinates_array , u[:len(coordinates_array)] , 'ro-')
	plt.show()


