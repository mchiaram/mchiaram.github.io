## @author CEE 513 @ Princeton University
#  @date 10/30/2018
#  @brief a simple example of how to use the implemented basis 

import matplotlib 
import matplotlib.pyplot as plt
import numpy as np
from numpy import polynomial
from lagrange_polynomials import *

# @brief formats a plot 
def format_plot():

	ax = plt.gca() 
	ax.set_xlabel(r'$x$')
	ax.set_ylabel(r'$f(x),f^p(x)$')
	ax.autoscale(enable=True, axis='x', tight=True)
	ax.set_xlim([-1.1,1.1])
	ax.set_ylim([-1,5])
	ax.grid('on')

# @param[in] f the function to be interpolated
# @param[in] nodes the nodes of the interpolation
# @param[in] poly_order the polynomial order of the interp
# @param[in] x the point where to evaluate the interpolation
# @return the value of the interpolation
def lagrange_interpolation( f, nodes, poly_order, x ):
	
	# Initialize the value of the interpolation
	fp_val = 0 

	# Sum up over all the basis 
	for i in range(0,poly_order+1):
		
		fp_val += ??? # <-- fill here

	return fp_val

if __name__ == "__main__":

	# The function we are trying to interpolate 
	f = lambda x: np.cos( np.pi*2*x )*0.5 + np.exp(x) + x**3
	
	# Some points for plotting 
	x = np.linspace(-1,1,200)

	# The order of polynomials to test
	min_order = 1
	max_order = 10

	# Test the interpolation for a orders
	# of polynomial interpolants
	for i in np.arange(min_order,max_order,1):
		
		# The polynomial order
		poly_order = i

		# Create equally spaced nodes
		nodes = np.arange(-1., 1. + 1e-9,  2./(poly_order) )

		# Create the approximation of our function

		# Evaluate exact function at nodes
		fp = lambda x: lagrange_interpolation( f, nodes, poly_order, x )

		# Plot the functions

		# The interpolated function
		plt.plot(x, fp(x),'--',c=(1-1./i,0,1./i),)
		plt.plot(nodes, fp(nodes),'o',c=(1-1./i,0,1./i),)

		# The exact function 
		plt.plot(x, f(x),'-k')

		# Make the plot pretty and show it
		format_plot()
		plt.draw()

		# Plot or show
		if False:
			plt.show()
		else:
			plt.savefig('lagrange_interp_%i.pdf'%poly_order)
			plt.close()
