# @param[in] nodes the coordinates of the nodes
# @param[in] index the index of the basis function
# @param[in] order the order of the polynomial basis
# @param[in] x the coordinate where to evaluate the basis
# @return the value of the index-th basis function at point x
def lagrange_basis(nodes, index, order, x): 
	
	# Check that we have the right number of nodes
	assert len(nodes) == (order + 1)

	# Create the initial value of the function
	ell = 1.

	# Loop over all nodes
	for i in range(0,order+1):
		
		# If this node is the same as the 
		# support node of the basis function, skip it
		if i == index:
			continue

		# Otherwise perform the multiplication
		ell *= ( x - nodes[i])/(nodes[index] - nodes[i]) # <- fill me here!!!
	
	return ell	

# @param[in] nodes the coordinates of the nodes
# @param[in] index the index of the basis function
# @param[in] order the order of the polynomial basis
# @param[in] x the coordinate where to evaluate the basis
# @return the value of the index-th basis function at point x
def d_lagrange_basis(nodes, index, order, x): 
	
	# Check that we have the right number of nodes
	assert len(nodes) == (order + 1)

	# Create the initial value of the function
	d_ell = 0

	# Loop over all nodes
	for i in range(0,order+1):
		
		# If this node is the same as the 
		# support node of the basis function, skip it
		if i == index:
			continue

		# Otherwise perform the multiplication
		d_ell += 1./(x - nodes[i]) # <- fill me here!!!
	
	# Multiply by the corresponding Lagrange basis
	d_ell *= lagrange_basis(nodes, index, order, x)
		
	return d_ell	

