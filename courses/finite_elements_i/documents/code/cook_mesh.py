"""
Creates th geometry for the cook problem.
"""
from dolfin import *
import numpy as np

def cook_mesh_random(size=100.0):
	# Create the empty mesh
	coarse_mesh = Mesh() 
	editor = MeshEditor()
	editor.open(coarse_mesh, 2, 2) 
	editor.init_vertices(13)  
	editor.init_cells( 16)
	editor.add_vertex(0, np.array([12, 20.25]))
	editor.add_vertex(1, np.array([0.0, 0.0]))
	editor.add_vertex(2, np.array([24.0, 22.0]))
	editor.add_vertex(3, np.array([24.0, 37.0]))     
	editor.add_vertex(4, np.array([0.0, 22.0]))
	editor.add_vertex(5, np.array([12.0,38.75 ]))
	editor.add_vertex(6, np.array([24.0, 52.0]))
	editor.add_vertex(7, np.array([0.0, 44.0]))
	editor.add_vertex(8, np.array([36.0, 50.25]))        
	editor.add_vertex(9, np.array([48.0, 52.0]))
	editor.add_vertex(10, np.array([48.0, 60.0]))
	editor.add_vertex(11, np.array([36.0, 38.75]))
	editor.add_vertex(12, np.array([48.0, 44.0]))
	editor.add_cell(0, np.array([0, 1, 2], dtype=np.uintp))
	editor.add_cell(1, np.array([0, 2, 3], dtype=np.uintp))
	editor.add_cell(2, np.array([0, 3, 4], dtype=np.uintp))
	editor.add_cell(3, np.array([0, 4, 1], dtype=np.uintp))
	editor.add_cell(4, np.array([5, 4, 3], dtype=np.uintp))
	editor.add_cell(5, np.array([5, 3, 6], dtype=np.uintp))
	editor.add_cell(6, np.array([5, 6, 7], dtype=np.uintp))
	editor.add_cell(7, np.array([5, 7, 4], dtype=np.uintp))
	editor.add_cell(8, np.array([8, 3, 9], dtype=np.uintp))
	editor.add_cell(9, np.array([8, 9, 10], dtype=np.uintp))
	editor.add_cell(10, np.array([8, 10, 6], dtype=np.uintp))
	editor.add_cell(11, np.array([8, 6, 3], dtype=np.uintp))
	editor.add_cell(12, np.array([11, 2, 12], dtype=np.uintp))
	editor.add_cell(13, np.array([11, 12, 9], dtype=np.uintp))
	editor.add_cell(14, np.array([11, 9, 3], dtype=np.uintp))
	editor.add_cell(15, np.array([11, 3, 2], dtype=np.uintp))
	editor.close()
	return coarse_mesh


def cook_mesh_right(size = 100.):
	# Creates a left biased triangular mesh for 
	# the cook problem
	# Create the empty mesh
	# The vertex coordinates of the mesh
	x_min = 0.
	x_max = 48.
	y_min_left = 0.
	y_min_right = 44.
	y_max_left = 44.
	y_max_right = 60.
	
	# Define the slopes of the line
	slope_bottom = (y_min_right-y_min_left)/(x_max-x_min)
	slope_top = (y_max_right-y_max_left)/(x_max-x_min)

	# Define some mesh parameters
	num_ver_x = 29 # Number of vertex along x
	num_ver_y = 29 # Number of vertex along y
	total_vertex = num_ver_x*num_ver_y # Total number of vertices

	num_cell_x = (num_ver_x-1)*2
	num_cell_y = num_ver_y-1

	num_cell = num_cell_x*num_cell_y

	coarse_mesh = Mesh() 
	editor = MeshEditor()
	editor.open(coarse_mesh, 2, 2) 
	editor.init_vertices(total_vertex)  
	editor.init_cells(num_cell)

	# Add the vertices
	for j in range(num_ver_y):
		for i in range(num_ver_x):
			editor.add_vertex(i+j*num_ver_x, np.array([x_min+i*((x_max-x_min)/(num_ver_x-1)), y_min_left+j*(y_max_left-y_min_left)/(num_cell_y)+i*(x_max-x_min)/(num_ver_x-1)*(slope_bottom-j*(slope_bottom-slope_top)/num_cell_y)]))
	i =0
	j=0
	# Add the cells
	for j in range(num_cell_y):
		for i in range(num_cell_x):
			editor.add_cell(i+num_cell_x*j, np.array([i-int(np.floor(i/(num_cell_x/2)))*(num_cell_x/2)+j*num_ver_y, (i+1)+int(np.floor(i/(num_cell_x/2)))+j*num_ver_x, i-int(np.floor(i/(num_cell_x/2)))*(num_cell_x/2)+j*num_ver_y+(num_ver_x+1)-int(np.floor(i)/(num_cell_x/2))], dtype=np.uintp))

	editor.close()
	return coarse_mesh

if __name__ == "__main__":
	mesh = cook_mesh_right()
	it = 0
	i = 0
	# Refine the mesh uniformly for a fixed number of times
	while i < it:
		cell_markers = CellFunction("bool", mesh)
		cell_markers.set_all(True)
		adapt(mesh,cell_markers)
		mesh=mesh.child()
		i+=1
	#plot(mesh, interactive=True)
	XDMFFile("mesh/cook_mesh.xdmf").write(mesh);
