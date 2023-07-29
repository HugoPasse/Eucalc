import numpy as np
import matplotlib.pyplot as plt
import sys
 
sys.path.append('../')

import embedded_cubical_complex as ecc

def kernel(x):
	return np.exp(x)

data = np.array([[1,0,0,0],
				[0,1,1,0],
				[1,1,1,1],
				[1,1,0,0],])

cplx = ecc.EmbeddedComplex(data,input_top_cells=False)	
cplx.impose_upper_star_filtration_from_vertices()
cplx.print_filtration()
print('')

cplx.preproc_radon_transform()

direction = np.array([2,2])
index = cplx.get_vector_index(direction)

ordv = cplx.get_ordinary_critical_values(index)
ordp = cplx.get_ordinary_critical_vertices(index)

h = 0
for i in range(len(ordp)):
	tmp = kernel(np.dot(cplx.get_vertex_embedding(ordp[i]),direction))*ordv[i]
	print('Point',ordp[i],'with value',ordv[i],' : ',end='')
	print('h =',h,'+',tmp)
	h -= tmp

print('HT(' + str(direction) + ')=',h)
print('HT :',cplx.compute_hybrid_transform('exp',[direction]))