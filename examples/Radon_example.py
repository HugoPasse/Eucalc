import numpy as np
import matplotlib.pyplot as plt
import sys
 
sys.path.append('../')

import embedded_cubical_complex as ecc


data = np.array([[1,0,0,0],
				[0,1,1,0],
				[1,1,1,1],
				[1,1,0,0],])

cplx = ecc.EmbeddedComplex(data,input_top_cells=False)
cplx.impose_upper_star_filtration_from_vertices()
cplx.print_filtration()

print('\n---------------------------------------------')
print('Pre-processing')
cplx.preproc_radon_transform()

directions = [[-1,-1],[-1,1],[1,-1],[1,1]]
for d in directions:
	i = cplx.get_vector_index(d)
	print('\nDirection :',d)
	print('Ordinary critical vertices  :',cplx.get_ordinary_critical_vertices(i))
	print('Ordinary critical values    :',cplx.get_ordinary_critical_values(i))
	print('Classical critical vertices :',cplx.get_classical_critical_vertices(i))
	print('Classical critical values   :',cplx.get_classical_critical_values(i))
print('')

direction = np.array([2,2])
radon = cplx.compute_radon_transform(direction)
attributes = radon.get_attributes()

print(attributes)

print('T  = ',["{0:0.3f}".format(a) for a in attributes[0]])
print('E  = ',["{0:0.3f}".format(a) for a in attributes[1]])

print('Tc = ',["{0:0.3f}".format(a) for a in attributes[2]])
print('Ec = ',["{0:0.3f}".format(a) for a in attributes[3]])



for i in range(len(attributes[0])-1):
	plt.hlines(attributes[1][i],attributes[0][i],attributes[0][i+1])
plt.scatter(attributes[2],attributes[3])

plt.title('Radon transform')
plt.ylabel(r'$R(\xi,t)$')
plt.xlabel(r't')
plt.show()