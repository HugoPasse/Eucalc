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
print('')

cplx.preproc_radon_transform()

direction = np.array([2,2])
ect = cplx.compute_euler_characteristic_transform(direction)

print([-2 + i for i in range(5)])
print(ect.vectorize(-2,2,5))