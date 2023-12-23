from minimatrix import *
import numpy as np

data=[[-73, -97, -56, -25, -87],
                [ 44, -81, -89, -63,  95],
                [ 63,  95,  70,  45, -65],
                [-88,  69,  48,  10, -44],
                [ 77,   2,  52, -91,  30]]
a=Matrix(data)
b=Matrix([[5,7],[2,6]])
na=np.array(data)
print(f'{(a**4).dot(a)=}')
print(np.linalg.matrix_power(na,4)@na)