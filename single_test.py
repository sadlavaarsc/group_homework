from minimatrix import *

a=Matrix([[3,4],[9,3]])
b=Matrix([[5,7],[2,6]])
print('a=',a)
print('b=',b)
print('a+b=',a+b)
print('a-b=',a-b)
print('a*b=',a.dot(b))
print('rank(a)=',a.rank())
print('det(a)=',a.det())
print('inv(a)=',a.inverse())