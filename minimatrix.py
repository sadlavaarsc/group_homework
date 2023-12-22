# Framework for IEEE course final project
# Fan Cheng, 2022

import random
from matrix_base import Matrix_base

class Matrix:
	def __init__(self, data=None, dim=None, init_value=0):
		if data != None:
			self.data={e.copy() for e in data}
			self.dim=(len(data),len(data[0]))
		elif dim!=None:
			super('Matrix',self).__init__(dim,init_value)
		else:
			raise ValueError('Wrong init paramter')

	# 重载copy防止copy出基类
	def copy(self):
		cp=Matrix(self.data)
		return cp

	def shape(self):
		return self.dim

	def reshape(self, newdim):
		if newdim[0]*newdim[1]!=self.dim[0]*self.dim[1]:
			raise ValueError('The new dimission isn\'t equal to the old.')
		# 先接一串，再分割
		t_lst=[]
		for e in self.data:
			t_lst.extend(e)
		t_data=[]
		for i in range(newdim[0]):
			t_data.append(t_data[i*newdim[1]:(i+1)*newdim[1]])
		return Matrix(t_data)

	'''基类中已经实现
	def dot(self, other):
		pass'''
	# 转置
	def T(self):
		n_data=[]
		for j in range(res.dim[1]):
			t_data=[]
			for i in range(res.dim[0]):
				t_data.append(self.data[i][j])
			n_data.append(t_data)
		return Matrix(n_data)

	def sum(self, axis=None): 
		#当然可以只实现一种然后让另一种转置，这里为了性能没有这么做
		match(axis):
			case 0:
				t_data=[]
				for j in range(self.dim[1]):
					s=0
					for i in range(self.dim[0]):
						s+=self.data[i][j]
					t_data.append(s)
				return Matrix(t_data)
			case 1:
				return Matrix([[sum(self.data[i])] for i in range(self.dim[0])])
			case _:
				#没有定义异常值的处理，这里选择全部捕获
				return Matrix([sum([sum(self.data[i]) for i in range(self.dim[0])])])
		

	def Kronecker_product(self, other):
		d0,d1=other.dim[0],other.dim[1]
		res=Matrix(dim=(self.dim[0]*d0,self.dim[1]*d1),init_value=0)
		for ai in range(self.dim[0]):
			for aj in range(self.dim[1]):
				for bi in range(d0):
					for bj in range(d1):
						res[ai*d0+bi][aj*d1+bj]=self.data[ai][aj]*other.data[bi][bj]
		return res
	
	def __getitem__(self, key):
		if not isinstance(key,turple) and len(key)!=2:
			raise ValueError('The index should be a turple whose lenghth is 2')
		#这里支持了step第三个切片
		def format_slice(id):
			if key[id].start==None:
				key[id].start=0
			if key[id].stop==None:
				key[id].stop=self.dim[id]
			if key[id].step==None:
				key[id].step=1
		if isinstance(key[0],slice) and isinstance(key[1],slice):
			format_slice(0)
			format_slice(1)
			t_lst=0
			for i in range(key[0].start,key[0].stop,key[0].step):
				tt_lst=[]
				for j in range(key[1].start,key[1].stop,key[1].step):
					tt_lst.append(self.data[i][j])
				t_lst.append(tt_lst)
			return Matrix(t_lst)
		else:
			return self.data[key[0]][key[1]]

	def __setitem__(self, key, value):
		if not isinstance(key,turple) and len(key)!=2:
			raise ValueError('The index should be a turple whose lenghth is 2')
		#这里支持了step第三个切片
		def format_slice(id):
			if key[id].start==None:
				key[id].start=0
			if key[id].stop==None:
				key[id].stop=self.dim[id]
			if key[id].step==None:
				key[id].step=1
		if isinstance(key[0],slice) and isinstance(key[1],slice):
			format_slice(0)
			format_slice(1)
			# 这里支持了把某一块全部设为一个值
			t_mat=None
			if isinstance(value,Matrix):
				t_mat=value
			else:
				t_mat=Matrix(dim=self.dim,init_value=value)
			for i in range(key[0].start,key[0].stop,key[0].step):
				for j in range(key[1].start,key[1].stop,key[1].step):
					self.data[i][j]=t_mat[i][j]
		else:
			if isinstance(value,Matrix):
				raise ValueError('The right value should be a number')
			self.data[key[0]][key[1]]=value

	def __pow__(self, n):
		if n<0 or int(n)!=n:
			raise ValueError('n should be a natural number')
		if self.dim[0]!=self.dim[1]:
			raise ValueError('The matrix should be a square matrix')
		# 快速幂
		res=Matrix.I(self.dim[0])
		t_m=self.copy()
		while n>0:
			if n&1:
				res=res.dot(t_m)
			t_m.dot(t_m)
			n>>=1
		return res

	def __add__(self, other):
		return self.add(other)

	def __sub__(self, other):
		return self.add(other.kmul(-1))

	def __mul__(self, other):
		if self.dim!=other.dim:
			raise ValueError('Dimisions of the two matrix must be the same.')
		res=Matrix(dim=(self.dim[0],self.dim[1]))
		for i in range(self.dim[0]):
			for j in range(self.dim[1]):
				res.data[i][j] = self.data[i][j]*other.data[i][j]
		return res

	def __len__(self):
		return self.dim[0]*self.dim[1]

	def __str__(self):
		r"""
		按照
		[[  0   1   4   9  16  25  36  49]
 		 [ 64  81 100 121 144 169 196 225]
 		 [256 289 324 361 400 441 484 529]]
 		的格式将矩阵表示为一个 字符串
 		！！！ 注意返回值是字符串
		"""
		res='['
		for i in range(self.dim[0]):
			res+='['
			for j in range(self.dim[1]):
				res+=f'\t{self.data[i][j]}'
			res+=']' if i==self.dim[0]-1 else']\n'
		return res

	def det(self):
		# 原理：初等变换不改变行列式的值
		# 这里直接对行和列分别高斯消元得到一个对角阵，再对对角阵进行对角元累乘
		if self.dim[0]!=self.dim[1]:
			raise ValueError('The Matrix should be a square matrix')
		# 两套消元
		t_m=self.copy()
		st=t_m.guass()
		t_m=t_m.T()
		st+=t_m.guass()
		# 累乘
		res=1
		for i in range(self.dim[0]):
			res*=t_n[i][i]
		return res*(-1)**st

	def inverse(self):
		# 用了经典的拼接原矩阵和单位矩阵进行高斯消元
		r"""
		计算非奇异方阵的逆矩阵。对于非方阵或奇异阵的情形应抛出异常。
		要求: 该函数应不改变 self 的内容; 该函数的时间复杂度应该不超过 O(n**3).
		提示: Gauss消元

		Returns:
			Matrix: 一个 Matrix 实例，表示逆矩阵
		"""
		if self.dim[0]!=self.dim[1]:
			raise ValueError('The matrix should be a square matrix')
		eye=self.I(self.dim[0])
		t_m=self.concatenate((self,eye),0)
		t_m.guass()
		# 先检验是否是奇异阵
		if t_m.get_not_zero_row_num()!=t_m.dim[0]:
			raise ValueError('The matrix is a singular matrix')
		else:
			return t_m[self.dim[0]:self.dim[0]*2][self.dim[1]:self.dim[1]*2]
		
	#求矩阵秩
	def rank(self):
		t=self.copy()
		t.gauss()
		return t.get_not_zero_row_num()

def I(n):
	res=Matrix(dim=(n,n))
	for i in range(n):
		res.data[i][i]=1
	return res

def narray(dim, init_value=1):
	return Matrix(dim=dim,init_value=init_value)

def arange(start, end, step=1):
	n=(end-start)//step
	res=Matrix(dim=(1,n))
	res.data=list(range(start,end,step))
	pass

def zeros(dim):
	return Matrix(dim=dim,init_value=0)

def zeros_like(matrix):
	return Matrix(dim=matrix.dim,init_value=0)

def ones(dim):
	return Matrix(dim=dim,init_value=1)

def ones_like(matrix):
	return Matrix(dim=matrix.dim,init_value=1)

def nrandom(dim):
	res=Matrix(dim=dim)
	for i in range(dim[0]):
		for j in range(dim[1]):
			res.data[i][j]=random.random()
	return res

def nrandom_like(matrix):
	dim=matrix.dim
	res=Matrix(dim=dim)
	for i in range(dim[0]):
		for j in range(dim[1]):
			res.data[i][j]=random.random()
	return res

def concatenate(items, axis=0):
	pass
	if axis==0:
		res=next(items).copy()
		for mat in items:
			if mat.dim[0]!=res.dim[0]:
				raise ValueError('The matrixs to be oncatenated should have proper size')
			for i in range(mat.dim[0]):
				res.data[i].extend(mat.data[i])
	else:
		res=next(items).copy()
		for mat in items:
			if mat.dim[1]!=res.dim[1]:
				raise ValueError('The matrixs to be oncatenated should have proper size')
			for e in mat.data:
				res.append(e)
			
#其实就是让实现一个类似装饰器的效果
#将给定函数进行向量化
def vectorize(func):
	def inner(mat):
		for i in range(mat.dim[0]):
			for j in range(mat.dim[1]):
				mat.data[i][j]=func(mat.data[i][j])
	return inner


if __name__ == "__main__":
	print("test here")
	data=[[5,2,4],[5,2,5],[0,3,4]];
	A=Matrix(data)
	det(A)
	pass