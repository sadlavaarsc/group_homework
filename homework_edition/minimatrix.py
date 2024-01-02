# Framework for IEEE course final project
# Fan Cheng, 2022
# Finished by LWT
# 正文（派生类）在136行
import random
import math

class Matrix_base(object):
    #tol控制的是某些运算（高斯消元)以及转字符串的时候的精度
    #经过单元测试，虽然不知道是否是numpy的问题，但在纯float的精度下四位小数是多次运算后的可靠性的保证
    def __init__(self, dim, init_value=0,tol=1e-4):
        self.data = [[init_value for _ in range(
            dim[1])] for _ in range(dim[0])]
        self.dim = dim
        self.tol=tol
	# 用来辅助输出的宽度获取
    def get_column_widths(self):
        widths = []
        for j in range(self.dim[1]):
            col_width = max(len(format_float_with_max_precision(self.data[i][j], int(abs(math.log10(self.tol)))))
                            for i in range(self.dim[0]))
            widths.append(col_width)
        return widths

    def __repr__(self):
        col_widths = self.get_column_widths()
        res = '['
        for i in range(self.dim[0]):
            if i != 0:
                res += ' '
            res += '['
            for j in range(self.dim[1]):
                if j != 0:
                    res += ' '
                res += format_float_with_max_precision(self.data[i][j], int(
                    abs(math.log10(self.tol)))).rjust(col_widths[j])
            res += ']' if i == self.dim[0] - 1 else ']\n'
        return res + ']'

    def add(self, mat):
        if self.dim != mat.dim:
            raise ValueError("Matrix dimensions do not match for addition.")
        result = Matrix_base(self.dim)
        for i in range(self.dim[0]):
            for j in range(self.dim[1]):
                result.data[i][j] = self.data[i][j] + mat.data[i][j]
        return result

    # 数乘
    def kmul(self, k):
        result = Matrix_base((self.dim[0], self.dim[1]))
        for i in range(self.dim[0]):
            for j in range(self.dim[1]):
                result.data[i][j] = self.data[i][j] * k
        return result

    def switch_rows(self, i, j):
        self.data[i], self.data[j] = self.data[j], self.data[i]

    def k_mul_row(self, i, k):
        self.data[i] = [elem * k for elem in self.data[i]]

    def row_plus_krow(self, i, j, k):
        self.data[i] = [elem1 + k * elem2 for elem1,
                        elem2 in zip(self.data[i], self.data[j])]

    def get_main_element_index(self, i):
        for j in range(self.dim[1]):
            if self.data[i][j] != 0:
                return j
        return -1

    def get_not_zero_row_num(self):
        count = 0
        for i in range(self.dim[0]):
            if any(self.data[i]):
                count += 1
        return count
    def format_row(self,id):
        for i in range(self.dim[0]):
            if abs(self.data[id][i])<self.tol:
                self.data[id][i]=0
    # 高斯约旦消元，不是纯高斯消元是为了后面进行一些操作（比如求行列式，判断维数
    # 会改变原来的矩阵，Normalized会尝试将主元归一（不归一是用来求行列式的）
    # tol参数控制容差，于0差距小于tol的数字都会被当作0
    # 函数返回一个值，表示行列式正负（每交换一次行行列式就会变号）
    def gauss(self,normalized=True):
        # 与手工消元略有不同，这里进行了微调
        # 粗消元
        res=1
        for i in range(self.dim[0]):
            main_elem_index = self.get_main_element_index(i)
            if main_elem_index == -1 or self.data[i][main_elem_index] == 0:
                continue

            # 将主元归一
            main_elem = self.data[i][main_elem_index]
            if normalized:
                self.k_mul_row(i, 1 / main_elem)
                #print(f'row {i} / {main_elem} result:{self}')
                main_elem=1
            for j in range(i + 1, self.dim[0]):
                if self.data[j][main_elem_index] != 0:
                    k = -self.data[j][main_elem_index]/main_elem
                    self.row_plus_krow(j, i, k)
                    self.format_row(j)
                    #print(f'row {j} + {k} * row {i} result:{self}')
        # 调整行顺序，这里选择重载一下比较函数，并追踪交换行的次数
        # 这里选择在每一行的一位开头额外增加一列储存主元位置
        swap_time=0
        def swap_key(row):
            nonlocal swap_time
            swap_time+=1
            return row[0]
        sort_data=sorted([[self.get_main_element_index(i)]+self.data[i] for i in range(self.dim[0])],key=swap_key)
        self.data=[e[1:] for e in sort_data]
        # 处理主元所在列
        for i in range(self.dim[0]-1,-1,-1):
            meId=self.get_main_element_index(i)
            for j in range(i-1,-1,-1):
                #print(f'{j=} {meId=} {self.data[j][meId]=}')
                #self.data[j][meId]=0
                k = -self.data[j][meId]/self.data[i][meId]
                self.row_plus_krow(j,i,k)
                self.format_row(j)
        return swap_time

    def copy(self):
        new_matrix = Matrix_base(self.dim)
        for i in range(self.dim[0]):
            for j in range(self.dim[1]):
                new_matrix.data[i][j] = self.data[i][j]
        return new_matrix

# 给str用的，截取小数点位数
def format_float_with_max_precision(number, max_precision):
    number_str = str(number)
    decimal_index = number_str.find('.')
    
    if decimal_index == -1:
        return number_str
    decimal_part = number_str[decimal_index + 1:]
    
    if len(decimal_part) > max_precision:
        return f"{number:.{max_precision}f}"
    else:
        return number_str

class Matrix(Matrix_base):
	# tol默认为1e-8，有必要可以自行更改
	def __init__(self, data=None, dim=None, init_value=0):
		if data != None:
			self.data=[e.copy() for e in data]
			self.dim=(len(data),len(data[0]))
			self.tol=1e-8
		elif dim!=None:
			super().__init__(dim,init_value)
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
			t_data.append(t_lst[i*newdim[1]:(i+1)*newdim[1]])
		self.dim=newdim
		return Matrix(t_data)

	def dot(self, mat):
		if self.dim[1] != mat.dim[0]:
			raise ValueError(
                "Number of columns in the first matrix must be equal to the number of rows in the second matrix.")
		result = Matrix(dim=(self.dim[0], mat.dim[1]))
		for i in range(self.dim[0]):
			for j in range(mat.dim[1]):
				for k in range(self.dim[1]):
					result.data[i][j] += self.data[i][k] * mat.data[k][j]
		return result

	# 转置
	def T(self):
		n_data=[]
		for j in range(self.dim[1]):
			t_data=[]
			for i in range(self.dim[0]):
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
				return sum([sum(self.data[i]) for i in range(self.dim[0])])
		

	def Kronecker_product(self, other):
		d0,d1=other.dim[0],other.dim[1]
		res=Matrix(dim=(self.dim[0]*d0,self.dim[1]*d1),init_value=0)
		for ai in range(self.dim[0]):
			for aj in range(self.dim[1]):
				for bi in range(d0):
					for bj in range(d1):
						res.data[ai*d0+bi][aj*d1+bj]=self.data[ai][aj]*other.data[bi][bj]
		return res
	
	def __getitem__(self, key):
		if not isinstance(key,tuple) or len(key)!=2:
			raise ValueError('The index should be a tuple whose lenghth is 2')
		#这里支持了step第三个切片
		start,stop,step=[0,0],[0,0],[0,0]
		def format_slice(id):
			nonlocal start,stop,step
			start[id]=key[id].start if key[id].start!=None else 0
			stop[id]=key[id].stop if key[id].stop!=None and key[id].stop<=self.dim[id] else self.dim[id]
			step[id]=key[id].step if key[id].step!=None else 1
		if isinstance(key[0],slice) and isinstance(key[1],slice):
			format_slice(0)
			format_slice(1)
			#print(f'{start=} {stop=} {step=}')
			t_lst=[]
			for i in range(start[0],stop[0],step[0]):
				tt_lst=[]
				for j in range(start[1],stop[1],step[1]):
					tt_lst.append(self.data[i][j])
				t_lst.append(tt_lst)
			return Matrix(t_lst)
		else:
			return self.data[key[0]][key[1]]

	def __setitem__(self, key, value):
		if not isinstance(key,turple) and len(key)!=2:
			raise ValueError('The index should be a turple whose lenghth is 2')
		#这里支持了step第三个切片
		start,stop,step=[0,0],[0,0],[0,0]
		def format_slice(id):
			nonlocal start,stop,step
			start[id]=key[id].start if key[id].start!=None else 0
			stop[id]=key[id].stop if key[id].stop!=None and key[id].stop<=self.dim[id] else self.dim[id]
			step[id]=key[id].step if key[id].step!=None else 1
		if isinstance(key[0],slice) and isinstance(key[1],slice):
			format_slice(0)
			format_slice(1)
			# 这里支持了把某一块全部设为一个值
			t_mat=None
			if isinstance(value,Matrix):
				t_mat=value
			else:
				t_mat=Matrix(dim=self.dim,init_value=value)
			for i in range(start[0],stop[0],step[0]):
				for j in range(start[1],stop[1],step[1]):
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
		res=I(self.dim[0])
		t_m=self.copy()
		while n>0:
			if n&1:
				res=res.dot(t_m)
			t_m=t_m.dot(t_m)
			n>>=1
		return res

	def __add__(self, other):
		return Matrix_base_to_Matrix(self.add(other))

	def __sub__(self, other):
		return Matrix_base_to_Matrix(self.add(other.kmul(-1)))

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

	# 基类中已经实现了__repr__，不过这里还是重写一个类似的功能
	def __str__(self):
		col_widths = self.get_column_widths()
		res = '['
		for i in range(self.dim[0]):
			if i != 0:
				res += ' '
			res += '['
			for j in range(self.dim[1]):
				if j != 0:
					res += ' '
				res += format_float_with_max_precision(self.data[i][j], int(abs(math.log10(self.tol)))).rjust(col_widths[j])
			res += ']' if i == self.dim[0] - 1 else ']\n'
		return res + ']'

	def det(self):
		# 原理：初等变换不改变行列式的值
		# 这里直接对行和列分别高斯消元得到一个对角阵，再对对角阵进行对角元累乘
		if self.dim[0]!=self.dim[1]:
			raise ValueError('The Matrix should be a square matrix')
		# 两套消元
		t_m=self.copy()
		st=t_m.gauss(False)
		t_m=t_m.T()
		st+=t_m.gauss(False)
		# 累乘
		res=1
		for i in range(self.dim[0]):
			res*=t_m.data[i][i]
		return res*(-1)**st

	def inverse(self):
		# 用了经典的拼接原矩阵和单位矩阵进行高斯消元
		if self.dim[0]!=self.dim[1]:
			raise ValueError('The matrix should be a square matrix')
		eye=I(self.dim[0])
		t_m=concatenate((self,eye),0)
		t_m.gauss()
		# 先检验是否是奇异阵
		if t_m.get_not_zero_row_num()!=t_m.dim[0]:
			raise ValueError('The matrix is a singular matrix')
		else:
			return t_m[:,self.dim[1]:self.dim[1]*2]
		
	#求矩阵秩
	def rank(self):
		t=self.copy()
		t.gauss()
		return t.get_not_zero_row_num()

def Matrix_base_to_Matrix(matbase):
	return Matrix(matbase.data)
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
	res.data=[list(range(start,end,step))]
	return res

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
		ite=iter(items)
		res=next(ite).copy()
		for mat in ite:
			if mat.dim[0]!=res.dim[0]:
				raise ValueError('The matrixs to be oncatenated should have proper size')
			for i in range(mat.dim[0]):
				res.data[i].extend(mat.data[i])
		res.dim=(len(res.data),len(res.data[0]))
		return res
	else:
		ite=iter(items)
		res=next(ite).copy()
		for mat in ite:
			if mat.dim[1]!=res.dim[1]:
				raise ValueError('The matrixs to be oncatenated should have proper size')
			for e in mat.data:
				res.append(e)
		res.dim=(len(res.data),len(res.data[0]))
		return res
			
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
	#det(A)
	pass