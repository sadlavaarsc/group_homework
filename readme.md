### 项目说明

编码：LWT

调试：LWT

(我也知道这个是团队作业，不过一边写一边找人人没找到代码写完了...)

项目包含文件与说明：

| 文件名                   | 说明                                                         |
| ------------------------ | ------------------------------------------------------------ |
| readme.md                | 本文件，对项目进行了必要说明                                 |
| 技术报告.pdf/技术报告.md | 技术报告，更加正式地表述了此项目的技术细节                   |
| minimatrix.py            | 库本体，实现了Matrix_base（基类）和Matrix（派生类）          |
| main.py                  | 作业要求的测试&展示文件，输出较长，有必要可以尝试用类似于python ./main.py > log.txt 这种方式存到文件中仔细观看（备注：输出到文件中时矩阵的右对齐可能会因为默认字体设置等等各种原因导致没有完全对齐，具体效果以直接在控制台输出的版本为准） |
| unit_test.py（附件）     | 个人测试时写的单元测试程序，对大部分函数（包括Kronecker_product）都进行了大规模单元测试并使用numpy进行对拍，需要请自取（需要pytest） |
| single_test.py（附件）   | 个人测试的时候用来调试一些特定出错数据的文件，目前状态为调试numpy在数据范围过大时溢出的问题，需要测试特定功能可以直接在此处编辑 |

[项目链接](https://github.com/sadlavaarsc/group_homework)

详细项目报告见``项目报告.pdf``

备注：

1. 考虑到精度误差问题，库中的高斯消元采用了容差对0进行判断，同时默认输出时不输出全部小数点位数（默认只输出四位小数，容差计算同理，可以在构造函数的tol处修改）
2. 题目在矩阵转文字部分格式没有太多说明，考虑到输出结果长度不一，先获取最宽列元素再统一向右对齐的手法
3. 原本要求的接口是``__str__``但是``__str__``在使用fstring输出（例如``print(f'{mat=}')``）时并不会正确输出，这里额外添加了``__repr__``
4. 由于代码的前半部分是用平板在在线ide上编写的，当时的ide采用了空格缩进，但后半段用电脑写的代码却是tab缩进，导致直接复制黏贴两个类的东西会出问题，得重新调整缩进，暂时没想到办法批量调整（悲）
4. 以下的文档是当初初步设计时的草案，虽然与实际代码略有不同，但已经基本反映了设计思路，如有必要可以参考具体代码，这里再复述一遍大体设计思路：由于各个函数相互依赖较为严重耦合度高，这里先实现一个Matrix_base类实现了最基本的矩阵行变换、高斯消元以及数据储存等等需要被后面复杂函数反复调用的功能，随后其派生出实现具体业务逻辑的耦合度相对低的Matrix类，在Matrix类中实现了各种具体功能，并通过其他架构对代码进行了单元测试。具体的架构并不是特别复杂，整个库没有用到什么特别的先进计算架构或者管线，这里就简单描述一下代码组成部分了。

### 需求明确

1. 基本需求：在既定接口上实现二维矩阵，并且完成测试

2. 具体需求：

```python

class Matrix:
	r"""
	自定义的二维矩阵类

	Args:
		data: 一个二维的嵌套列表，表示矩阵的数据。即 data[i][j] 表示矩阵第 i+1 行第 j+1 列处的元素。
			  当参数 data 不为 None 时，应根据参数 data 确定矩阵的形状。默认值: None
		dim: 一个元组 (n, m) 表示矩阵是 n 行 m 列, 当参数 data 为 None 时，根据该参数确定矩阵的形状；
			 当参数 data 不为 None 时，忽略该参数。如果 data 和 dim 同时为 None, 应抛出异常。默认值: None
		init_value: 当提供的 data 参数为 None 时，使用该 init_value 初始化一个 n 行 m 列的矩阵，
					即矩阵各元素均为 init_value. 当参数 data 不为 None 时，忽略该参数。 默认值: 0
	
	Attributes:
		dim: 一个元组 (n, m) 表示矩阵的形状
		data: 一个二维的嵌套列表，表示矩阵的数据

	Examples:
		>>> mat1 = Matrix(dim=(2, 3), init_value=0)
		>>> print(mat1)
		>>> [[0 0 0]
			 [0 0 0]]
		>>> mat2 = Matrix(data=[[0, 1], [1, 2], [2, 3]])
		>>> print(mat2)
		>>> [[0 1]
			 [1 2]
			 [2 3]]
	"""
	def __init__(self, data=None, dim=None, init_value=0):
		# self.data
		# self.dim
		pass

	def shape(self):
		r"""
		返回矩阵的形状 dim
		"""
		pass

	def reshape(self, newdim):
		r"""
		将矩阵从(m,n)维拉伸为newdim=(m1,n1)
		该函数不改变 self
		
		Args:
			newdim: 一个元组 (m1, n1) 表示拉伸后的矩阵形状。如果 m1 * n1 不等于 self.dim[0] * self.dim[1],
					应抛出异常
		
		Returns:
			Matrix: 一个 Matrix 类型的返回结果, 表示 reshape 得到的结果
		"""
		pass

	def dot(self, other):
		r"""
		矩阵乘法：矩阵乘以矩阵
		按照公式 A[i, j] = \sum_k B[i, k] * C[k, j] 计算 A = B.dot(C)

		Args:
			other: 参与运算的另一个 Matrix 实例
		
		Returns:
			Matrix: 计算结果
		
		Examples:
			>>> A = Matrix(data=[[1, 2], [3, 4]])
			>>> A.dot(A)
			>>> [[ 7 10]
				 [15 22]]
		"""
		pass

	def T(self):
		r"""
		矩阵的转置

		Returns:
			Matrix: 矩阵的转置

		Examples:
			>>> A = Matrix(data=[[1, 2], [3, 4]])
			>>> A.T()
			>>> [[1 3]
				 [2 4]]
			>>> B = Matrix(data=[[1, 2, 3], [4, 5, 6]])
			>>> B.T()
			>>> [[1 4]
				 [2 5]
				 [3 6]]
		"""
		pass 

	def sum(self, axis=None): 
		r"""
		根据指定的坐标轴对矩阵元素进行求和

		Args:
			axis: 一个整数，或者 None. 默认值: None
				  axis = 0 表示对矩阵进行按列求和，得到形状为 (1, self.dim[1]) 的矩阵
				  axis = 1 表示对矩阵进行按行求和，得到形状为 (self.dim[0], 1) 的矩阵
				  axis = None 表示对矩阵全部元素进行求和，得到形状为 (1, 1) 的矩阵
		
		Returns:
			Matrix: 一个 Matrix 类的实例，表示求和结果

		Examples:
			>>> A = Matrix(data=[[1, 2, 3], [4, 5, 6]])
			>>> A.sum()
			>>> [[21]]
			>>> A.sum(axis=0)
			>>> [[5 7 9]]
			>>> A.sum(axis=1)
			>>> [[6]
				 [15]]
		"""
		pass

	def copy(self):
		r"""
		返回matrix的一个备份

		Returns:
			Matrix: 一个self的备份
		"""
		pass

	def Kronecker_product(self, other):
		r"""
		计算两个矩阵的Kronecker积，具体定义可以搜索，https://baike.baidu.com/item/克罗内克积/6282573

		Args:
			other: 参与运算的另一个 Matrix

		Returns:
			Matrix: Kronecke product 的计算结果
		"""
		pass
	
	def __getitem__(self, key):
		r"""
		实现 Matrix 的索引功能，即 Matrix 实例可以通过 [] 获取矩阵中的元素（或子矩阵）

		x[key] 具备以下基本特性：
		1. 单值索引
			x[a, b] 返回 Matrix 实例 x 的第 a 行, 第 b 列处的元素 (从 0 开始编号)
		2. 矩阵切片
			x[a:b, c:d] 返回 Matrix 实例 x 的一个由 第 a, a+1, ..., b-1 行, 第 c, c+1, ..., d-1 列元素构成的子矩阵
			特别地, 需要支持省略切片左(右)端点参数的写法, 如 x 是一个 n 行 m 列矩阵, 那么
			x[:b, c:] 的语义等价于 x[0:b, c:m]
			x[:, :] 的语义等价于 x[0:n, 0:m]

		Args:
			key: 一个元组，表示索引

		Returns:
			索引结果，单个元素或者矩阵切片

		Examples:
			>>> x = Matrix(data=[
						[0, 1, 2, 3],
						[4, 5, 6, 7],
						[8, 9, 0, 1]
					])
			>>> x[1, 2]
			>>> 6
			>>> x[0:2, 1:4]
			>>> [[1 2 3]
				 [5 6 7]]
			>>> x[:, :2]
			>>> [[0 1]
				 [4 5]
				 [8 9]]
		"""
		pass

	def __setitem__(self, key, value):
		r"""
		实现 Matrix 的赋值功能, 通过 x[key] = value 进行赋值的功能

		类似于 __getitem__ , 需要具备以下基本特性:
		1. 单元素赋值
			x[a, b] = k 的含义为，将 Matrix 实例 x 的 第 a 行, 第 b 处的元素赋值为 k (从 0 开始编号)
		2. 对矩阵切片赋值
			x[a:b, c:d] = value 其中 value 是一个 (b-a)行(d-c)列的 Matrix 实例
			含义为, 将由 Matrix 实例 x 的第 a, a+1, ..., b-1 行, 第 c, c+1, ..., d-1 列元素构成的子矩阵 赋值为 value 矩阵
			即 子矩阵的 (i, j) 位置赋值为 value[i, j]
			同样地, 这里也需要支持如 x[:b, c:] = value, x[:, :] = value 等省略写法
		
		Args:
			key: 一个元组，表示索引
			value: 赋值运算的右值，即要赋的值

		Examples:
			>>> x = Matrix(data=[
						[0, 1, 2, 3],
						[4, 5, 6, 7],
						[8, 9, 0, 1]
					])
			>>> x[1, 2] = 0
			>>> x
			>>> [[0 1 2 3]
				 [4 5 0 7]
				 [8 9 0 1]]
			>>> x[1:, 2:] = Matrix(data=[[1, 2], [3, 4]])
			>>> x
			>>> [[0 1 2 3]
				 [4 5 1 2]
				 [8 9 3 4]]
		"""
		pass

	def __pow__(self, n):
		r"""
		矩阵的n次幂，n为自然数
		该函数应当不改变 self 的内容

		Args:
			n: int, 自然数

		Returns:
			Matrix: 运算结果
		"""
		pass

	def __add__(self, other):
		r"""
		两个矩阵相加
		该函数应当不改变 self 和 other 的内容

		Args:
			other: 一个 Matrix 实例
		
		Returns:
			Matrix: 运算结果
		"""
		pass

	def __sub__(self, other):
		r"""
		两个矩阵相减
		该函数应当不改变 self 和 other 的内容

		Args:
			other: 一个 Matrix 实例
		
		Returns:
			Matrix: 运算结果
		"""
		pass

	def __mul__(self, other):
		r"""
		两个矩阵 对应位置 元素  相乘
		注意 不是矩阵乘法dot
		该函数应当不改变 self 和 other 的内容

		Args:
			other: 一个 Matrix 实例
		
		Returns:
			Matrix: 运算结果

		Examples:
			>>> Matrix(data=[[1, 2]]) * Matrix(data=[[3, 4]])
			>>> [[3 8]]
		"""
		pass


	def __len__(self):
		r"""
		返回矩阵元素的数目

		Returns:
			int: 元素数目，即 行数 * 列数
		"""
		pass

	def __str__(self):
		r"""
		按照
		[[  0   1   4   9  16  25  36  49]
 		 [ 64  81 100 121 144 169 196 225]
 		 [256 289 324 361 400 441 484 529]]
 		的格式将矩阵表示为一个 字符串
 		！！！ 注意返回值是字符串
		"""
		pass

	def det(self):
		r"""
		计算方阵的行列式。对于非方阵的情形应抛出异常。
		要求: 该函数应不改变 self 的内容; 该函数的时间复杂度应该不超过 O(n**3).
		提示: Gauss消元
		
		Returns:
			一个 Python int 或者 float, 表示计算结果
		"""
		pass

	def inverse(self):
		r"""
		计算非奇异方阵的逆矩阵。对于非方阵或奇异阵的情形应抛出异常。
		要求: 该函数应不改变 self 的内容; 该函数的时间复杂度应该不超过 O(n**3).
		提示: Gauss消元

		Returns:
			Matrix: 一个 Matrix 实例，表示逆矩阵
		"""
		pass

	def rank(self):
		r"""
		计算矩阵的秩
		要求: 该函数应不改变 self 的内容; 该函数的时间复杂度应该不超过 O(n**3).
		提示: Gauss消元

		Returns:
			一个 Python int 表示计算结果
		"""
		pass

def I(n):
	'''
	return an n*n unit matrix
	'''

def narray(dim, init_value=1): # dim (,,,,,), init为矩阵元素初始值
	r"""
	返回一个matrix，维数为dim，初始值为init_value
	
	Args:
		dim: Tuple[int, int] 表示矩阵形状
		init_value: 表示初始值，默认值: 1

	Returns:
		Matrix: 一个 Matrix 类型的实例
	"""
	#return Matrix(dim, None, init_value)

def arange(start, end, step):
	r"""
	返回一个1*n 的 narray 其中的元素类同 range(start, end, step)

	Args:
		start: 起始点(包含)
		end: 终止点(不包含)
		step: 步长

	Returns:
		Matrix: 一个 Matrix 实例
	"""
	pass

def zeros(dim):
	r"""
	返回一个维数为dim 的全0 narray

	Args:
		dim: Tuple[int, int] 表示矩阵形状

	Returns:
		Matrix: 一个 Matrix 类型的实例
	"""
	pass

def zeros_like(matrix):
	r"""
	返回一个形状和matrix一样 的全0 narray

	Args:
		matrix: 一个 Matrix 实例
	
	Returns:
		Matrix: 一个 Matrix 类型的实例

	Examples:
		>>> A = Matrix(data=[[1, 2, 3], [2, 3, 4]])
		>>> zeros_like(A)
		>>> [[0 0 0]
			 [0 0 0]]
	"""
	pass

def ones(dim):
	r"""
	返回一个维数为dim 的全1 narray
	类同 zeros
	"""
	pass

def ones_like(matrix):
	r"""
	返回一个维数和matrix一样 的全1 narray
	类同 zeros_like
	"""
	pass

def nrandom(dim):
	r"""
	返回一个维数为dim 的随机 narray
	参数与返回值类型同 zeros
	"""
	pass

def nrandom_like(matrix):
	r"""
	返回一个维数和matrix一样 的随机 narray
	参数与返回值类型同 zeros_like
	"""
	pass

def concatenate(items, axis=0):
	r"""
	将若干矩阵按照指定的方向拼接起来
	若给定的输入在形状上不对应，应抛出异常
	该函数应当不改变 items 中的元素

	Args:
		items: 一个可迭代的对象，其中的元素为 Matrix 类型。
		axis: 一个取值为 0 或 1 的整数，表示拼接方向，默认值 0.
			  0 表示在第0维即行上进行拼接
			  1 表示在第1维即列上进行拼接
	
	Returns:
		Matrix: 一个 Matrix 类型的拼接结果

	Examples:
		>>> A, B = Matrix([[0, 1, 2]]), Matrix([[3, 4, 5]])
		>>> concatenate((A, B))
		>>> [[0 1 2]
			 [3 4 5]]
		>>> concatenate((A, B, A), axis=1)
		>>> [[0 1 2 3 4 5 0 1 2]]
	"""
	pass

def vectorize(func):
	r"""
	将给定函数进行向量化
	
	Args:
		func: 一个Python函数
	
	Returns:
		一个向量化的函数 F: Matrix -> Matrix, 它的参数是一个 Matrix 实例 x, 返回值也是一个 Matrix 实例；
		它将函数 func 作用在 参数 x 的每一个元素上
	
	Examples:
		>>> def func(x):
				return x ** 2
		>>> F = vectorize(func)
		>>> x = Matrix([[1, 2, 3],[2, 3, 1]])
		>>> F(x)
		>>> [[1 4 9]
			 [4 9 1]]
		>>> 
		>>> @vectorize
		>>> def my_abs(x):
				if x < 0:
					return -x
				else:
					return x
		>>> y = Matrix([[-1, 1], [2, -2]])
		>>> my_abs(y)
		>>> [[1, 1]
			 [2, 2]]
	"""
	pass
```

从上面可以发现，绝大多数函数对其他函数没有依赖，存在依赖链主要有：

```python
1. __pow__->__mul__
2. det->gauss->初等行变换
3. inverse->gauss->初等行变换
4. rank->gauss->初等行变换
```

### 架构设计

备注：以下设计为初期设计的草稿，仅供参考，具体接口以代码为准

主要实现一个基本矩阵类Matrix_base，实现如下基本功能：

```python
class Matrix_base(object):
    # 与子类相同
    def __init__(self, data=None, dim=None, init_value=0):
        # 题目提到的属性
        self.data=data.copy()# 储存数据
        self.dim=dim#维数
        pass
    # 加乘
   	def add(self,mat):
        pass
    # 这个是严格的矩阵乘法
    def dot(self,mat):
        pass
    # 初等行变换
    # 交换行
    def switch_rows(self,i,j):
        pass
    # 某行 × k倍
   	def k_mul_row(self,i,k):
        pass
    # 第i行+第j行*k倍
    def row_plus_krow(self,i,j,k):
        pass
    # 获取第i行主元位置
    def get_main_element_index(self,i):
        pass
    # 获取非零行数
    def get_not_zero_row_num(self):
        pass
    # 高斯消元，会修改self
    def gauss(self):
        pass
    # copy
    def copy(self):
        pass
```

在此基础上实现一个Matrix_base的派生类Matrix，分工如下：

```python

class Matrix:
# A
	def __init__(self, data=None, dim=None, init_value=0):
		# self.data
		# self.dim
		pass
# B
	def shape(self):
		pass
# B
	def reshape(self, newdim):
		pass
# B
	def dot(self, other):
		pass
# B
	def T(self):
		pass 
# B
	def sum(self, axis=None): 
		pass
# A
	def copy(self):
		pass
# B
	def Kronecker_product(self, other):
		pass
	# A
	def __getitem__(self, key):
		pass
# A
	def __setitem__(self, key, value):
		pass
# A
	def __pow__(self, n):
		pass
# B
	def __add__(self, other):
		pass
# B
	def __sub__(self, other):
		pass
# B
	def __mul__(self, other):
		pass

# B
	def __len__(self):
		pass
# B
	def __str__(self):
		pass
# A
	def det(self):
		pass
# A 
	def inverse(self):
		pass
# A
	def rank(self):
		pass
# B
def I(n):
    pass
# B
def narray(dim, init_value=1):
    pass
# B
def arange(start, end, step):
	pass
# B
def zeros(dim):
	pass
# B
def zeros_like(matrix):
	pass
# B
def ones(dim):
	pass
# B
def ones_like(matrix):
	pass
# B
def nrandom(dim):
	pass
# B
def nrandom_like(matrix):
	pass
# A
def concatenate(items, axis=0):
	pass
# A
def vectorize(func):
	pass
```

A，B是描述当初准备找第二个人合作时的分工，很遗憾没来得及找到
