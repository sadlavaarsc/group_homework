class Matrix_base(object):
    def __init__(self, dim, init_value=0):
        self.data = [[init_value for _ in range(
            dim[1])] for _ in range(dim[0])]
        self.dim = dim

    def add(self, mat):
        if self.dim != mat.dim:
            raise ValueError("Matrix dimensions do not match for addition.")
        result = Matrix_base(self.dim)
        for i in range(self.dim[0]):
            for j in range(self.dim[1]):
                result.data[i][j] = self.data[i][j] + mat.data[i][j]
        return result

    def dot(self, mat):
        if self.dim[1] != mat.dim[0]:
            raise ValueError(
                "Number of columns in the first matrix must be equal to the number of rows in the second matrix.")
        result = Matrix_base((self.dim[0], mat.dim[1]))
        for i in range(self.dim[0]):
            for j in range(mat.dim[1]):
                for k in range(self.dim[1]):
                    result.data[i][j] += self.data[i][k] * mat.data[k][j]
        return result
    # 数乘

    def kmul(self, k):
        result = Matrix_base((self.dim[0], mat.dim[1]))
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
    # 高斯约旦消元，不是纯高斯消元是为了后面进行一些操作（比如求行列式，判断维数
    # 会改变原来的矩阵，Normalized会尝试将主元归一（不归一是用来求行列式的）
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
                main_elem=1
                #print(f'row {i} / {main_elem} result:{self}')
            for j in range(i + 1, self.dim[0]):
                if self.data[j][main_elem_index] != 0:
                    k = -self.data[j][main_elem_index]/main_elem
                    self.row_plus_krow(j, i, k)
                    print(f'row {j} + {k} * row {i} result:{self}')
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
                self.data[j][meId]=0
        return swap_time

    def copy(self):
        new_matrix = Matrix_base(self.dim)
        for i in range(self.dim[0]):
            for j in range(self.dim[1]):
                new_matrix.data[i][j] = self.data[i][j]
        return new_matrix

    def __str__(self):
        return str(self.data)



a,b=Matrix_base((3,3)),Matrix_base((2,2))
a.data=[[5,2,4],[5,2,5],[0,3,4]]
b.data=[[2,3],[2,3]]
print(a.gauss(False))
print(a)
