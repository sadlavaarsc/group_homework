# Test code for IEEE course final project
# Fan Cheng, 2024

import minimatrix as mm

#Test your code here!
# The following code is only for your reference
# Please write the test code yourself

if __name__ == "__main__":
    print("Testing the Matrix library")

    # 创建一个3x3的矩阵A
    data_A = [[1, 2, 3], [6, 5, 4], [7, 8, 9]]
    A = mm.Matrix(data_A)

    # 创建一个3x3的矩阵B
    data_B = [[5, 2, 4], [5, 2, 5], [0, 3, 4]]
    B = mm.Matrix(data_B)

    # 显示A和B的形状
    print("Shape of A:", A.shape())
    print("Shape of B:", B.shape())

    # 计算A的转置
    print("Transpose of A:")
    print(A.T())
    # 计算A的方幂
    print("A ** cls5:")
    print(A**5)

    # 计算A和B的乘积
    C = A.dot(B)
    print("A dot B:")
    print(C)

    # 计算A和B的和
    D = A + B
    print("A + B:")
    print(D)

    # 计算A和B的差
    E = A - B
    print("A - B:")
    print(E)

    # 计算A和B的按元素乘积
    E = A * B
    print("A * B:")
    print(E)

    # 计算A的行列式
    det_A = A.det()
    print("Determinant of A:", det_A)

    # 计算A的逆矩阵
    inv_A = A.inverse()
    print("Inverse of A:")
    print(inv_A)

    # 计算A的秩
    rank_A = A.rank()
    print("Rank of A:", rank_A)

    # 使用arange()函数生成0到24的矩阵
    m24 = mm.arange(0, 24)
    print("Original matrix m24:")
    print(m24)

    print("\nReshaping m24 to [3,8]:")
    print(m24.reshape([3, 8]))

    print("\nReshaping m24 to [24,1]:")
    print(m24.reshape([24, 1]))

    print("\nReshaping m24 to [4,6]:")
    print(m24.reshape([4, 6]))

    print("\n3. Testing zeros() and zeros_like() functions")

    zero_mat = mm.zeros([3, 3])
    print("Zero matrix:")
    print(zero_mat)

    zero_like_m24 = mm.zeros_like(m24)
    print("\nZero matrix like m24:")
    print(zero_like_m24)

    print("\n4. Testing ones() and ones_like() functions")

    ones_mat = mm.ones([3, 3])
    print("Ones matrix:")
    print(ones_mat)

    ones_like_m24 = mm.ones_like(m24)
    print("\nOnes matrix like m24:")
    print(ones_like_m24)

    print("\n5. Testing nrandom() and nrandom_like() functions")

    random_mat = mm.nrandom([3, 3])
    print("Random matrix:")
    print(random_mat)

    random_like_m24 = mm.nrandom_like(m24)
    print("\nRandom matrix like m24:")
    print(random_like_m24)

    print("\n6. Testing least squares problem using custom Matrix class")

    m = 1000
    n = 100

    # 生成随机矩阵和向量
    X = mm.nrandom((m, n))
    w = mm.nrandom((n,1))
    e = mm.nrandom((m, 1))
    # 零均值化
    average=e.sum()/m
    e=e-mm.Matrix(dim=(m,1),init_value=average)
    
    #算Y
    Y=X.dot(w)+e 
    # 最小二乘估计
    gw=(X.T().dot(X)).inverse().dot(X.T()).dot(Y)
    gY=X.dot(gw)+e 

    print("Estimated w and actual w:")
    print(f'{gw=}\n{w=}')
    print("Estimated Y and actual y:")
    print(f'{gY=}\n{Y=}')

    # 比较结果
    diff=gw-w
    print("\nDifference between estimated w and actual w:")
    print(f'{diff}')
    print("\nAverage difference:")
    print(f'{diff.sum()/m}')