import pytest
import numpy as np
import minimatrix as mm
import random
import time

#自动的多次测试装饰器
def multi_test(test_time):
    def out_wrapper(func):
        def wrapper(*argv,**kargv):
            print(f'{func.__name__} test start')
            for i in range(test_time):
                func(*argv,**kargv)
            print(f'{test_time} tests completed')
            return None
        return wrapper
    return out_wrapper
            

def timer(func):
    def wrapper(*argv,**kargv):
        t1=time.time()
        result=func(*argv,**kargv)
        t2=time.time()
        print(f'{func.__name__} cost {t2-t1}s')
        return result
    return wrapper
# 判断两种不同类型的矩阵是否相等，这里采用容差检测解决精度问题
def is_equal(a:np.array,b:mm.Matrix,tol=1e-4):
    return np.all(np.abs(a - b.data) < tol)
# 生成指定范围内的整数随机矩阵，返回一个tuple
# 第一个元素为minimatrix，第二个为numpy.array
def rand_mat(dim=(5,5),rge=(-100,100)):
    data=[]
    for i in range(dim[0]):
        row=[]
        for j in range(dim[1]):
            row.append(random.randint(rge[0],rge[1]))
        data.append(row)
    mmat=mm.Matrix(data)
    nmat=np.array(data)
    return mmat,nmat

# 接下来就是正经测试
@multi_test(test_time=10000)
def test_plus():
    ma,na=rand_mat()
    mb,nb=rand_mat()
    assert is_equal(na+nb,ma+mb)

@multi_test(test_time=10000)
def test_sub():
    ma,na=rand_mat()
    mb,nb=rand_mat()
    assert is_equal(na-nb,ma-mb)

@multi_test(test_time=10000)
def test_mul():
    ma,na=rand_mat()
    mb,nb=rand_mat()
    assert is_equal(na*nb,ma*mb)

@multi_test(test_time=10000)
def test_dot():
    ma,na=rand_mat()
    mb,nb=rand_mat()
    assert is_equal(na@nb,ma.dot(mb))

@multi_test(test_time=10000)
def test_T():
    ma,na=rand_mat()
    assert is_equal(na.T,ma.T())

@multi_test(test_time=10000)
def test_inv():
    ma,na=rand_mat()
    assert is_equal(np.linalg.inv(na),ma.inverse()),f'{np.linalg.inv(na)=}\n{ma.inverse().data=}'

@multi_test(test_time=10000)
def test_rank():
    ma,na=rand_mat()
    assert np.linalg.matrix_rank(na)==ma.rank()

@multi_test(test_time=10000)
def test_kron():
    ma,na=rand_mat()
    mb,nb=rand_mat()
    assert is_equal(np.kron(na,nb),ma.Kronecker_product(mb))
@multi_test(test_time=10000)
def test_sum():
    ma,na=rand_mat()
    assert is_equal(np.sum(na,axis=1),ma.sum(1).T()) \
    and is_equal(np.sum(na,axis=0),na.sum(0)) ,\
    f'{np.sum(na,axis=1)=}\n{ma.sum(1).T().data=}'
@multi_test(test_time=10000)
def test_pow():
    # 太大的数据Numpy会溢出，这里就算了
    ma,na=rand_mat(rge=(-10,10))
    n=random.randint(0,5)
    assert is_equal(np.linalg.matrix_power(na,n),ma**n),\
        f'{na=}\n{np.linalg.matrix_power(na,n)=}\n{(ma**n).data=}\n{n=}'