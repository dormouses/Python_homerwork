# -*- coding: utf-8 -*-

class Error(Exception):
    def __init__(self, _value):
        self.value = _value
    def __str__(self):
        return self.value

class Gauss:
    def __init__(self, _a, _b):
        n=len(_b)
        n1=len(_a)
        if (n!=n1):
            raise Error('Размеры вектора и матрицы несовместимы')
        self.a=_a
        for i in range (n):
            self.a[i].append(-_b[i])
        self.x = [0] * n
        #прямой ход
        for i in range(n):
            tmp=self.a[i][i]
            for j in range(n, i-1, -1):
                self.a[i][j]/=tmp
            for j in range(i+1, n):
                tmp=self.a[j][i]
                for k in range(n, i-1, -1):
                    self.a[j][k]-=tmp*self.a[i][k]

        #обратный ход
        self.x[n-1]=self.a[n-1][n]
        for i in range(n-2, -1, -1):
            self.x[i]=self.a[i][n]
            for j in range (i+1, n):
                self.x[i]-=self.a[i][j]*self.x[j]

    def __getitem__(self, key):
        if (key<0) or (key >=len(self.a)):
            raise Error('Недопустимое значение')
        return self.x[key]
    
    def __setitem__(self, key, value):
         raise Error('Нельзя изменить ответ')
    

try:
    n = int (raw_input('Введите количество уравнений '))
    a = [None] * n;
    for i in range(n):
        a[i] = list(map(float, raw_input().split()))
    m = int (raw_input('Ведите длину вектора '))
    b = list(map(float, raw_input().split()))
 
    ans =  Gauss (a, b)
    
    for i in range (len(b)):
        print ans[i]
    
    ans[0]=10

except Error as ex:
    print 'ой'
    print ex