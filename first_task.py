from models import Matrix

eq_cnt, var_cnt = [int(i) for i in input().split()]
tmp = []
data = []
for i in range(eq_cnt):
    tmp = [int(i) for i in input().split()]
    data.append(tmp)

matrix = Matrix(data=data)
print("A:")
print(str(matrix))
matrix.lup()
print("A after LUP:")
print(str(matrix))
print("L:")
print(str(matrix.get_l()))
print("U:")
print(str(matrix.get_u()))
print("P:")
print(str(matrix.p))
