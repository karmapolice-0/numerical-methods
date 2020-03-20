import copy
from numbers import Number


class Matrix:
    EPS = 10e-6
    ROUND_DIG = 13

    def __init__(self, data):
        self.data = data
        self.rows = len(self.data)
        self.columns = len(self.data[0])
        self.p = []
        self.init_p()

    def __deepcopy__(self, memodict={}):
        if memodict is None:
            memodict = {}
        new = Matrix(copy.deepcopy(self.data))
        return new

    def __len__(self):
        return self.rows

    def __getitem__(self, index):
        return self.data[index]

    def __setitem__(self, key, value):
        self.data[key] = value

    def __str__(self):
        tmp = copy.deepcopy(self)
        for i in range(len(tmp)):
            for j in range(len(tmp[i])):
                if abs(tmp[i][j]) <= self.EPS:
                    tmp[i][j] = 0
        s = [[str(e) for e in row] for row in tmp]
        lens = [max(map(len, col)) for col in zip(*s)]
        fmt = '\t'.join('{{:{}.{}}}'.format(x + 2, x if x > self.ROUND_DIG else self.ROUND_DIG) for x in lens)
        table = [fmt.format(*row) for row in s]
        return '\n'.join(table)

    def __round__(self, n=14):
        for i in range(self.rows):
            for j in range(self.columns):
                if abs(self.data[i][j]) < self.EPS:
                    self.data[i][j] = 0
                else:
                    self.data[i][j] = round(self.data[i][j], n)

    def __eq__(self, o: object) -> bool:
        if isinstance(o, Matrix):
            if len(self.data) != len(o.data) or len(self.data[0]) != len(o.data[0]):
                return False
            else:
                return self.data == o.data
        return False

    def __add__(self, other):
        if not isinstance(other, Matrix):
            raise AttributeError("Unsupported operand type for +: '{}' and '{}'".format(self.__class__, type(other)))
        new_data = []
        if len(self.data) != len(other.data) or len(self.data[0]) != len(other.data[0]):
            return None
        for i in range(self.rows):
            tmp = []
            for j in range(self.columns):
                tmp.append(self.data[i][j] + other.data[i][j])
            new_data.append(tmp)
        new_matrix = Matrix(data=new_data)
        return new_matrix

    def __sub__(self, other):
        if not isinstance(other, Matrix):
            raise AttributeError("Unsupported operand type for +: '{}' and '{}'".format(self.__class__, type(other)))
        new_data = []
        if len(self.data) != len(other.data) or len(self.data[0]) != len(other.data[0]):
            return None

        for i in range(0, self.rows):
            tmp = []
            for j in range(0, self.columns):
                tmp.append(self.data[i][j] - other.data[i][j])
            new_data.append(tmp)
        new_matrix = Matrix(data=new_data)
        return new_matrix

    # Add Strassen algo for multiplication later
    def __mul__(self, other):
        if not isinstance(other, Matrix) and not isinstance(other, Number) and not isinstance(other, list):
            raise AttributeError("Unsupported operand type for +: '{}' and'{}'".format(self.__class__, type(other)))
        new_matrix = []

        if isinstance(other, Number):
            for row in self.data:
                tmp = []
                for column in row:
                    tmp.append(column * other)
                new_matrix.append(tmp)
            return Matrix(data=new_matrix)

        if type(other) == Matrix or type(other) == list:
            if len(other) != self.columns:
                raise AttributeError("Unsupported matrix dimensions for operation: '{}'x'{}'".format(len(other),
                                                                                                     len(other[0])))
            else:
                new_matrix = []
                for i in range(0, self.rows):
                    tmp = []
                    for j in range(len(other[0])):
                        tmp_sum = 0
                        for k in range(self.columns):
                            tmp_sum += self.data[i][k] * other[k][j]
                        tmp.append(tmp_sum)
                    new_matrix.append(tmp)
                return Matrix(data=new_matrix)

    def __rmul__(self, other):
        if not isinstance(other, Matrix) and not isinstance(other, Number) and not isinstance(other, list):
            raise AttributeError("Unsupported operand type for +: '{}' and '{}'".format(self.__class__, type(other)))
        new_matrix = []

        if isinstance(other, Number):
            for row in self.data:
                tmp = []
                for column in row:
                    tmp.append(column * other)
                new_matrix.append(tmp)
            return Matrix(data=new_matrix)

        if type(other) == Matrix or type(other) == list:
            if len(other[0]) != self.rows:
                raise AttributeError("Unsupported matrix dimensions for operation: '{}'x'{}'".format(len(other),
                                                                                                     len(other[0])))
            else:
                new_matrix = []
                for i in range(len(other)):
                    tmp = []
                    for j in range(self.columns):
                        tmp_sum = 0
                        for k in range(len(other[0])):
                            tmp_sum += other[i][k] * self[k][j]
                        tmp.append(tmp_sum)
                    new_matrix.append(tmp)
                return Matrix(data=new_matrix)

    def __invert__(self):
        # LUP decomposition
        self.lup()
        inversed_data = []
        for i_tmp in range(self.rows):
            inversed_data.append([])

        for i in range(self.rows):
            b = []
            for j in range(self.rows):
                if i == j:
                    b.append([1])
                else:
                    b.append([0])
            b_mat = Matrix(b)
            self.forward_substitution(b_mat)
            self.back_substitution(b_mat)

            for k in range(b_mat.rows):
                for v in range(b_mat.columns):
                    inversed_data[k].append(b_mat[k][v])
        inversed_mat = Matrix(data=inversed_data)
        return inversed_mat

    def init_p(self):
        # Init p
        for i in range(self.rows):
            tmp = []
            for j in range(self.columns):
                if i == j:
                    tmp.append(1)
                else:
                    tmp.append(0)
            if i < len(self.p):
                self.p[i] = tmp
            else:
                self.p.append(tmp)

    def transpose(self):
        if self.rows >= self.columns:

            for i in range(self.rows):
                for j in range(self.columns):
                    if i == j:
                        break
                    tmp = self.data[i][j]
                    if i >= self.columns:
                        self.data[j].append(tmp)
                    else:
                        self.data[i][j] = self.data[j][i]
                        self.data[j][i] = tmp
            for i in range(self.columns, self.rows):
                self.data.remove(self.data[i])
        else:

            for i in range(self.rows, self.columns):
                self.data.append([])
            for i in range(self.rows):
                new_row = []
                for j in range(i + 1, self.columns):
                    if j < self.rows:
                        tmp = self.data[i][j]
                        self.data[i][j] = self.data[j][i]
                        self.data[j][i] = tmp
                    else:
                        self.data[j].append(self.data[i][j])
                        self.data[i].remove(self.data[i][j])
        tmp = self.rows
        self.rows = self.columns
        self.columns = tmp

    def getelem(self, x, y):
        try:
            return self.data[x][y]
        except IndexError:
            return None

    def switch_row(self, first, second):
        if first >= self.rows or second >= self.rows:
            return None
        try:
            tmp = copy.deepcopy(self.data[first])
            self.data[first] = self.data[second]
            self.data[second] = tmp
        except IndexError:
            return None
        try:
            tmp = copy.deepcopy(self.p[first])
            self.p[first] = self.p[second]
            self.p[second] = tmp
        except IndexError:
            pass

    def switch_column(self, first, second):
        if first >= self.columns or second >= self.columns:
            return None
        for i in range(0, self.rows):
            try:
                tmp = self.data[i][first]
                self.data[i][first] = self.data[i][second]
                self.data[i][second] = tmp
            except IndexError:
                return None

    def lu(self):
        if self.rows != self.columns:
            raise AttributeError("Only square matrices are supported for LU decomposition")

        for i in range(self.rows - 1):
            if abs(self.data[i][i]) <= self.EPS:
                raise ArithmeticError("Pivot element equals zero.")
            for j in range(i + 1, self.rows):
                self.data[j][i] /= self.data[i][i]
                for k in range(i + 1, self.rows):
                    self.data[j][k] -= self.data[j][i] * self.data[i][k]

    def lup(self):
        if self.rows != self.columns:
            raise AttributeError("Only square matrices are supported for LUP decomposition.")
        self.init_p()

        for i in range(self.rows - 1):
            max_idx = i
            for k in range(i + 1, self.columns):
                # find max pivot in column
                if abs(self.data[k][i]) > abs(self.data[max_idx][i]):
                    max_idx = k
            if abs(self.data[max_idx][i]) <= self.EPS:
                raise ArithmeticError("Pivot element equals zero.")

            if max_idx != i:
                self.switch_row(i, max_idx)

            for j in range(i + 1, self.rows):
                self.data[j][i] /= self.data[i][i]
                for k in range(i + 1, self.rows):
                    self.data[j][k] -= self.data[j][i] * self.data[i][k]

    def get_u(self):
        u_mat = copy.deepcopy(self)
        for i in range(self.rows):
            for j in range(self.columns):
                if j < i:
                    u_mat[i][j] = 0
        return u_mat

    def get_l(self):
        l_mat = copy.deepcopy(self)
        for i in range(self.rows):
            for j in range(self.columns):
                if j > i:
                    l_mat[i][j] = 0
                elif j == i:
                    l_mat[i][j] = 1
        return l_mat

    def forward_substitution(self, b):
        if not isinstance(b, Matrix):
            raise AttributeError("Unsupported parameter type '{}' for parameter b, b has to be Matrix".format(type(b)))
        if self.rows != self.columns:
            raise AttributeError("Only square matrices are supported for LU decomposition")
        b.data = (self.p * b).data
        for i in range(self.rows - 1):  # maybe columns???
            for j in range(i + 1, self.rows):
                b.data[j][0] -= self.data[j][i] * b[i][0]

    def back_substitution(self, b):
        if not isinstance(b, Matrix):
            raise AttributeError("Unsupported parameter type '{}' for parameter b, b has to be Matrix".format(type(b)))
        if self.rows != self.columns:
            raise AttributeError("Only square matrices are supported for LU decomosition")
        for i in range(self.rows - 1, -1, -1):
            if abs(self.data[i][i]) < self.EPS:
                raise AttributeError("Matrix is singular.")
            b.data[i][0] /= self.data[i][i]
            for j in range(i):
                b.data[j][0] -= self.data[j][i] * b.data[i][0]

    def det(self):
        lup = self.lup()
        d = 1
        for i in range(self.rows):
            d = d * lup.lu[i][i]
        d = d * lup.sgn
        return d

    def outputfile(self, filename):
        with open(filename, "w") as f:
            print(str(self), file=f)

    @classmethod
    def fromtextfile(cls, file):
        data = []
        with open(file) as f:
            i = 0
            for line in f:
                tmp = []
                j = 0
                line_split = line.split()
                for line_char in line_split:
                    try:
                        value_tmp = float(line_char)
                        tmp.append(value_tmp)
                    except ValueError:
                        tmp.append(0)
                    j += 1
                data.append(tmp)
                i += 1
        matrix = Matrix(data=data)
        return matrix

    @classmethod
    def get_identity_matrix(cls, n):
        data = []
        for i in range(n):
            tmp = []
            for j in range(n):
                if j == i:
                    tmp.append(1)
                else:
                    tmp.append(0)
            data.append(tmp)
        matrix = Matrix(data=data)
        return matrix
