#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sympy as sp
import abc


class HankelMatrix:
    """
    Simple class for holding methods like evaluating minors or formatting submatrices.
    """
    def __init__(self, matrix_symbol=None, element_symbol=None):
        if matrix_symbol is None:
            matrix_symbol = sp.IndexedBase('H')

        self.matrix_symbol = matrix_symbol

        if element_symbol is None:
            element_symbol = matrix_symbol.copy()
            element_symbol.label.name = element_symbol.label.name.lower()

        self.element_symbol = element_symbol

    def minor(self, rows, columns):
        """
        Calculate the minor at the specified rows and columns
        :param rows: array of row indices
        :param columns: array of column indices
        :return: determinant expression
        """
        return self.submatrix(rows, columns).det()

    def submatrix(self, rows, columns):
        """
        Calculate the Hankel submatrix at the specified rows and columns.
        """

        return sp.Matrix([[self.element_symbol[r + c] for c in columns] for r in rows])

    def format_minor(self, rows, columns, short=True):
        """
        Make latex for the minor.

        Either short one with just rows-columns specification,
        or long one with the values themselves (when short=False).
        """

        if short:
            mat = sp.Matrix([rows, columns])
            return sp.latex(self.matrix_symbol) + sp.latex(mat)
        else:
            mat = self.submatrix(rows, columns)
            return r'\det ' + sp.latex(mat)


class Conjecture(abc.ABC):
    """
    Abstract class with some methods for testing and formatting hypotheses.
    """
    def __init__(self):
        self.valid = None

    @abc.abstractmethod
    def check(self):
        """
        Verify the conjecture.
        :return: Boolean result of the verification.
        """
        return NotImplemented

    def format(self, short):
        """
        Format LaTeX code stating the (in)validity of conjecture.
        :param short: use short or long form of minors
        :return: tex code in str
        """
        if self.valid is None:
            self.valid = self.check()

        return self._format(short)

    @abc.abstractmethod
    def _format(self, short):
        """
        Format LaTeX code stating the (in)validity of conjecture.
        The real work is done here. The `format` method is a wrapper.
        :param short: use short or long form of minors
        :return: tex code in str
        """
        return NotImplemented


class Conjecture_1plus2equals3(Conjecture):
    """
    The hypothesis states that (minor1 + minor2 == minor3).
    """
    def __init__(self, matrix, minor1, minor2, minor3):
        super().__init__()
        self.matrix = matrix
        self.minor1 = minor1
        self.minor2 = minor2
        self.minor3 = minor3

    def check(self):
        det1 = self.matrix.minor(*self.minor1)
        det2 = self.matrix.minor(*self.minor2)
        det3 = self.matrix.minor(*self.minor3)
        self.valid = (det1 + det2).equals(det3)
        return self.valid

    def _format(self, short=True):

        tex_eq = r'=' if self.valid else r'\neq'

        tex1 = self.matrix.format_minor(*self.minor1, short=short)
        tex2 = self.matrix.format_minor(*self.minor2, short=short)
        tex3 = self.matrix.format_minor(*self.minor3, short=short)
        return '{} + {} {} {}'.format(tex1, tex2, tex_eq, tex3)


class Conjecture_increment(Conjecture):
    """
    The hypothesis states that sum of minors with one-at-a-time augmented indices by some value
    is independent whether rows are augmented either columns.
    """
    def __init__(self, matrix, minor, delta):
        super().__init__()
        self.matrix = matrix
        self.base_minor = minor
        self.delta = delta

        # Obtain the minors for summation
        self.row_minors = []
        self.column_minors = []
        rows, columns = self.base_minor
        for i in range(len(rows)):
            incremented_rows = rows[:]
            incremented_rows[i] = rows[i] + self.delta
            self.row_minors.append((incremented_rows, columns))
            incremented_columns = columns[:]
            incremented_columns[i] = columns[i] + self.delta
            self.column_minors.append((rows, incremented_columns))

    def check(self):
        expr = 0
        for minor in self.row_minors:
            expr += self.matrix.minor(*minor)

        for minor in self.column_minors:
            expr -= self.matrix.minor(*minor)

        return expr.equals(0)

    def _format(self, short=True):

        tex_eq = r'=' if self.valid else r'\neq'

        row_sum = ' + '.join([sp.latex(self.matrix.format_minor(*minor, short=short)) for minor in self.row_minors])
        column_sum = ' + '.join([sp.latex(self.matrix.format_minor(*minor, short=short)) for minor in self.column_minors])
        return '{} {} {}'.format(row_sum, tex_eq, column_sum)


class ConjectureLogger:
    """
    Class for logging the conjectures.
    """
    def __init__(self, matrix):
        """
        :param matrix: The HankelMatrix to state conjectures about.
        """
        self.matrix = matrix
        self.short_results = []
        self.long_results = []

    def log(self, conjecture_class, *args):
        """
        Test the conjecture about minors and log the results.
        """

        conjecture = conjecture_class(self.matrix, *args)
        conjecture.check()
        
        self.short_results.append(conjecture.format(True))
        self.long_results.append(conjecture.format(False))

    def save(self, filename, short=True):
        """
        Save the logged results to the file.
        """

        with open(filename, 'w') as file:
            results = self.short_results if short else self.long_results
            file.writelines('\\begin{dmath}' + line + '\\end{dmath}\n' for line in results)


def main():
    matrix = HankelMatrix()
    with open('sample_hankel.tex', 'w') as file:
        file.write('H = ')
        file.write(sp.latex(matrix.submatrix(range(5), range(5))))

    # Now check some conjectures
    logger = ConjectureLogger(matrix)
    i, j, k, l, m, n, p, q, t = sp.symbols('i j k l m n p q t', cls=sp.Idx)
    logger.log(Conjecture_1plus2equals3, ([0, 1, 2], [0, j, i - 1]), ([0, 1, i], [0, 1, j]),
               ([0, 1, i - 1], [0, 1, j + 1]))
    logger.log(Conjecture_1plus2equals3, ([0, 1, 3], [0, j, i - 1]), ([0, 1, i], [0, 2, j]),
               ([0, 2, i - 1], [0, 1, j + 1]))
    logger.log(Conjecture_1plus2equals3, ([0, 1, l + 1], [0, j, i - 1]), ([0, 1, i], [0, l, j]),
               ([0, l, i - 1], [0, 1, j + 1]))
    logger.log(Conjecture_1plus2equals3, ([k, 1, l + 1], [0, j, i - 1]), ([k, 1, i], [0, l, j]),
               ([0, l, i - 1], [k, 1, j + 1]))
    logger.log(Conjecture_1plus2equals3, ([0, 1, 2, 3], [0, 1, j, i - 1]),
               ([0, 1, 2, i], [0, 1, 2, j]), ([0, 1, 2, i - 1], [0, 1, 2, j + 1]))
    logger.log(Conjecture_1plus2equals3, ([0, 1, 2, 3, 4], [0, 1, 2, j, i - 1]),
               ([0, 1, 2, 3, i], [0, 1, 2, 3, j]), ([0, 1, 2, 3, i - 1], [0, 1, 2, 3, j + 1]))
    logger.log(Conjecture_increment, ([i, j], [l, m]), t)
    logger.log(Conjecture_increment, ([i, j, k], [l, m, n]), t)
    logger.log(Conjecture_increment, ([i, j, k, l], [m, n, p, q]), t)
    logger.save('payload_short.tex', True)
    logger.save('payload_long.tex', False)


if __name__ == '__main__':
    main()
