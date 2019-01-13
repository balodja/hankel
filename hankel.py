#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sympy as sp


def hankel(n, m, symbols=None):
    """
    Construct the n by m Hankel matrix with variables from symbols array.
    
    If symbols is not present or None then use x_{n + m}.
    """
    if symbols is None:
        symbols = ['x_{}'.format(i) for i in range(n + m)]
    return sp.Matrix(n, m, lambda i, j: symbols[i + j])


def matrix_minor(matrix, rows, columns):
    """
    Calculate determinant of specified submatrix leaving prescribed rows and columns.
    """
    return matrix[rows, columns].det()


def check_hypothesis(matrix, zero_minor, minor1, minor2):
    """
    Check the hypothesis that minor1 equals minor2 under the condition of zero_minor = 0.
    
    Actually the function tests whether (minor1 - minor2 == ±zero_minor).
    """
    det_zero = matrix_minor(matrix, *zero_minor)
    det1 = matrix_minor(matrix, *minor1)
    det2 = matrix_minor(matrix, *minor2)
    return (det_zero - det1 + det2).equals(0) or (det_zero + det1 - det2).equals(0)


class HypothesisLogger():
    """
    Log all the tested hypotheses, and format them to specified tex-file later.
    """

    def __init__(self, n, m):
        self.matrix_name = 'M'
        self.matrix = hankel(n, m)
        self.short_results = []
        self.long_results = []
        
        matrix_itself = '$$ ' + self.matrix_name + ' = ' + sp.latex(self.matrix) + ' $$'
        self.short_results.append(matrix_itself)
        self.long_results.append(matrix_itself)
    
    def format_minor(self, minor, short=True):
        """
        Make latex for the minor.
        
        Either short one with just rows-columns specification,
        or long one with the values themselves (when short=False).
        """
        
        rows, columns = minor
        if short:
            mat = sp.Matrix([rows, columns])
            return self.matrix_name + sp.latex(mat)
        else:
            mat = self.matrix[rows, columns]
            return r'\det ' + sp.latex(mat)
    
    def log_hypothesis(self, zero_minor, minor1, minor2):
        """
        Test the hypothesis about minors and log the results.
        """
        
        result = check_hypothesis(self.matrix, zero_minor, minor1, minor2)
        tex_iff = '\Leftrightarrow' if result else '\nLeftrightarrow'
        
        def format_results(short):
            tex_zero = self.format_minor(zero_minor, short)
            tex1 = self.format_minor(minor1, short)
            tex2 = self.format_minor(minor2, short)
            template = '$$ {} = 0 {} {} = {} $$'
            return template.format(tex_zero, tex_iff, tex1, tex2)
        
        self.short_results.append(format_results(True))
        self.long_results.append(format_results(False))

    def save(self, filename, short=True):
        """
        Save the logged results to the file.
        """

        with open(filename, 'w') as file:
            results = self.short_results if short else self.long_results
            file.writelines(line + '\n' for line in results)


if __name__ == '__main__':
    logger = HypothesisLogger(7, 7)
    logger.log_hypothesis(([0, 1, 2], [0, 2, 3]), ([0, 1, 4], [0, 1, 2]), ([0, 1, 3], [0, 1, 3]))
    logger.log_hypothesis(([0, 1, 2], [0, 2, 4]), ([0, 1, 5], [0, 1, 2]), ([0, 1, 4], [0, 1, 3]))
    logger.log_hypothesis(([0, 1, 2], [0, 2, 5]), ([0, 1, 6], [0, 1, 2]), ([0, 1, 5], [0, 1, 3]))
    logger.log_hypothesis(([0, 1, 2], [0, 3, 4]), ([0, 1, 5], [0, 1, 3]), ([0, 1, 4], [0, 1, 4]))
    logger.save('payload_short.tex', True)
    logger.save('payload_long.tex', False)
