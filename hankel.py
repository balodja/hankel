#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sympy as sp


def hankel_submatrix(symbol, rows, columns):
    """
    Calculate the Hankel submatrix at the specified rows and columns.
    """

    return sp.Matrix([[symbol[r + c] for c in columns] for r in rows])


def check_conjecture(symbol, minor1, minor2, minor12):
    """
    Check the conjecture that minor1 plus minor2 equals minor12 for Hankel matrix.
    
    Actually the function tests whether (minor1 + minor2 - minor12 == 0).
    """
    det12 = hankel_submatrix(symbol, *minor12).det()
    det1 = hankel_submatrix(symbol, *minor1).det()
    det2 = hankel_submatrix(symbol, *minor2).det()
    return (det1 + det2 - det12).equals(0)


class ConjectureLogger():
    """
    Log all the tested hypotheses, and format them to specified tex-file later.
    """

    def __init__(self):
        self.matrix_symbol = sp.IndexedBase('H')
        self.element_symbol = sp.IndexedBase('h')
        self.short_results = []
        self.long_results = []
    
    def format_minor(self, minor, short=True):
        """
        Make latex for the minor.
        
        Either short one with just rows-columns specification,
        or long one with the values themselves (when short=False).
        """
        
        rows, columns = minor
        if short:
            mat = sp.Matrix([rows, columns])
            return sp.latex(self.matrix_symbol) + sp.latex(mat)
        else:
            mat = hankel_submatrix(self.element_symbol, rows, columns)
            return r'\det ' + sp.latex(mat)
    
    def log_conjecture(self, minor1, minor2, minor12):
        """
        Test the conjecture about minors and log the results.
        """
        
        result = check_conjecture(self.element_symbol, minor1, minor2, minor12)
        tex_eq = r'=' if result else r'\neq'
        
        def format_results(short):
            tex1 = self.format_minor(minor1, short)
            tex2 = self.format_minor(minor2, short)
            tex12 = self.format_minor(minor12, short)
            template = '$$ {} + {} {} {} $$'
            return template.format(tex1, tex2, tex_eq, tex12)
        
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
    # Just one submatrix for clearance
    sample_matrix = hankel_submatrix(sp.IndexedBase('h'), range(5), range(5))
    with open('sample_hankel.tex', 'w') as file:
        file.write('H = ')
        file.write(sp.latex(sample_matrix))

    # Now check some conjectures
    logger = ConjectureLogger()
    i, j, k, l = sp.symbols('i j k l', cls=sp.Idx)
    logger.log_conjecture(([0, 1, 2], [0, j, i - 1]), ([0, 1, i], [0, 1, j]), ([0, 1, i - 1], [0, 1, j + 1]))
    logger.log_conjecture(([0, 1, 3], [0, j, i - 1]), ([0, 1, i], [0, 2, j]), ([0, 2, i - 1], [0, 1, j + 1]))
    logger.log_conjecture(([0, 1, l + 1], [0, j, i - 1]), ([0, 1, i], [0, l, j]), ([0, l, i - 1], [0, 1, j + 1]))
    logger.log_conjecture(([k, 1, l + 1], [0, j, i - 1]), ([k, 1, i], [0, l, j]), ([0, l, i - 1], [k, 1, j + 1]))
    logger.log_conjecture(([0, 1, 2, 3], [0, 1, j, i - 1]),
                          ([0, 1, 2, i], [0, 1, 2, j]), ([0, 1, 2, i - 1], [0, 1, 2, j + 1]))
    logger.log_conjecture(([0, 1, 2, 3, 4], [0, 1, 2, j, i - 1]),
                          ([0, 1, 2, 3, i], [0, 1, 2, 3, j]), ([0, 1, 2, 3, i - 1], [0, 1, 2, 3, j + 1]))
    logger.save('payload_short.tex', True)
    logger.save('payload_long.tex', False)
