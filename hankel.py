#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sympy as sp


def hankel_submatrix(symbol, rows, columns):
    """
    Calculate the Hankel submatrix at the specified rows and columns.
    """

    return sp.Matrix([[symbol[r + c] for c in columns] for r in rows])


def check_conjecture(symbol, zero_minor, minor1, minor2):
    """
    Check the conjecture that minor1 equals minor2 under the condition of zero_minor = 0. For Hankel matrix.
    
    Actually the function tests whether (minor1 - minor2 == ±zero_minor).
    """
    det_zero = hankel_submatrix(symbol, *zero_minor).det()
    det1 = hankel_submatrix(symbol, *minor1).det()
    det2 = hankel_submatrix(symbol, *minor2).det()
    return (det_zero - det1 + det2).equals(0) or (det_zero + det1 - det2).equals(0)


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
    
    def log_conjecture(self, zero_minor, minor1, minor2):
        """
        Test the conjecture about minors and log the results.
        """
        
        result = check_conjecture(self.element_symbol, zero_minor, minor1, minor2)
        tex_iff = r'\Leftrightarrow' if result else r'\nLeftrightarrow'
        
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
    # Just one submatrix for clearance
    sample_matrix = hankel_submatrix(sp.IndexedBase('h'), range(5), range(5))
    with open('sample_hankel.tex', 'w') as file:
        file.write('H = ')
        file.write(sp.latex(sample_matrix))

    # Now check some conjectures
    logger = ConjectureLogger()
    i, j, k, l = sp.symbols('i j k l', cls=sp.Idx)
    logger.log_conjecture(([0, 1, 2], [0, j, i - 1]), ([0, 1, i], [0, 1, j]), ([0, 1, i - 1], [0, 1, j + 1]))
    logger.log_conjecture(([0, 3, 1], [0, j, i - 1]), ([0, 1, i], [0, 2, j]), ([0, 2, i - 1], [0, 1, j + 1]))
    logger.log_conjecture(([0, l + 1, 1], [0, j, i - 1]), ([0, 1, i], [0, l, j]), ([0, l, i - 1], [0, 1, j + 1]))
    logger.save('payload_short.tex', True)
    logger.save('payload_long.tex', False)
