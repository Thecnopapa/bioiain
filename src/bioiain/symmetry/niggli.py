import math
from ..utilities import *
from ..utilities.exceptions import *
from ..utilities.maths import *

from .cell import UnitCell

from math import radians, cos, copysign, sqrt, acos, degrees

# Implementation from:
# https://journals.iucr.org/a/issues/1976/02/00/a12875/a12875.pdf





def cell_to_niggli(cell, angles_are_radians=False, verbose=False, return_matrix=False, eps=0.01, max_iter=100):
    def equal(a, b):
        return abs(a-b) < eps

    def greater(a, b):
        return a-b > eps

    def lower(a, b):
        return b-a > eps

    def sign(x):
        if x >= 0:
            return 1
        elif x < 0:
            return -1
        return ZeroDivisionError()

    log(2, "Calculating niggli cell...")
    a = cell.a
    b = cell.b
    c = cell.c

    alpha = cell.alpha
    beta = cell.beta
    gamma = cell.gamma

    if not angles_are_radians:
        alpha = radians(alpha)
        beta = radians(beta)
        gamma = radians(gamma)
    if verbose:
        log("header", f"a={a:.3f} b={b:.3f} c={c:.3f}, alpha={alpha:.3f}, beta={beta:.3f}, gamma={gamma:.3f}")

    # 0 ################################################################################################################

    A = a*a
    B = b*b
    C = c*c

    ξ = 2*b*c*cos(alpha) # ξ (Xi)
    η = 2*a*c*cos(beta) # η (Eta)
    ζ = 2*a*b*cos(gamma) # ζ (Zeta)
    if verbose:
        log("header", f"A={A:.3f} B={B:.3f} C={C:.3f}, ξ={ξ:.3f}, η={η:.3f}, ζ={ζ:.3f}")
        print()

    # Loop #############################################################################################################
    done = False
    n = 0
    while not (done or n >= max_iter):
        #raise Exception("not ready yet")
        n+=1
        # 1 ############################################################################################################

        if greater(A, B) or ( equal(A,B) and greater(abs(ξ), abs(η)) ):
            if verbose:
                log(1, "1", UnitCell(A, B, C, ξ, η, ζ))
            A, ξ, B, η = B, η, A, ξ#


        # 2 ############################################################################################################

        if greater(B, C) or ( equal(B, C) and greater(abs(η), abs(ζ)) ):
            if verbose:
                log(1, "2", UnitCell(A, B, C, ξ, η, ζ))
            B, η, C, ζ =  C, ζ, B, η
            continue

        # 3 ############################################################################################################

        if greater(ξ*η*ζ, 0):
            if verbose:
                log(1, "3", UnitCell(A, B, C, ξ, η, ζ))
            ξ, η, ζ = abs(ξ), abs(η), abs(ζ)

        # 4 ############################################################################################################

        if not greater(ξ*η*ζ, 0):
            if verbose:
                log(1, "4", UnitCell(A, B, C, ξ, η, ζ))
            ξ, η, ζ = -abs(ξ), -abs(η), -abs(ζ)

        # 5 ############################################################################################################

        if greater(abs(ξ), B) or ( equal(ξ,B) and lower(2*η, ζ) ) or ( equal(ξ, -B) and lower(ζ, 0) ):
            if verbose:
                log(1, "5", UnitCell(A, B, C, ξ, η, ζ))
            C = B+C-(ξ*sign(ξ))
            η = η-(ζ*sign(ξ))
            ξ = ξ-((2*B)*sign(ξ))
            continue

        # 6 ############################################################################################################
        if greater(abs(η), A) or ( equal(η, A) and lower(2*ξ, ζ) ) or ( equal(η, -A) and lower(ζ, 0) ):
            if verbose:
                log(1, "6", UnitCell(A, B, C, ξ, η, ζ))
            C = A+C-(η*sign(η))
            ξ = ξ-(ζ*sign(η))
            η = η-(2*A*sign(η))
            continue

        # 7 ############################################################################################################

        if greater(abs(ζ), A) or ( equal(ζ, A) and lower(2*ξ, η) ) or ( equal(ζ, -A) and lower(η,  0) ):
            if verbose:
                log(1, "7", UnitCell(A, B, C, ξ, η, ζ))
            B = A+B-(ζ*sign(ζ))
            ξ = ξ-(η*sign(ζ))
            ζ = ζ-(2*A*sign(ζ))
            continue

        # 8 ############################################################################################################

        if lower(ξ+η+ζ+A+B, 0) or ( equal(ξ+η+ζ+A+B, 0) and greater((2*A)+(2*η)+ζ, 0) ):
            if verbose:
                log(1, "8", UnitCell(A, B, C, ξ, η, ζ))
            C = A+B+C+ξ+η+ζ
            ξ = (2*B)+ξ+ζ
            η = (2*A)+η+ζ
            continue

        done = True
    if verbose:
        log(1, ">", UnitCell(A, B, C, ξ, η, ζ))
        log(1, f"N-iterations: {n}")
        print()
    ####################################################################################################################

    if return_matrix:
        return np.array([[A, B, C], [ξ/2, η/2, ζ/2 ]])


    a, b, c = sqrt(A), sqrt(B), sqrt(C)
    alpha = degrees(acos(ξ/2/b/c))
    beta = degrees(acos(η/2/a/c))
    gamma = degrees(acos(ζ/2/a/b))
    niggli_cell = UnitCell(a, b, c, alpha, beta, gamma)
    return niggli_cell

