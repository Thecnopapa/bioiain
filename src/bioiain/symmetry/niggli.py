import math
from ..utilities import *
from ..utilities.exceptions import *
from ..utilities.maths import *

from .cell import UnitCell

from math import radians, cos, copysign

# Implementation from:
# https://journals.iucr.org/a/issues/1976/02/00/a12875/a12875.pdf

def cell_to_niggli(cell, angles_are_radians=False):
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

    log("header", f"a={a:.3f} b={b:.3f} c={c:.3f}, alpha={alpha:.3f}, beta={beta:.3f}, gamma={gamma:.3f}")

    # 0 ################################################################################################################

    A = a*a
    B = b*b
    C = c*c

    ξ = -2*b*c*cos(alpha) # ξ (Xi)
    η = -2*a*c*cos(beta) # η (Eta)
    ζ = -2*a*b*cos(gamma) # ζ (Zeta)

    log("header", f"A={A:.3f} B={B:.3f} C={C:.3f}, ξ={ξ:.3f}, η={η:.3f}, ζ={ζ:.3f}")
    print()

    # Loop #############################################################################################################
    done = False
    n = 0
    while not (done or n > 1000):
        n+=1
        # 1 ############################################################################################################
        log(1, "1", UnitCell(A, B, C, ξ, η, ζ))
        if (A > B) or ( (A==B) and (abs(ξ) > abs(η)) ):
            A, ξ, B, η = B, η, A, ξ

        # 2 ############################################################################################################
        log(1, "2", UnitCell(A, B, C, ξ, η, ζ))
        if (B > C) or ( (B==C) and (abs(η) > abs(ζ)) ):
            B, η, C, ζ =  C, ζ, B, η
            continue

        # 3 ############################################################################################################
        log(1, "3", UnitCell(A, B, C, ξ, η, ζ))
        #if ξ*η*ζ > 0:
        if ξ*η*ξ > 0: # ??? typo?
            ξ, η, ζ = abs(ξ), abs(η), abs(ζ)

        # 4 ############################################################################################################
        log(1, "4", UnitCell(A, B, C, ξ, η, ζ))
        if ξ*η*ζ <= 0:
            ξ, η, ζ = -abs(ξ), -abs(η), -abs(ζ)

        # 5 ############################################################################################################
        log(1, "5", UnitCell(A, B, C, ξ, η, ζ))
        if (abs(ξ) > B) or ( (ξ == B) and (2*η < ζ) ) or ( (ξ == -B) and (ζ < 0) ):
            C = copysign(B+C-ξ, ξ)
            η = copysign(η-ζ, ξ)
            ξ = copysign(ξ-(2*B), ξ)
            continue

        # 6 ############################################################################################################
        log(1, "6", UnitCell(A, B, C, ξ, η, ζ))
        if (abs(η) > B) or ( (η == B) and (2*ξ < ζ) ) or ( (η == -B) and (ζ < 0) ):
            C = copysign(A+C-η, η)
            ξ = copysign(ξ-ζ, η)
            η = copysign(η-(2*A), η)
            continue

        # 7 ############################################################################################################
        log(1, "7", UnitCell(A, B, C, ξ, η, ζ))
        if (abs(ζ) > B) or ( (ζ == B) and (2*ξ < η) ) or ( (ζ == -B) and (η < 0) ):
            B = copysign(A+B-ζ, ζ)
            ξ = copysign(ζ-η, ζ)
            ζ = copysign(ξ-(2*A), ζ)
            continue

        # 8 ############################################################################################################
        log(1, "8", UnitCell(A, B, C, ξ, η, ζ))
        if (ξ+η+ζ+A+B < 0) or ( (ξ+η+ζ+A+B == 0) and ((2*A)+(2*η)+ζ > 0) ):
            C = A+B+C+ξ+η+ζ
            ξ = (2*B)+ξ+ζ
            η = (2*A)+η+ζ
            continue
        done = True


    ####################################################################################################################
    niggli_cell = UnitCell(A, B, C, ξ, η, ζ)

    return niggli_cell

