import math
from ..utilities import *
from ..utilities.exceptions import *
from ..utilities.maths import *

from .cell import UnitCell

from math import radians, cos, copysign

# Implementation from:
# https://journals.iucr.org/a/issues/1976/02/00/a12875/a12875.pdf


def sign(x):
    if x >= 0:
        return 1
    elif x < 0:
        return -1
    return ZeroDivisionError()

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

    ξ = 2*b*c*cos(alpha) # ξ (Xi)
    η = 2*a*c*cos(beta) # η (Eta)
    ζ = 2*a*b*cos(gamma) # ζ (Zeta)

    log("header", f"A={A:.3f} B={B:.3f} C={C:.3f}, ξ={ξ:.3f}, η={η:.3f}, ζ={ζ:.3f}")
    print()

    # Loop #############################################################################################################
    done = False
    n = 0
    while not (done or n >= 6):
        n+=1
        # 1 ############################################################################################################

        if (A > B) or ( (round(A)==round(B)) and (abs(ξ) > abs(η)) ):
            log(1, "1", UnitCell(A, B, C, ξ, η, ζ))
            A, ξ, B, η = B, η, A, ξ#


        # 2 ############################################################################################################

        if (B > C) or ( (round(B)==round(C)) and (abs(η) > abs(ζ)) ):
            log(1, "2", UnitCell(A, B, C, ξ, η, ζ))
            B, η, C, ζ =  C, ζ, B, η
            continue

        # 3 ############################################################################################################

        if ξ*η*ζ > 0:
        #if ξ*η*ξ > 0: # ??? typo?
            log(1, "3", UnitCell(A, B, C, ξ, η, ζ))
            ξ, η, ζ = abs(ξ), abs(η), abs(ζ)

        # 4 ############################################################################################################

        if ξ*η*ζ <= 0:
            log(1, "4", UnitCell(A, B, C, ξ, η, ζ))
            ξ, η, ζ = -abs(ξ), -abs(η), -abs(ζ)

        # 5 ############################################################################################################

        if (abs(ξ) > B) or ( (round(ξ) == round(B)) and (2*η < ζ) ) or ( (round(ξ) == round(-B)) and (ζ < 0) ):
            log(1, "5", UnitCell(A, B, C, ξ, η, ζ))
            # C = copysign(B+C-ξ, ξ)
            # η = copysign(η-ζ, ξ)
            # ξ = copysign(ξ-(2*B), ξ)
            C = B+C-(ξ*sign(ξ))
            η = η-(ζ*sign(ξ))
            ξ = ξ-((2*B)*sign(ξ))
            continue

        # 6 ############################################################################################################

        if (abs(η) > B) or ( (round(η) == round(B)) and (2*ξ < ζ) ) or ( (round(η) == round(-B)) and (ζ < 0) ):
            log(1, "6", UnitCell(A, B, C, ξ, η, ζ))
            C = A+C-(η*sign(η))
            ξ = ξ-(ζ*sign(η))
            η = η-(2*A*sign(η))
            continue

        # 7 ############################################################################################################

        if (abs(ζ) > A) or ( (round(ζ) == round(A)) and (2*ξ < η) ) or ( (round(ζ) == round(-A)) and (η < 0) ):
            log(1, "7", UnitCell(A, B, C, ξ, η, ζ))
            B = A+B-(ζ*sign(ζ))
            ξ = ξ-(η*sign(ζ))
            ζ = ζ-(2*A*sign(ζ))
            continue

        # 8 ############################################################################################################

        if (ξ+η+ζ+A+B < 0) or ( (round(ξ+η+ζ+A+B) == 0) and ((2*A)+(2*η)+ζ > 0) ):
            log(1, "8", UnitCell(A, B, C, ξ, η, ζ))
            C = A+B+C+ξ+η+ζ
            ξ = (2*B)+ξ+ζ
            η = (2*A)+η+ζ
            continue

        done = True
    log(1, ">", UnitCell(A, B, C, ξ, η, ζ))
    ####################################################################################################################
    niggli_cell = UnitCell(A, B, C, ξ, η, ζ)

    return niggli_cell

