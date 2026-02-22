"""
Modulo: Curvas integrais

Objetivo:
- Definir funcoes determinantes a partir de funcoes pre-estabelecidas
pela biblioteca Functions.
- Definir funcao norma N
- Definir as componentes do campo

"""

import Functions as fun
from numpy import sqrt

#Definicao da funcao Dp
def D_p(u, v, z):
    #Variaveis da multiplicacao da diagonal principal da matriz jacobiana 2x2
    fwv = fun.fwv(u, v, z)
    goz = fun.foc(u, v, z)

    #Variaveis da multiplicacao da diagonal secundaria da matriz jacobiana 2x2
    gov = fun.fov(u, v, z)
    fwz = fun.fwc(u, v, z)
    lbdz = fun.lbdc(u, v, z)

    #Diagonais da matriz jacobiana
    Diag_p = fwv*goz
    Diag_s = (gov - lbdz)*fwz

    return Diag_p - Diag_s

#Definicao da funcao Dq
def D_q(u, v, z):
    #Variaveis da multiplicacao da diagonal principal da matriz jacobiana 2x2
    fwz = fun.fwc(u, v, z)
    gou = fun.fou(u, v, z)

    #Variaveis da multiplicacao da diagonal secundaria da matriz jacobiana 2x2
    fwu = fun.fwu(u, v, z)
    lbdz = fun.lbdc(u, v, z)
    goz = fun.foc(u, v, z)

    #Diagonais da matriz jacobiana
    Diag_p = fwz*gou
    Diag_s = (fwu - lbdz)*goz

    return Diag_p - Diag_s

#Definicao da funcao Dr
def D_r(u, v, z):
    #Variaveis da multiplicacao da diagonal principal da matriz jacobiana 2x2
    fwu = fun.fwu(u, v, z)
    gov = fun.fov(u, v, z)
    lbdz = fun.lbdc(u, v, z)
    
    #Variaveis da multiplicacao da diagonal secundaria da matriz jacobiana 2x2
    fwv = fun.fwv(u, v, z)
    gou = fun.fou(u, v, z)

    #Diagonais da matriz jacobiana
    Diag_p = fwu*gov - (fwu + gov) * lbdz + lbdz**2
    Diag_s = fwv*gou

    return Diag_p - Diag_s

#Definicao da funcao norma N
def N(u, v, z):
    return sqrt((D_r(u, v, z))**2 + (D_p(u, v, z))**2 + (D_q(u, v, z))**2)

#Definicao dos componentes do campo  |SEPARADOS|
def P(u, v, z):
    return D_p(u, v, z)/N(u, v, z)

def Q(u, v, z):
    return D_q(u, v, z)/ N(u, v, z)

def R(u, v, z):
    return D_r(u, v, z)/N(u, v, z)

#Componentes do campo JUNTOS
def campos(u, v, z):
    return P(u, v, z), Q(u, v, z), R(u, v, z)