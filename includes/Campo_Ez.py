"""
Modulo: Funcoes determinates

Objetivo:
- Definir funcoes determinantes a partir de funcoes pre-estabelecidas
pela biblioteca Functions.
- Definir funcao norma N
- Definir as componentes do campo

Com isto concluido, poderemos passar para proxima etapa: integrar as componentes do campo por um metodo de integracao  
"""

import includes.Functions as fun

from numpy import sqrt, cos

#Redefinindo funcoes para nao haver erro circular
#def az(z):
#    return (cos(z))

def az(z):
    return 1/(1+z**2)

def lambdz(u, v, z, alpha):
    return(fun.fw(u, v, z)/(u + alpha*az(z)))

#Definicao da funcao Dp
def D_p(u, v, z, alpha):
    #Variaveis da multiplicacao da diagonal principal da matriz jacobiana 2x2
    fwv = fun.fwv(u, v, z)
    goz = fun.foc(u, v, z)

    #Variaveis da multiplicacao da diagonal secundaria da matriz jacobiana 2x2
    gov = fun.fov(u, v, z)
    fwz = fun.fwc(u, v, z)
    lbdz = lambdz(u, v, z, alpha)

    #Diagonais da matriz jacobiana
    Diag_p = fwv*goz
    Diag_s = (gov - lbdz)*fwz

    return Diag_p - Diag_s

#Definicao da funcao Dq
def D_q(u, v, z, alpha):
    #Variaveis da multiplicacao da diagonal principal da matriz jacobiana 2x2
    fwz = fun.fwc(u, v, z)
    gou = fun.fou(u, v, z)

    #Variaveis da multiplicacao da diagonal secundaria da matriz jacobiana 2x2
    fwu = fun.fwu(u, v, z)
    lbdz = lambdz(u, v, z, alpha)
    goz = fun.foc(u, v, z)

    #Diagonais da matriz jacobiana
    Diag_p = fwz*gou
    Diag_s = (fwu - lbdz)*goz

    return Diag_p - Diag_s

#Definicao da funcao Dr
def D_r(u, v, z, alpha):
    #Variaveis da multiplicacao da diagonal principal da matriz jacobiana 2x2
    fwu = fun.fwu(u, v, z)
    gov = fun.fov(u, v, z)
    lbdz = lambdz(u, v, z, alpha)
    
    #Variaveis da multiplicacao da diagonal secundaria da matriz jacobiana 2x2
    fwv = fun.fwv(u, v, z)
    gou = fun.fou(u, v, z)

    #Diagonais da matriz jacobiana
    Diag_p = fwu*gov - (fwu + gov) * lbdz + lbdz**2
    Diag_s = fwv*gou

    return Diag_p - Diag_s

#Definicao da funcao norma N
def N(u, v, z, alpha):
    return sqrt((D_r(u, v, z, alpha))**2 + (D_p(u, v, z, alpha))**2 + (D_q(u, v, z, alpha))**2)

#Componentes do campo SEPARADOS
def P(u, v, z, alpha): #Componente multiplo de i
    return D_p(u, v, z, alpha)/N(u, v, z, alpha)

def Q(u, v, z, alpha): #Componente multiplo de j
    return D_q(u, v, z, alpha)/ N(u, v, z, alpha)

def R(u, v, z, alpha): #Componente multiplo de k
    return D_r(u, v, z, alpha)/N(u, v, z, alpha)

#Componentes do campo JUNTOS
def campos(u, v, z, alpha):
    return [P(u, v, z, alpha), Q(u, v, z, alpha), R(u, v, z, alpha)]