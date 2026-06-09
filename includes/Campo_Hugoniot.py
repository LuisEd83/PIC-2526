"""
Modulo: Definicao de campo de Hugoniot 


Objetivo:
-> Definir as equacoes do campo de Hugoniot
"""

import includes.Functions as fun
import includes.Auxiliar_Functions as af

from numpy import isclose, sqrt

def h0(
        z,              #Concentracao variavel entre [0,1]
        z0             #Concentracao inicial
      )-> float:
    
    if(isclose(z, z0)):
        return af.az(z) #Retorna (d/dt)(a(z)) se z = z0
    
    return (af.a(z) - af.a(z0))/(z - z0)

def sig(
        alpha,          #Variavel de controle
        z,              #Variavel pertencente ao intervalo [0,1]
        u0,             #u inicial
        v0,             #v inicial
        z0              #z inicial
        )-> float:      #retorno esperado
    """
    [Explicacao] 
    """
    num = fun.fw(u0, v0, z0)        #Numerador
    den = u0 + alpha * h0(z, z0)    #Denominador

    return num/den

##############################
def F(
    u0, #Variavel do ponto fixo U0
    v0, #Variavel do ponto fixo U0
    z0, #Variavel do ponto fixo U0
    u,  #Variavel do ponto U1
    v,  #Variavel do ponto U1
    z   #Variavel do ponto U1
):
    #Define-se os termos (para simplicidade):
    firstT  = (fun.fw(u, v, z) - fun.fw(u0, v0, z0)) * (v - v0)
    secondT = (fun.fo(u, v, z) - fun.fo(u0, v0, z0)) * (u - u0)

    return firstT - secondT

def G(
    alpha,  #Varivavel de controle
    u0,     #Varivavel do ponto fixo U0
    v0,     #Varivavel do ponto fixo U0
    z0,     #Varivavel do ponto fixo U0
    u,      #Varivavel do ponto U1
    v,      #Varivavel do ponto U1
    z       #Varivavel do ponto U1
):
    #Define-se os termos (para simplicidade):
    firstT  = (fun.fw(u, v, z) - fun.fw(u0, v0, z0)) * (u0 * (z - z0) + alpha * (af.a(z) - af.a(z0)))
    secondT = fun.fw(u0, v0, z0) * (z - z0) * (u - u0)

    return firstT - secondT

##############################
def Fu(
       u0, #u inicial
       v0, #v inicial
       z0, #z inicial
       u, #Variaveis
       v, #Variaveis
       z  #Variaveis
       )-> float:
    
    """
    [Explicacao] - Esta eh a derivada pacial de F em relacao a variavel u
    """

    #Definindo o primeiro termo:
    firstT  = fun.fwu(u, v, z) * (v - v0)
    #Definindo o segundo termo:
    secondT = fun.fou(u, v, z) * (u - u0)
    #Definindo o terceiro termo:
    thirdT  = fun.fo(u, v, z) - fun.fo(u0, v0, z0)

    return firstT - secondT - thirdT

def Fv(
       u0, #u inicial
       v0, #v inicial
       z0, #z inicial
       u, #Variaveis
       v, #Variaveis
       z  #Variaveis
       )-> float:
    
    """
    [Explicacao] - Esta eh a derivada pacial de F em relacao a variavel v
    """

    #Definindo o primeiro termo:
    firstT  = fun.fwv(u, v, z) * (v - v0)
    #Definindo o segundo termo:
    secondT = fun.fov(u, v, z) * (u - u0)
    #Definindo o terceiro termo:
    thirdT  = fun.fw(u, v, z) - fun.fw(u0, v0, z0)

    return firstT - secondT + thirdT

def Fz(
       u0, #u inicial
       v0, #v inicial
       u, #Variaveis
       v, #Variaveis
       z  #Variaveis
       )-> float:
    
    """
    [Explicacao] - Esta eh a derivada pacial de F em relacao a variavel z
    """

    #Definindo o primeiro termo:
    firstT  = fun.fwc(u, v, z) * (v - v0)
    #Definindo o segundo termo:
    secondT = fun.foc(u, v, z) * (u - u0)

    return firstT - secondT

def gradF(u0, v0, z0, u, v, z) -> list:
    return [Fu(u0, v0, z0, u, v, z), Fv(u0, v0, z0, u, v, z), Fz(u0, v0, u, v, z)]

##############################
def Gu(
       alpha,   #Variavel de controle
       u0,      #u inicial
       v0,      #v inicial
       z0,      #z inicial
       u,       #Variaveis
       v,       #Variaveis
       z        #Variaveis
       )-> float:
    
    """
    [Explicacao] - Esta eh a derivada pacial de G em relacao a variavel u
    """

    #Definindo o primeiro termo:
    firstT  = fun.fwu(u, v, z) * (u0 * (z - z0) + alpha * (af.a(z) - af.a(z0)))
    
    #Definindo o segundo termo:
    secondT = fun.fw(u0, v0, z0) * (z - z0)

    return firstT - secondT

def Gv(
       alpha,   #Variavel de controle
       u0,      #u inicial
       z0,      #z inicial
       u,       #Variaveis
       v,       #Variaveis
       z        #Variaveis
       )-> float:
    

    """
    [Explicacao] - Esta eh a derivada pacial de G em relacao a variavel v
    """

    #Definindo o primeiro termo:
    firstT  = fun.fwv(u, v, z) * (u0*(z - z0) + alpha * (af.a(z) - af.a(z0)))

    return firstT 

def Gz(
       alpha,   #Variavel de controle
       u0,      #u inicial
       v0,      #v inicial
       z0,      #z inicial
       u,       #Variaveis
       v,       #Variaveis
       z        #Variaveis
       )-> float:
    
    """
    [Explicacao] - Esta eh a derivada pacial de G em relacao a variavel z
    """

    #Definindo o primeiro termo:
    firstT  = fun.fwc(u, v, z) * (u0*(z - z0) + alpha * (af.a(z) - af.a(z0)))

    #Definindo o segundo termo:
    secondT = (fun.fw(u, v, z) - fun.fw(u0, v0, z0)) * (u0 + alpha * af.az(z))

    #Defininfo o terceiro termo
    thirdT = fun.fw(u0, v0, z0) * (u - u0)

    return firstT + secondT - thirdT

def gradG(alpha, u0, v0, z0, u, v, z):
    return [Gu(alpha, u0, v0, z0, u, v, z), Gv(alpha, u0, z0, u, v, z), Gz(alpha, u0, v0, z0, u, v, z)]

##############################

"""
Determinantes montados para solucionar o sistema (atualmente (39))
"""
def Det_u(alpha, u0, v0, z0, u, v, z):
    #definindo a diagonal principal
    diag_p = Fv(u0, v0, z0, u, v, z) * Gz(alpha, u0, v0, z0, u, v, z)

    #Definindo a diagonal secundaria
    diag_s = Fz(u0, v0, u, v, z) * Gv(alpha, u0, z0, u, v, z)

    return diag_p - diag_s

def Det_v(alpha, u0, v0, z0, u, v, z):
    #definindo a diagonal principal
    diag_p = Fu(u0, v0, z0, u, v, z) * Gz(alpha, u0, v0, z0, u, v, z)

    #Definindo a diagonal secundaria
    diag_s = Fz(u0, v0, u, v, z) * Gu(alpha, u0, v0, z0, u, v, z)

    return diag_p - diag_s

def Det_z(alpha, u0, v0, z0, u, v, z):
    #Definindo a diagonal principal
    diag_p = Fu(u0, v0, z0, u, v, z) * Gv(alpha, u0, z0, u, v, z)
    
    #Definindo a diagonal secundaria
    diag_s = Fv(u0, v0, z0, u, v, z) * Gu(alpha, u0, v0, z0, u, v, z)

    return diag_p - diag_s

#####################
def Norma(alpha, u0, v0, z0, u, v, z):
    return sqrt((Det_u(alpha, u0, v0, z0, u, v, z)**2) + 
                (Det_v(alpha, u0, v0, z0, u, v, z)**2) + 
                (Det_z(alpha, u0, v0, z0, u, v, z)**2))

def Hug1(alpha, u0, v0, z0, u, v, z):
    return Det_u(alpha, u0, v0, z0, u, v, z)/Norma(alpha, u0, v0, z0, u, v, z)

def Hug2(alpha, u0, v0, z0, u, v, z):
    return -Det_v(alpha, u0, v0, z0, u, v, z)/Norma(alpha, u0, v0, z0, u, v, z)

def Hug3(alpha, u0, v0, z0, u, v, z):
    return Det_z(alpha, u0, v0, z0, u, v, z)/Norma(alpha, u0, v0, z0, u, v, z)


def CampoHug(alpha, u0, v0, z0, u, v, z):
    return [Hug1(alpha, u0, v0, z0, u, v, z),
            Hug2(alpha, u0, v0, z0, u, v, z),
            Hug3(alpha, u0, v0, z0, u, v, z)]