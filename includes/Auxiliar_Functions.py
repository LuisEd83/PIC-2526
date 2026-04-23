"""
Modulo: Funcoes auxiliares

Objetivo:
- Armazenar as funcoes que serao utilizadas no codigo e selecionar as necessarias

"""

import includes.Numericals_methods as nm
import includes.Functions as fun

import numpy as np

#Definindo funcao para determinar se o usuário quer ver ou nao os pontos dentro do prisma
def transparencia():
    x = 1 #x = 1 para ver apenas os pontos dentro do prisma; x = 0 e para ver todos os pontos
    return x

#Definindo uma funcao responsavel por permitir a plotagem dos ramos relacionados a curva de
#nivel lambda_s = lambda_z
def branchSlow(): 
    s = 1 #s = 1 para permitir que seja plotado o ramo
    return s

#Definindo uma funcao responsavel por permitir a plotagem dos ramos relacionados a curva de
#nivel lambda_f = lambda_z
def branchFast(): 
    f = 1 #s = 1 para permitir que seja plotado o ramo
    return f

def points_curv(): #Esta funcao habilita os pontos na curva integral.
    x = 1   #x = 1 para habilitar
    return x

##############################################################
#Derivada da funcao continua e crescente no intervalo [0,1]
def az(z): #Da/dz
    return np.cos(z) #Derivada de a(z) em relacao a z

def lambdz(u, v, z, alpha):
    return(fun.fw(u, v, z)/(u + alpha*az(z)))

#derivada de lambdz
def Dlambdz(u, v, z, alpha):
    from includes.Campo_Ez import campos

    #Calculando gradiente:
    pc = (fun.Du(u, v, z) * (u + alpha*az(z)) - fun.fw(u, v, z))/((u+alpha*az(z))**2)                       #Primeira componente
    sc = (fun.Dv(u, v, z))/(u + alpha*az(z))                                                                #Segunda componente
    tc = (fun.Dc(u, v, z) * (u + alpha*az(z)) + fun.fw(u, v, z) * (alpha*np.sin(z)))/((u + alpha*az(z))**2) #Terceira componente

    #Extraindo campos (componentes de direção)
    camps = campos(u, v, z, alpha)

    #O return sera o produto intero entre o gradiente e a direcao da curva integral (que eh, no fim das contas, os campos)
    return (pc*camps[0] + sc*camps[1] + tc*camps[2])

def color_point(colors : list, Point, alpha): #Determina a cor do ponto a partir da derivada do Lambda_z
    if(Dlambdz(*Point, alpha) > 0):
        colors.append('r')
    elif(Dlambdz(*Point, alpha) < 0):
        colors.append('b')
    else:
        colors.append('white')

##############################################################
#Definindo uma constante para a funcao posterior
sqr3 = np.sqrt(3)

#Como está sendo utilizado o diagrama ternario
#Então:
def if_PointInEq(Point : list): #Retorna 1 se verdadeiro
    #Extraindo pontos
    u, v, z = Point

    #______USANDO O TEOREMA DE VIVIANI______#

    #Definindo as distancias entre o ponto P e as arestas do triangulo equilatero
    h1 = np.abs(v)
    h2 = (np.abs(-sqr3 * u + v))/2
    h3 = (np.abs(sqr3 * u + v - sqr3))/2

    #Altura total
    h = sqr3/2

    #Relizando comparacao
    if(np.abs(h1 + h2 + h3 - h) <= 1e-6):
        return 1
    
    return 0

def if_PointInRet(Point : list): #Retorna 1 se verdadeiro
    #Extraindo pontos
    u, v, z = Point

    if((0 <= u <= 1) and (0 <= v <= 1) and (0 <= u + v <= 1)):
        return 1

    return 0

#Definindo uma funcao para concatenar os pontos para h > 0 e h < 0
def Array_Concatenated(alpha, Point : list, integ_config : list):
    array_ph = nm.Euler_method(alpha, Point, integ_config) #Array dos pontos tal que h > 0
    array_mh = nm.Euler_method(alpha, Point, [-integ_config[0], integ_config[1]]) #Array dos pontos tal que h < 0

    #Retirando o ponto inicial de array_mp
    array_mh = array_mh[~np.all(array_mh == Point, axis = 1)]

    #Concatenando os arrays na ordem: h > 0 e h < 0
    return np.concatenate((array_ph, array_mh))

#Definindo um filtro para armazenar os pontos que pertencam ao dominio do prisma:
def Points_Filter(array, bar, N):
    #Definindo uma lista para armazenar os pontos ditos corretos
    Array_pc = []

    if(bar): #bar tem relacao com a funcao baricentrica 
        for i in range(2*N+1):
            if(if_PointInEq(array[i])):
                Array_pc.append(array[i])
    else:
        for i in range(2*N+1):
            if(if_PointInRet(array[i])):
                Array_pc.append(array[i])
    
    return np.array(Array_pc, float)