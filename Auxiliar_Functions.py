"""
Modulo: Funcoes auxiliares

Objetivo:
- Armazenar as funcoes que serao utilizadas no codigo e selecionar as necessarias

"""

import Euler_Integration as ei
import Functions as fun

import numpy as np

#Derivada da funcao continua e crescente no intervalo [0,1]
def az(z): #Da/dz
    return np.cos(z) #Derivada de a(z) em relacao a z

def lambdz(u, v, z, alpha):
    return(fun.fw(u, v, z)/(u + alpha*az(z)))

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
    if((h1 + h2 + h3 - h == 0) and (0.0 <= z <= 1.0)):
        return 1
    
    return 0

def if_PointInRet(Point : list): #Retorna 1 se verdadeiro
    #Extraindo pontos
    u, v, z = Point

    if((0 <= u <= 1) and (0 <= v <= 1) and (0 <= u + v <= 1) and (0 <= z <= 1)):
        return 1

    return 0

#Definindo uma funcao para concatenar os pontos para h > 0 e h < 0
def Array_Concatenated(alpha, Point : list, integ_config : list):
    array_ph = ei.Euler_method(alpha, Point, integ_config) #Array dos pontos tal que h > 0
    array_mh = ei.Euler_method(alpha, Point, [-integ_config[0], integ_config[1]]) #Array dos pontos tal que h < 0

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