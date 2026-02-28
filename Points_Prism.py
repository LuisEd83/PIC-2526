"""
Modulo: plotagem

Objetivos:
- Plotar os pontos no interior do dominio do prisma;
- Concatenar os pontos com h > 0 e h < 0;
- Deletar os pontos que fiquem no exterior do prisma.
"""

#Bibliotecas  
import Euler_Integration as ei
import Determinants_Functions as df

import numpy as np
import matplotlib.pyplot as plot

#Definindo constantes:
h = 0.001 #Passo do metodo de Euler
N = 1000 #Numero de realizacao de passo

#Como o prisma é definido por:
#P = {0 ≤ u ≤ 1, 0 ≤ v ≤ 1, 0 ≤ u + v ≤ 1, 0 ≤ z ≤ 1}
#Então:
def if_PointInPrism(Point : list):
    #Extraindo pontos
    u, v, z = Point

    #Relizando comparacao
    if((0.0 <= u <= 1.0) and (0.0 <= v <= 1.0) and (0.0 <= u + v <= 1.0) and (0.0 <= z <= 1.0)):
        return 1
    
    return 0

#Definindo uma funcao para concatenar os pontos para h > 0 e h < 0
def Array_Concatenated(Point : list):
    array_ph = ei.Euler_method(Point, [h, N]) #Array dos pontos tal que h > 0
    array_mh = ei.Euler_method(Point, [-h, N]) #Array dos pontos tal que h < 0

    #Retirando o ponto inicial de array_mp
    array_mh = np.delete(array_mh, Point)

    #Concatenando os arrays na ordem: h > 0 e h < 0
    return np.concatenate((array_ph, array_mh))

#Definindo um filtro para armazenar os pontos que pertencam ao dominio do prisma:
def Points_Filter(array : np.array):
    #Definindo uma lista para armazenar os pontos ditos corretos
    Array_pc = []

    for i in range(2*N+1):
        if(if_PointInPrism(array[i])):
            Array_pc.append(array[i])
    
    return np.array(Array_pc, float)