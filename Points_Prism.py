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
    #Extraindo ponto
    u, v, z = Point

    #Relizando comparacao
    if((0.0 <= u <= 1.0) and (0.0 <= v <= 1.0) and (0.0 <= u + v <= 1.0) and (0.0 <= z <= 1.0)):
        return 1
    
    return 0

#Definindo uma funcao para concatenar os pontos para h > 0 e h < 0
def Array_Concatenated(Point : list):
    array_ph = ei.Euler_method(Point, [h, N])
    array_mh = ei.Euler_method(Point, [-h, N])

    #Concatenando os arrays na ordem: h > 0 e h < 0
    return np.concatenate((array_ph, array_mh))