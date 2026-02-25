"""
Modulo: Integracao pelo metodo de Euler

Objetivo:
- Utilizar o metodo de Euler para integrar os campos a partir de um ponto escolhido
"""

import Determinants_Functions as df
import numpy as np


#Definindo uma funcao que implementa o metodo de Euler para integracao
def Euler_method(point : list, integ_config : list):
    #Extraindo ponto inicial
    u0 = point[0]
    v0 = point[1]
    z0 = point[2]

    #Extraindo configuracao de integracao
    h = integ_config[0]
    N = integ_config[1]
