"""
Modulo: Grafico de autovalores


Objetivo:
- Plotar um único grafico;
- No grafico haverao 2 eixos: x-axis = tempo, y-axis = valores dos autovalores.
- Ainda no grafico haverao 3 curvas

- Extrair os pontos da curva dentro do prisma para encontrar todos os valores dos autovalores.

"""

"""
1) Flags: Plotar todos os pontos // 
"""

from Auxiliar_Functions import lambdz, transparencia
from _Branches import Branches_point, Branches_point_colors
from Functions import lbdaf, lbdas

import numpy as np
import matplotlib.pyplot as plot

def lamb_graph(alpha, Point : list, integ_config : list):
    #Importando biblioteca para nao haver erro durante a lida do arquivo Functions.py
    import Functions as fun
    
    #Extraio as branches:
    branches = Branches_point(alpha, Point, integ_config) #Ramos (conjunto de pontos)

    #Variavel de escala (0 a 1)
    k = []
    c = 0  #Variavel de correcao de erro
    tam = sum(len(branchei) for branchei in branches) #Tamanho total das branches
    for i in range(len(branches)):
        for _ in range(len(branches[i])):
            c += 1
            k.append(c/tam)                           #Reescala para o intervalo [0,1]

    #Variaveis que irao armazenar os pontos
    Lambda_S = []
    Lambda_F = []
    Lambda_Z = []

    for i in range(len(branches)):            #Varre os indices das branches
            for j in range(len(branches[i])): #Varre a i-esima branch
                Lambda_S.append(lbdas(*branches[i][j]))
                Lambda_F.append(lbdaf(*branches[i][j]))
                Lambda_Z.append(lambdz(*branches[i][j], alpha))
    
    #Inicializo a figura
    plot.figure()

    plot.plot(k, Lambda_S, color = 'b', linestyle = '-', label = 'Lambda S')        #Plotagem das variaveis
    plot.plot(k, Lambda_F, color = 'r', linestyle = '-', label = 'Lambda F')        #Plotagem das variaveis
    plot.plot(k, Lambda_Z, color = 'magenta', linestyle = '-', label = 'Lambda S')  #Plotagem das variaveis

    plot.grid(True)
    plot.legend()
    plot.xlabel("Eixo X - Passos")
    plot.ylabel("Eixo Y - Valores numéricos")
    plot.show()