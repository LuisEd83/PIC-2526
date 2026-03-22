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

from Auxiliar_Functions import lambdz, transparencia, if_PointInEq, if_PointInRet
from _Branches import Branches_point
from Functions import lbdaf, lbdas
from Inicia import baricentrica

import numpy as np
import matplotlib.pyplot as plot

def pointInP(Point):
    if(not transparencia()):
        return 1
    if(baricentrica()):
        return (if_PointInEq(Point) and (0.0 <= Point[2] <= 1.0))
    else:
        return (if_PointInRet(Point) and (0.0 <= Point[2] <= 1.0))

def lamb_graph(alpha, Point : list, integ_config : list):
    #Importando biblioteca para nao haver erro durante a lida do arquivo Functions.py
    import Functions as fun
    
    def ponto_igual(p1, p2, tol=1e-9):
        return np.allclose(p1, p2, atol=tol) 

    #Extraio as branches:
    branches = Branches_point(alpha, Point, integ_config) #Ramos (conjunto de pontos)

    #Limpando os pontos iguais nas branches:
    j = 1
    for i in range(len(branches)):                          #Varre as brances
        if(j == len(branches)):
            break
        if(ponto_igual(branches[i][-1], branches[j][0])):
             branches[j] = branches[j][1:]
        j += 1
            

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

    if(pointInP(Point) or (not transparencia)):
        #Variavel para armazenar a branch e o indice na branch
        local_point = []
        for i, branch in enumerate(branches):
            for j, p in enumerate(branch):
                if(ponto_igual(p, Point)):     #Localiza o ponto dentre as branches
                    #i = branch; j = indice na branch i
                    local_point.append(i)
                    local_point.append(j)

        dist = 0  #Variavel de distancia do primeiro ponto ate o Point 
        if(local_point[0] == 0):
            dist += local_point[1]
        else:
            for i in range(local_point[0]):
                dist += len(branches[i])
            dist += local_point[1] #Corrige a distancia 

        dist /= tam   #Redimensiona a distancia para o intervalo [0, 1]

    for i in range(len(branches)):            #Varre os indices das branches
        for j in range(len(branches[i])):     #Varre a i-esima branch
            Lambda_S.append(lbdas(*branches[i][j]))
            Lambda_F.append(lbdaf(*branches[i][j]))
            Lambda_Z.append(lambdz(*branches[i][j], alpha))
    
    #Inicializo a figura
    plot.figure()

    if(pointInP(Point) or (not transparencia)):
        plot.axvline(x = dist, color = 'k', linestyle = '--', label = 'Ponto inicial')  #Plotagem da localização do Point

    plot.plot(k, Lambda_S, color = 'b', linestyle = '-', label = 'Lambda S')        #Plotagem das variaveis relacionadas a Lambda_S
    plot.plot(k, Lambda_F, color = 'r', linestyle = '-', label = 'Lambda F')        #Plotagem das variaveis relacionadas a Lambda_F
    plot.plot(k, Lambda_Z, color = 'magenta', linestyle = '-', label = 'Lambda Z')  #Plotagem das variaveis relacionadas a Lambda_Z

    plot.grid(True)
    plot.legend()
    plot.xlabel("Eixo X - Passos")
    plot.ylabel("Eixo Y - Valores numéricos")
    plot.show()