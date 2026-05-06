"""
Modulo: Grafico de autovalores


Objetivo:
- Plotar um único grafico;
- No grafico haverao 2 eixos: x-axis = tempo, y-axis = valores dos autovalores.
- Ainda no grafico haverao 3 curvas

- Extrair os pontos da curva dentro do prisma para encontrar todos os valores dos autovalores.

"""

from includes.Auxiliar_Functions import lambdz, transparencia, if_PointInEq, if_PointInRet, read_points
from includes._Branches import Branches_point
from includes.Functions import lbdaf, lbdas
from includes.Inicia import baricentrica

import numpy as np
import matplotlib.pyplot as plot

def pointInP(Point):
    if(not transparencia()):
        return 1
    if(baricentrica()):
        return (if_PointInEq(Point) and (0.0 <= Point[2] <= 1.0))
    else:
        return (if_PointInRet(Point) and (0.0 <= Point[2] <= 1.0))

def merge_continuous_branches(branches, tol=1e-9):
    """Agrupa branches contíguas em segmentos maiores para plotagem contínua."""
    def ponto_igual(p1, p2):
        return np.allclose(p1, p2, atol=tol)

    segments = []
    current = list(branches[0])

    for i in range(1, len(branches)):
        if ponto_igual(current[-1], branches[i][0]):
            #Contínua: anexa (sem repetir o ponto de junção)
            current.extend(branches[i][1:])
        else:
            #Descontinuidade: fecha o segmento atual e abre um novo
            segments.append(current)
            current = list(branches[i])

    segments.append(current)  #fecha o último segmento
    return segments

def lambGraphPlot(ax, Point, alpha, segments, k, dist):
    """Plota o grafico de autovalores em um ax ja existente."""

    if dist is not None:
        ax.axvline(x=dist, color='k', linestyle='--', label='Ponto inicial')

    k_offset = 0
    for i, segment in enumerate(segments):
        n = len(segment)
        k_branch = k[k_offset : k_offset + n]
        k_offset += n

        LS = [lbdas(*p) for p in segment]
        LF = [lbdaf(*p) for p in segment]
        LZ = [lambdz(*p, alpha) for p in segment]

        label_s = 'Lambda S' if i == 0 else '_nolegend_'
        label_f = 'Lambda F' if i == 0 else '_nolegend_'
        label_z = 'Lambda Z' if i == 0 else '_nolegend_'

        ax.plot(k_branch, LS, color='b',       linestyle='-', label=label_s)
        ax.plot(k_branch, LF, color='r',       linestyle='-', label=label_f)
        ax.plot(k_branch, LZ, color='magenta', linestyle='-', label=label_z)

    ax.grid(True)
    ax.legend()
    ax.set_xlabel("Eixo X - Passos")
    ax.set_ylabel("Eixo Y - Valores numéricos")

    #Titulo para identificar qual o i-esimo ponto inicial que gera o grafico
    label = f'P0 = ({Point[0]:.2f}, {Point[1]:.2f}, {Point[2]:.2f})'
    ax.set_title(label) 

def lamb_graph(alpha, integ_config : list):
    def ponto_igual(p1, p2, tol=1e-9):
        return np.allclose(p1, p2, atol=tol) 

    #Variavel que armazena os pontos iniciais
    iniPoints = read_points()

    for Point in iniPoints:
        #Extraio as branches:
        branches = Branches_point(alpha, Point, integ_config) #Ramos (conjunto de pontos)
        segments = merge_continuous_branches(branches)        #Branches de ramos contínuos

        #Limpando os pontos iguais nas branches:
        j = 1
        for i in range(len(branches)):                          #Varre as branches
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

        fig, ax = plot.subplots()
        lambGraphPlot(ax, Point, alpha, segments, k, dist)
    
    #Plotando todas as figuras de uma so vez
    plot.show()