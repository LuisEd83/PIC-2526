"""
Modulo: Grafico de autovalores


Objetivo:
- Plotar um único grafico;
- No grafico haverao 2 eixos: x-axis = tempo, y-axis = valores dos autovalores.
- Ainda no grafico haverao 3 curvas

- Extrair os pontos da curva dentro do prisma para encontrar todos os valores dos autovalores.

"""

from includes.Auxiliar_Functions import lambdz, transparencia, if_PointInEq, if_PointInRet, read_points, colorPoint
from includes._Branches import Branches_point
from includes.Functions import lbdaf, lbdas
from includes.Inicia import baricentrica

import numpy as np
import matplotlib.pyplot as plot
import matplotlib.lines as mlines

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

    #Plot do valor dos autovalores em relação ao Ponto inicial (se estiver no prisma)
    if dist is not None:
        lambZPoint = lambdz(*Point, alpha)
        lambSPoint = lbdas(*Point)
        lambFPoint = lbdaf(*Point)

        ax.plot(dist, lambZPoint, color = 'k', marker = 'o', zorder = 3) #Formato (U0, lamb(U0))
        ax.plot(dist, lambSPoint, color = 'k', marker = 'o', zorder = 3) #Formato (U0, lamb(U0))
        ax.plot(dist, lambFPoint, color = 'k', marker = 'o', zorder = 3) #Formato (U0, lamb(U0))

    k_offset = 0

    for i, segment in enumerate(segments):
        n = len(segment)
        k_branch = k[k_offset : k_offset + n]
        k_offset += n

        LS = [lbdas(*p) for p in segment]
        LF = [lbdaf(*p) for p in segment]
        LZ = [lambdz(*p, alpha) for p in segment]

        colors = [colorPoint(alpha, p) for p in segment]
        #formato de colors => [color_LZ, color_LS, color_LF]

        for j in range(len(segment) - 1):
            color_j = colors[j]

            #Lambda Z -> colors[0]
            ax.plot(k_branch[j:j+2], LZ[j:j+2], color=color_j[0], linestyle='-', linewidth = 1)

            #Lambda S -> colors[1]
            ax.plot(k_branch[j:j+2], LS[j:j+2], color=color_j[1], linestyle = '--', marker = '*', markersize = 4, linewidth = 1)

            #Lambda F -> colors[2]
            ax.plot(k_branch[j:j+2], LF[j:j+2], color=color_j[2], linestyle=':', linewidth = 1)


    #Plotagem da legenda do gráfico
    legend_z = mlines.Line2D([], [], color='k', linestyle='-',  label='Lambda Z')
    legend_s = mlines.Line2D([], [], color='k', linestyle='--', marker = '*', label='Lambda S')
    legend_f = mlines.Line2D([], [], color='k', linestyle='--',  label='Lambda F')

    ax.legend(handles=[legend_z, legend_s, legend_f])

    ax.grid(True)
    ax.set_xlabel("Eixo X - integral curve parametrization")
    ax.set_ylabel("Eixo Y - Numerical values")

    #Titulo para identificar qual o i-esimo ponto inicial que gera o grafico
    label = f'P0 = ({Point[0]:.2f}, {Point[1]:.2f}, {Point[2]:.2f}) with alpha = {alpha}'
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
                

        #Conta exatamente quantos pontos restaram apos a limpeza
        total_pontos = sum(len(branchei) for branchei in branches)

        #Cria o vetor k variando EXATAMENTE de 0 a 1 distribuído igualmente (uma lista do intervalo)
        k = np.linspace(0, 1, total_pontos).tolist()

        dist = None
        if pointInP(Point) or (not transparencia):
            #Variavel para armazenar o indice absoluto do ponto na malha sequencial
            indice_absoluto = 0
            ponto_encontrado = False
            
            for branch in branches:
                for p in branch:
                    if ponto_igual(p, Point):
                        ponto_encontrado = True
                        break
                    indice_absoluto += 1 #Vai incrementando enquanto nao acha
                if ponto_encontrado:
                    break

            #Mapeia o indice do ponto exatamente na mesma escala do vetor k
            #Se o total de pontos é N, os indices vão de 0 a N-1
            dist = indice_absoluto / (total_pontos - 1)

        fig, ax = plot.subplots()
        lambGraphPlot(ax, Point, alpha, segments, k, dist)
    
    #Plotando todas as figuras de uma so vez
    plot.show()

