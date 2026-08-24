"""
Modulo: Grafico de autovalores


Objetivo:
- Plotar um único grafico;
- No grafico haverao 2 eixos: x-axis = tempo, y-axis = valores dos autovalores.
- Ainda no grafico haverao 3 curvas

- Extrair os pontos da curva dentro do prisma para encontrar todos os valores dos autovalores.

"""

from includes.Auxiliar_Functions import lambdz, transparencia, if_PointInRet, colorPoint
from includes.Functions import lbdaf, lbdas
from includes.Inicia import baricentrica
from includes.Campo_Hugoniot import sig

import numpy as np
import matplotlib.pyplot as plot
import matplotlib.lines as mlines

def pointInP(Point):
    if(not transparencia()):
        return 1

    return (if_PointInRet(Point) and (0.0 <= Point[2] <= 1.0))

def merge_continuous_branches(
        alpha,          #Variavel de controle
        branches,       #Ramos (lista de arrays [u,v,z] por ramo) - SEM pontos repetidos entre eles
        h,              #Tolerância de proximidade para decidir se dois ramos são contíguos (ESCALAR)
        sigmas = None   #Lista paralela de sigmas por ramo (Hugoniot); None para Rarefacao
):
    """Agrupa branches contíguas (por proximidade, já que não compartilham pontos)
       em segmentos maiores para plotagem contínua.
       Se sigmas for fornecido, retorna (segments, sig_segments) com a MESMA
       fusão aplicada em paralelo aos sigmas. Caso contrário, retorna só segments.
    """
    def near_Ponto(alpha, p1, p2):
        if(colorPoint(alpha, p1)[0] == colorPoint(alpha, p2)[0]):
            return True
        return False

    segments = []
    current = list(branches[0])

    if sigmas is not None:
        sig_segments = []
        current_sig = list(sigmas[0])

        for i in range(1, len(branches)):
            if near_Ponto(alpha, current[-1], branches[i][0]):
                #Contígua: como não há ponto repetido, anexa o ramo INTEIRO
                current.extend(branches[i])
                current_sig.extend(sigmas[i])
            else:
                #Descontinuidade real: fecha o segmento atual e abre um novo
                segments.append(current)
                sig_segments.append(current_sig)
                current = list(branches[i])
                current_sig = list(sigmas[i])

        segments.append(current)
        sig_segments.append(current_sig)
        return segments, sig_segments

    else:
        for i in range(1, len(branches)):
            if near_Ponto(alpha, current[-1], branches[i][0]):
                current.extend(branches[i])
            else:
                segments.append(current)
                current = list(branches[i])

        segments.append(current)
        return segments

def lambGraphPlot(ax, Point, alpha, segments, k, dist, isHugoniot):
    """Plota o grafico de autovalores em um ax ja existente."""

    #Plot do valor dos autovalores em relacao ao Ponto inicial (se estiver no prisma)
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
            ax.plot(k_branch[j:j+2], LZ[j:j+2], color='green', linestyle='-', linewidth = 1.5)

            #Lambda S -> colors[1]
            ax.plot(k_branch[j:j+2], LS[j:j+2], color='blue', linestyle = '-', markersize = 4, linewidth = 1.5)

            #Lambda F -> colors[2]
            ax.plot(k_branch[j:j+2], LF[j:j+2], color='red', linestyle='-', linewidth = 1.5)

            # if(j%61 == 0):
            #     ax.plot(k_branch[j], LZ[j], linestyle='None', marker='*', markersize=4, color='green')

            #     ax.plot(k_branch[j], LS[j], linestyle='None', marker='s', markersize=4, color='blue')
                
            #     ax.plot(k_branch[j], LF[j], linestyle='None', marker='^', markersize=4, color='red')


    #Plotagem da legenda do grafico
    legend_z = mlines.Line2D([], [], color='green', linestyle='-',  label='Lambda Z')
    legend_s = mlines.Line2D([], [], color='blue', linestyle='-', label='Lambda S')
    legend_f = mlines.Line2D([], [], color='red', linestyle='-',  label='Lambda F')
    
    handles = [legend_z, legend_s, legend_f]
    if(isHugoniot):
        legend_sig = mlines.Line2D([], [], color='k', linestyle='-', label='Sigma')
        handles.append(legend_sig)

    ax.legend(handles=handles)

    ax.grid(True)
    if(isHugoniot):
        ax.set_xlabel("Hugoniot branch parametrization")
        ax.set_ylabel("Numerical values")
    else:
        ax.set_xlabel("Integral curve parametrization")
        ax.set_ylabel("Numerical values")

    #Titulo para identificar qual o i-esimo ponto inicial que gera o grafico
    label = f'U0 = ({Point[0]:.2f}, {Point[1]:.2f}, {Point[2]:.2f}) with alpha = {alpha}'
    ax.set_title(label) 

def sigGraphPlot(
        ax, 
        Point, 
        alpha, 
        sig_segments,
        k,
        dist
):
    if dist is not None:
        sigmaPoint = sig(alpha, Point[2], *Point)
    
        ax.plot(dist, sigmaPoint, color = 'k', marker = 'o', zorder = 3) #Formato (U0, sig(alpha, U0))

    k_offset = 0

    for i, sig_branch in enumerate(sig_segments):
        n = len(sig_branch)
        k_branch = k[k_offset : k_offset + n]
        k_offset += n

        #Garante que os valores de sigma sejam floats simples (1D)
        sig_branch = np.array(sig_branch, dtype=float)

        # Plota o ramo de sigma de uma so vez (muito mais rapido e sem erros de fatiamento)
        ax.plot(k_branch, sig_branch, color='k', linestyle='--', linewidth=2.25, dashes=(10, 10))

        #PLotagem dos marcadores
        #ax.plot(k_branch[::47], sig_branch[::47], linestyle='None', marker='x', markersize=4, color= 'k')

def lamb_graph(
        alpha        : float,
        Point        : list,
        integ_config : list,
        branches     : list,
        isHugoniot = False
):
    h, N = integ_config   #Desempacota o passo escalar 'h' usado na fusao de branches

    def ponto_igual(p1, p2, tol=1e-9):
        return np.allclose(p1, p2, atol=tol) 

    sigs = None
    if(isHugoniot):
        points_branches = []
        sig_branches = []
        for branch in branches:
            points = np.array([p[0] for p in branch])
            sigmas = np.array([p[1] for p in branch])
            points_branches.append(points)
            sig_branches.append(sigmas)
 
        branches = points_branches
        sigs = sig_branches

    segments = None
    seg_sigs = None
    if(isHugoniot):
        segments, seg_sigs = merge_continuous_branches(alpha, branches, h, sigs)
    else:
        segments = merge_continuous_branches(alpha, branches, h)

    #Conta exatamente quantos pontos restaram apos a limpeza
    total_pontos = sum(len(branch) for branch in branches)

    #Cria o vetor k variando EXATAMENTE de 0 a 1 distribuído igualmente (uma lista do intervalo)
    k = np.linspace(0, 1, total_pontos).tolist()   

    dist = None
    if pointInP(Point) or (not transparencia):
        indice_absoluto = 0
        ponto_encontrado = False
        
        for branch in branches:
            for p in branch:
                if ponto_igual(p, Point):
                    ponto_encontrado = True
                    break
                indice_absoluto += 1
            if ponto_encontrado:
                break

        if ponto_encontrado:
            dist = indice_absoluto / (total_pontos - 1)
        #se nao encontrado, dist continua None -> nao plota o marcador

    fig, ax = plot.subplots()
    lambGraphPlot(ax, Point, alpha, segments, k, dist, isHugoniot)

    if(isHugoniot):
        sigGraphPlot(ax, Point, alpha, seg_sigs, k, dist)
    
    #Plotando todas as figuras de uma so vez
    #plot.show()
