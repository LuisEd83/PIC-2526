"""
Modulo: Vetores do Lambda_Z

Objetivos:
-> Criar um algoritmo capaz de plotar vetores que indicam o sentido de crescimento do Lambda_Z
na curva integral.
"""

import includes._Branches as b
import includes.Functions as fun
import includes.Inicia as ini

import numpy as np

def plotVecsLambdaZWB(ax, alpha, branches, delta):
    import includes.Auxiliar_Functions as af

    def is_continuous(p_last, p_first, tol=1e-2):
        return np.allclose(p_last, p_first, atol=tol)

    #Agrupa branches contínuas em segmentos
    segments = []
    current_segment = list(branches[0])

    for i in range(len(branches) - 1):
        if is_continuous(branches[i][-1], branches[i+1][0]):
            #Continua: concatena a proxima branch no segmento atual
            current_segment.extend(branches[i+1])
        else:
            #Descontinuidade: salva o segmento atual e começa um novo
            segments.append(current_segment)
            current_segment = list(branches[i+1])
    
    segments.append(current_segment)  #adiciona o ultimo segmento

    #Plota os vetores para cada segmento isoladamente
    for segment in segments:
        flat    = af.remove_consecutive_duplicates(segment)
        sampled = af.stride_sample_symmetric(flat, delta)

        for i in range(len(sampled) - 1):
            if (not(i % 2)):
                d = np.array(sampled[i+1]) - np.array(sampled[i])
                if af.lambdz(*sampled[i], alpha) > af.lambdz(*sampled[i+1], alpha):
                    af.plotCone(ax, sampled[i], -d, af.colorPoint(alpha, sampled[i]))
                else:
                    af.plotCone(ax, sampled[i], d, af.colorPoint(alpha, sampled[i]))

def plotVecsLambdaZ(alpha, Point, integConfig, delta):
    import includes.Auxiliar_Functions as af
    import matplotlib.pyplot as plt

    def is_continuous(p_last, p_first, tol=1e-2):
        return np.allclose(p_last, p_first, atol=tol)

    #__________IICIALIZANDO AMBIENTE 3D__________#

    #mapeia para baricentrica ou não
    mp = ini.baricentrica()

    #Define as concentracoes minima e maxima
    zmin, zmax = ini.concentrations() 

    #Inicializando figura
    ax = ini.ambiente3d()
    ax.view_init(elev=30., azim=-130.) #Initial Camera Position

    if(af.transparencia()):
        #__________INICIALIZANDO PLOTAGEM DO PRISMA__________#
        
        #Vertices do triangulo
        Gw, Go = 0, 0
        Ww, Wo = 1, 0
        Ow, Oo = 0, 1

        # Mapeia para coordenadas baricentricas de mp = 1
        G1, G2 = fun.map(Gw, Go, mp)
        W1, W2 = fun.map(Ww, Wo, mp)
        O1, O2 = fun.map(Ow, Oo, mp)

        # Triangulo no plano cmin
        ax.plot([G1,W1], [G2, W2], [zmin,zmin], 'k')
        ax.plot([G1,O1], [G2, O2], [zmin,zmin], 'k')
        ax.plot([W1,O1], [W2, O2], [zmin,zmin], 'k')

        #Triangulo no plano c = cmax
        ax.plot([G1,W1], [G2, W2], [zmax,zmax], 'k')
        ax.plot([G1,O1], [G2, O2], [zmax,zmax], 'k')
        ax.plot([W1,O1], [W2, O2], [zmax,zmax], 'k')

        #Arestas verticais
        ax.plot([G1,G1], [G2,G2], [zmin,zmax], 'k')
        ax.plot([W1,W1], [W2,W2], [zmin,zmax], 'k')
        ax.plot([O1,O1], [O2,O2], [zmin,zmax], 'k')

        #Identificacao dos eixos
        ax.set_xlabel('$u$')
        ax.set_ylabel('$v$')
        ax.set_zlabel('$z$')

    #Inicializa as branches
    branches = b.Branches_point(alpha, Point, integConfig)

    #Realiza o plot dos vetores sobre a curva
    #Agrupa branches contínuas em segmentos
    segments = []
    current_segment = list(branches[0])

    for i in range(len(branches) - 1):
        if is_continuous(branches[i][-1], branches[i+1][0]):
            #Continua: concatena a proxima branch no segmento atual
            current_segment.extend(branches[i+1])
        else:
            #Descontinuidade: salva o segmento atual e começa um novo
            segments.append(current_segment)
            current_segment = list(branches[i+1])
    
    segments.append(current_segment)  #adiciona o ultimo segmento

    #Plota os vetores para cada segmento isoladamente
    for segment in segments:
        flat    = af.remove_consecutive_duplicates(segment)
        sampled = af.stride_sample_symmetric(flat, delta)

        for i in range(len(sampled) - 1):
            if (not(i % 2)):
                d = np.array(sampled[i+1]) - np.array(sampled[i])
                if af.lambdz(*sampled[i], alpha) > af.lambdz(*sampled[i+1], alpha):
                    af.plotCone(ax, sampled[i], -d, af.colorPoint(alpha, sampled[i]))
                else:
                    af.plotCone(ax, sampled[i], d, af.colorPoint(alpha, sampled[i]))

    plt.show()
