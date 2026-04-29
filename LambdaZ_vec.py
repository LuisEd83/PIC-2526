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

def plotVec(ax, alpha, p0, p1): #Suponto p1 o sucessor de p0
    import includes.Auxiliar_Functions as af 
    du, dv, dz = p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]   #deslocamentos (vetor que aponta de p0 para p1)
    if(af.lambdz(*p1, alpha) > af.lambdz(*p0, alpha)):
        ax.quiver(p0[0], p0[1], p0[2], du, dv, dz, color = "r")
    else:
        ax.quiver(p1[0], p1[1], p1[2], -du, -dv, -dz, color = "b")

def stride_sample_gen(points, step): #Seleciona os pontos com uma distancia (na lista) fixa
    for i in range(0, len(points), step):
        yield points[i]

def remove_consecutive_duplicates(points, tol=1e-9):
    unique = [points[0]]
    for p in points[1:]:
        if not np.allclose(p, unique[-1], atol=tol):
            unique.append(p)
    return unique

def plotVecsLambdaZWB(ax, alpha, branches, delta):
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
        flat = remove_consecutive_duplicates(segment)
        sampled = list(stride_sample_gen(flat, delta))

        for i in range(len(sampled) - 1):
            plotVec(ax, alpha, sampled[i], sampled[i+1])

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
        flat = remove_consecutive_duplicates(segment)
        sampled = list(stride_sample_gen(flat, delta))

        for i in range(len(sampled) - 1):
            plotVec(ax, alpha, sampled[i], sampled[i+1])
    
    plt.show()
