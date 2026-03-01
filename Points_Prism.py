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
import Inicia as ini
import Functions as fun

import numpy as np
import matplotlib.pyplot as plt

#Definindo constantes:
h = 0.01 #Passo do metodo de Euler
N = 1000 #Numero de realizacao de passo

#Definindo o ponto inicial para teste
U0 = [0.5, 0.2, 0.1]

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
    array_mh = array_mh[~np.all(array_mh == Point, axis = 1)]
    
    print(f"PH shape: {array_ph.shape}")
    print(f"MH shape: {array_mh.shape}")

    #Concatenando os arrays na ordem: h > 0 e h < 0
    return np.concatenate((array_ph, array_mh))

#Definindo um filtro para armazenar os pontos que pertencam ao dominio do prisma:
def Points_Filter(array):
    #Definindo uma lista para armazenar os pontos ditos corretos
    Array_pc = []

    for i in range(2*N+1):
        if(if_PointInPrism(array[i])):
            Array_pc.append(array[i])
    
    return np.array(Array_pc, float)

#Definindo a funcao principal de plot
def Prism_plot(U0 : list):
    #__________INICIALIZANDO PLOTAGEM DO PRISMA__________#

    #mapeia para baricentrica ou não
    mp = ini.baricentrica()

    #Define as concentracoes minima e maxima
    zmin, zmax = ini.concentrations() 

    #Inicializando figura
    ax = ini.ambiente3d()
    ax.view_init(elev=30., azim=-130.) #Initial Camera Position
    
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

    #__________INICIALIZANDO PLOTAGEM DA CURVA__________#
    
    #Colecao de pontos dentro do prisma:
    points_concatened = Array_Concatenated(U0) #Extrai todos os pontos
    array_points = Points_Filter(points_concatened) #Filtra os pontos

    #plot dos pontos
    for i in range(len(array_points)):
        ax.plot(*array_points[i], marker = ".", linestyle = "-", color = "b")

    #Inicia plotagem
    plt.show()

Prism_plot(U0)