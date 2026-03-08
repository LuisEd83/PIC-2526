"""
Modulo: plotagem

Objetivos:
- Plotar os pontos no interior do dominio do prisma;
- Concatenar os pontos com h > 0 e h < 0;
- Deletar os pontos que fiquem no exterior do prisma.
"""

#Bibliotecas
import Auxiliar_Functions as af
import Functions as fun

import numpy as np
import matplotlib.pyplot as plt

#Definindo constantes:
h = 0.01 #Passo do metodo de Euler
N = 1000 #Numero de realizacao de passo
alpha = 0.3 #Variável de controle

sqr3 = np.sqrt(3) 

#Definindo o ponto inicial para teste
U0 = [0.3, 0.2, 0.1]
U1 = [0.5, 0.15, 0.1]

#Definindo a funcao principal de plot
def Prism_plot(Point : list):
    #Importando biblioteca para nao haver erro durante a lida do arquivo Functions.py
    import Inicia as ini

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
    points_concatened_C = af.Array_Concatenated(alpha, Point, [h, N]) #Extrai todos os pontos
    array_points_c = af.Points_Filter(points_concatened_C, mp, N) #Filtra os pontos

    #plot dos pontos
    for i in range(len(array_points_c)):
        ax.plot(*array_points_c[i], 'b.-')

    #Inicia plotagem
    plt.show()

Prism_plot(U1)