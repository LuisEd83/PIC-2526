"""
Modulo: plotagem

Objetivos:
- Plotar os pontos no interior do dominio do prisma;
- Concatenar os pontos com h > 0 e h < 0;
- Deletar os pontos que fiquem no exterior do prisma.
"""

#Bibliotecas
import includes.Functions as fun
import includes._Branches as b
import includes.Auxiliar_Functions as af 
import Autovalues_graph as ag
import Coincidence_surfaces as cs
import Integral_curve as ic

import numpy as np
import matplotlib.pyplot as plt

#Definindo constantes:
h = 0.01        #Passo do metodo de Euler
N = 500        #Numero de realizacao de passo
alpha = 0.0    #Variável de controle
Num_z = 20      #Numero de curvas de nivel

#Definindo a funcao principal de plot
def main():
    #Importando biblioteca para nao haver erro durante a lida do arquivo Functions.py
    import includes.Inicia as ini

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

    #Variavel que armazena os pontos iniciais:
    iniPoints = af.read_points()

    #Itera sobre cada ponto em iniPoints
    for Point in iniPoints:
        #Variavel que armazena as branches
        branches = b.Branches_point(alpha, Point, [h, N])
        colors = b.Branches_point_colors(alpha, branches)

        #__________INICIALIZANDO PLOTAGEM DA CURVA INTEGRAL__________#
        ic.plotIntegralCurve(ax, Point, alpha, branches, colors)
        
    #__________INICIALIZANDO PLOTAGEM DE GRAFICO DAS VELOCIDADES__________#
    ag.lamb_graph(alpha, [h, N])
    
    #__________INICIALIZANDO PLOTAGEM DA SUPERFICIE DE COINCIDENCIA__________#
    base = cs.argmin_branch(branches)
    auto_branches = b.Branches_auto([h, N], Num_z, alpha, base - 0.5, altura = 1)
    auto_colors = b.Branches_auto_colors()

    cs.plotCoicindenceCurves(ax, auto_branches, auto_colors)

    #Inicia plotagem
    plt.show()

main()