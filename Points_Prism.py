"""
Modulo: plotagem

Objetivos:
- Plotar os pontos no interior do dominio do prisma;
- Concatenar os pontos com h > 0 e h < 0;
- Deletar os pontos que fiquem no exterior do prisma.
"""

#Bibliotecas
import Functions as fun
import _Branches as b
import Auxiliar_Functions as af

import numpy as np
import matplotlib.pyplot as plt

#Definindo constantes:
h = 0.01        #Passo do metodo de Euler
N = 1000        #Numero de realizacao de passo
alpha = 0.01    #Variável de controle
Num_z = 20      #Numero de curvas de nivel

sqr3 = np.sqrt(3) 

#Definindo o ponto inicial para teste
U0 = [0.3, 0.2, 0.1]

U1 = [0.5, 0.1, 0.1] #Ponto importante

U2 = [0.2, 0.4, 0.1] #Ponto inicial fora do prisma

def if_pointIn(bar, Point : list):
    if(bar):
        return af.if_PointInEq(Point)
    return af.if_PointInRet(Point)

#Definindo a funcao principal de plot
def Prism_plot(Point : list):
    #Importando biblioteca para nao haver erro durante a lida do arquivo Functions.py
    import Inicia as ini

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

    #__________INICIALIZANDO PLOTAGEM DA CURVA RELACIONADA AOS CAMPOS__________#
    
    branches = b.Branches_point(alpha, Point, [h, N]) #Branches
    colors = b.Branches_point_colors(alpha, branches) #Cores das branches

    for i in range(len(branches)): #Esse laco irah criar as conexoes entre os pontos
        ax.plot(branches[i][:, 0], #Conjunto da componente X da i-esima branch 
                branches[i][:, 1], #Conjunto da componente Y da i-esima branch
                branches[i][:, 2], #Conjunto da componente Z da i-esima branch
                color = colors[i], #Cor da i-iesima branch
                linestyle = '-')   #Tipo de linha (no caso sera a linha continua)

    #Plotando o ponto inicial
    ax.plot(*Point, 'ko')

    #__________INICIALIZANDO PLOTAGEM DA CURVA RELACIONADA AOS AUTOVALORES__________#
    auto_branches = b.Branches_auto([h, N], Num_z)
    auto_colors = b.Branches_auto_colors()

    for i in range(len(auto_branches)):
        for j in range(len(auto_branches[i])):  # ← já deve estar assim
            seg = np.array(auto_branches[i][j])
            ax.plot(seg[:, 0],
                    seg[:, 1],
                    seg[:, 2],
                    color = auto_colors[j],
                    linestyle = '-')


    #Inicia plotagem
    plt.show()

Prism_plot(U2)