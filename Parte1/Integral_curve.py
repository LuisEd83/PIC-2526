"""
Modulo: Curva Integral

Objetivos:
-> Plotar curvas integrais usando o metodo de euler
"""

import includes._Branches as b
import includes.Auxiliar_Functions as af
import includes.Functions as fun

import Parte1.LambdaZ_vec as lv
import matplotlib.pyplot as plot

#Possui as branches como parametro 
#Eh utilizada no Points_prism.py
def plotIntegralCurve(ax, Point, alpha, branches, colors):
    for i in range(len(branches)): #Esse laco irah criar as conexoes entre os pontos
        ax.plot(branches[i][:, 0], #Conjunto da componente X da i-esima branch 
                branches[i][:, 1], #Conjunto da componente Y da i-esima branch
                branches[i][:, 2], #Conjunto da componente Z da i-esima branch
                color = colors[i][0], #Cor da i-iesima branch
                linestyle = '-')   #Tipo de linha (no caso sera a linha continua)
                #zorder = 4)        #Ordem de prioridade

    
    if(alpha != 0):
        lv.plotVecsLambdaZWB(ax, alpha, branches, 5) #Plotagem dos vetores
    #ax.plot(Point[0], Point[1], Point[2], marker = '.', color = 'k', zorder = 5)
    ax.scatter(*Point, color='black', s=40, depthshade=False)

#Nao possui as branches como parametro (Calculo interno)
def integralCurve(alpha, configCurve):
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

    #Variavel que armazena os pontos iniciais
    iniPoints = af.read_points()

    for Point in iniPoints: 
        branches = b.Branches_point(alpha, Point, configCurve) #Branches
        colors = b.Branches_point_colors(alpha, branches) #Cores das branches

        plotIntegralCurve(ax, Point, alpha, branches, colors)

    plot.show()

#integralCurve(0.1, [0.01, 500])
#integralCurve(0.01, [0.6, 0.3, 0.2], [0.01, 500])