"""
Modulo: Curva Integral

Objetivos:
-> Plotar curvas integrais usando o metodo de euler
"""

import includes._Branches as b
import includes.Auxiliar_Functions as af
import includes.Functions as fun

import LambdaZ_vec as lv
import matplotlib.pyplot as plt

#Funcao que desenha um triangulo
def triangulo(
        ax,    #ax para plotagem
        mp     #Baricentrica
):
    #Vertices do triangulo
    Gw, Go = 0, 0
    Ww, Wo = 1, 0
    Ow, Oo = 0, 1

    G1, G2 = fun.map(Gw, Go, mp)
    W1, W2 = fun.map(Ww, Wo, mp)
    O1, O2 = fun.map(Ow, Oo, mp)

    ##Ponto umbilico no plano cl
    ## Precisa inicializar as viscosidades
    #Uwl, Uol = muwl/(muwl + muo + mug), muwl/(muwl + muo + mug)
    #Uwl, Uol = fun.map(Uwl, Uol, mp)

    ####################################################

    #O triangulo
    ax.plot([G1, W1], [G2, W2], 'k-') 
    ax.plot([G1, O1], [G2, O2], 'k-') 
    ax.plot([W1, O1], [W2, O2], 'k-') 

#Possui as branches como parametro 
#Eh utilizada no Points_prism.py
def plotIntegralCurve(ax, Point, alpha, branches, colors):
    from includes.Inicia import baricentrica

    if(baricentrica()):
        for branch in branches:
            for p in branch:
                u, v, z = p
                u, v = fun.map(u, v, baricentrica())
                p[0] = u
                p[1] = v

    for i in range(len(branches)): #Esse laco irah criar as conexoes entre os pontos
        ax.plot(branches[i][:, 0], #Conjunto da componente X da i-esima branch 
                branches[i][:, 1], #Conjunto da componente Y da i-esima branch
                branches[i][:, 2], #Conjunto da componente Z da i-esima branch
                color = colors[i][0], #Cor da i-iesima branch
                linestyle = '-')   #Tipo de linha (no caso sera a linha continua)
                #zorder = 4)        #Ordem de prioridade

    # if(alpha != 0):
    #     lv.plotVecsLambdaZWB(ax, alpha, branches, 30) #Plotagem dos vetores
    #Plot do ponto fixo
    x_bar, y_bar = fun.map(Point[0], Point[1], baricentrica())
    Point_3d = [x_bar, y_bar, Point[2]]

    ax.scatter(*Point_3d, color='black', s=40, depthshade=False)

def integralCurveProjections(Point, branches, colors):
    from includes.Inicia import baricentrica

    fig1, ax_UZ = plt.subplots(figsize = (5,5)) 
    fig2, ax_VZ = plt.subplots(figsize = (5,5))
    fig3, ax_UV = plt.subplots(figsize = (5,5))

    #Mapeamento do ponto inicial (Point)
    u0, v0, z0 = Point
    u0, v0 = fun.map(u0, v0, baricentrica())
    
    for i in range(len(branches)):
        #Projecao U x Z
        ax_UZ.plot(branches[i][:, 0],
                   branches[i][:, 2],
                   color = colors[i][0],
                   linestyle = '-')
        ax_UZ.scatter(u0, v0, color = 'black', s = 10)
        
        #Projecao V x Z 
        ax_VZ.plot(branches[i][:, 1],
                   branches[i][:, 2],
                   color = colors[i][0],
                   linestyle = '-')
        ax_VZ.scatter(v0, z0, color = 'black', s = 10)

        #Projecao U x V
        ax_UV.plot(branches[i][:, 0],
                   branches[i][:, 1],
                   color = colors[i][0],
                   linestyle = '-')
        ax_UV.scatter(u0, v0, color = 'black', s = 10)

    ax_UZ.set_xlabel("u")
    ax_UZ.set_ylabel("z")
    ax_UZ.set_title("Projection u vs z")
    ax_UZ.set_xlim(0,1)
    ax_UZ.set_ylim(0,1)
    ax_UZ.grid(True)

    ax_VZ.set_xlabel("v")
    ax_VZ.set_ylabel("z")
    ax_VZ.set_title("Projection v vs z")
    ax_VZ.set_xlim(0,1)
    ax_VZ.set_ylim(0,1)
    ax_VZ.grid(True)

    ax_UV.set_xlabel("u")
    ax_UV.set_ylabel("v")
    ax_UV.set_title("Projection u vs v")
    triangulo(ax_UV, baricentrica())
    ax_UV.set_xlim(0,1)
    ax_UV.set_ylim(0,1)
    ax_UV.grid(True)

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

    fig1, ax_UZ = plt.subplots(figsize = (8,5)) 
    fig2, ax_VZ = plt.subplots(figsize = (8,5))

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
        integralCurveProjections(Point, branches, colors)

    plt.show()