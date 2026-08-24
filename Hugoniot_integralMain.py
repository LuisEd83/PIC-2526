import matplotlib.pyplot as plt
import includes.Auxiliar_Functions as af
import includes.Functions as fun
import includes._Branches as b

import Autovalues_graph as ag
import Hugoniot_IntegralCurve as hic
import Coincidence_surfaces as cs

#Variavel do grafico                         #
iValue = 0.01                                #
fValue = 0.99                                #
eMask = True                                 #
factor = 2                                   #
Resol = 500 * factor                         #
Num_z = 20

#Para curva
alpha = 0.1

#######################################################
#Para as funções feitas à mão
h, N = 0.0005 , 3500
integ_config = [h, N]
#######################################################

def mainHugoniot():
    import includes.Inicia as ini

    #mapeia para baricentrica ou não
    mp = ini.baricentrica()

    #Define as concentracoes minima e maxima
    zmin, zmax = 0.0, 1.0

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
        fixed_point = Point

        #Variavel que armazena as branches
        branchesImplicitas, colors, branches = b.Branches_Hugoniot(alpha, fixed_point, integ_config, iValue, fValue, Resol, eMask)

        #__________INICIALIZANDO PLOTAGEM DO GRAFICO DOS AUTOVALORES E DO SIGMA__________#
        ag.lamb_graph(alpha, fixed_point, [h, N], branches, True)

        # #__________INICIALIZANDO PLOTAGEM DA CURVA HUGONIOT E DA CURVA IMPLICITA__________#
        hic.plotHugoniotCurve(ax, alpha, fixed_point, branches, branchesImplicitas, colors)

        #__________INICIALIZANDO PLOTAGEM DAS PROJECOES CURVA HUGONIOT__________#
        #hic.plotHugoniotProjecoes(fixed_point, branches, branchesImplicitas, colors)


    #__________INICIALIZANDO PLOTAGEM DA SUPERFICIE DE COINCIDENCIA__________#
    # auto_colors, auto_branches = b.Branches_auto([h, N], Num_z, alpha, 0.05, altura = 1.0)

    # cs.plotCoicindenceCurves(ax, auto_branches, auto_colors)

    plt.show()




mainHugoniot()