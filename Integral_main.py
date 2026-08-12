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
import includes.Campo_Hugoniot_Teste as cht
from includes.Campo_Hugoniot import F, G


import numpy as np
import matplotlib.pyplot as plt

#Definindo constantes:
h = 0.001        #Passo do metodo de Euler
N = 1500        #Numero de realizacao de passo
alpha = 0.0    #Variável de controle
Num_z = 20      #Numero de curvas de nivel

#Definindo a funcao principal de plot
def main():
    #Importando biblioteca para nao haver erro durante a lida do arquivo Functions.py
    import includes.Inicia as ini
    global h

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
    u0, v0, z0 = iniPoints[0][0], iniPoints[0][1], iniPoints[0][2]

    #Itera sobre cada ponto em iniPoints
    for Point in iniPoints:
        #Variavel que armazena as branches
        branches = b.Branches_point(alpha, Point, [h, N])
        colors = b.Branches_point_colors(alpha, branches)

        #__________INICIALIZANDO PLOTAGEM DA CURVA INTEGRAL__________#
        ic.plotIntegralCurve(ax, Point, alpha, branches, colors)

        #__________INICIALIZANDO PLOTAGEM DAS PROJECOES__________#
        #ic.integralCurveProjections(Point, branches, colors)

    arrayPos = cht.HugonioutEuler_method(alpha, *iniPoints, [h, N], tol = 0.01)
    arrayNeg = cht.HugonioutEuler_method(alpha, *iniPoints, [-h, N], tol = 0.01)

    # 1. Extrai apenas a parte [u, v, z] de cada elemento da li sta de pontos
    PointsPos = np.array([item[0] for item in arrayPos])
    PointsNeg = np.array([item[0] for item in arrayNeg])

   # arrayNeg = arrayNeg[1:]
    arrayPos = arrayPos[1:]

    sigP = np.array([item[1] for item in arrayPos])
    # sigN = np.array([item[1] for item in arrayNeg])
    sigP = sigP[1:]
    # sigN = sigN[1:]

    # 3. Agora a concatenação vai funcionar perfeitamente, pois ambos são matrizes Nx3 puras
    # points = np.concatenate((PointsPos, PointsNeg))

    # print("====================================================")
    # print("HUGONIOT")
    # print(f"Tamanho do array de pontos: {len(points)}")
    # for i in range(len(points) - 1):
    #     ax.scatter(points[i][0], points[i][1], points[i][2], color='black', marker='.', zorder = 10)
    #     if(i%10 == 0):
    #         print(f"Valor de F e G no ponto: F = {F(u0, v0, z0, points[i][0], points[i][1], points[i][2])} | G = {G(alpha, u0, v0, z0, points[i][0], points[i][1], points[i][2])}")
    # print("====================================================")
    

    #ax.plot(PointsNeg[:,0], PointsNeg[:,1], PointsNeg[:,2], linestyle = '-', color = 'magenta')
    ax.plot(PointsPos[:,0], PointsPos[:,1], PointsPos[:,2], linestyle = '-', color = 'orange')

    # ax.scatter(PointsPos[-1][0], PointsPos[-1][1], PointsPos[-1][2], color = 'red', marker = '.')

    # print("==========================================")
    # ax.scatter(PointsNeg[-1][0], PointsNeg[-1][1], PointsNeg[-1][2], color = 'blue', marker = '.')
    # u, v, z = PointsNeg[-1]
    # u0, v0, z0 = iniPoints[0][0], iniPoints[0][1], iniPoints[0][2] 
    # print(f"Ultimo ponto da integração errada: ({u, v, z})")
    # print(f"Diferenças absolutas: ({abs(u - u0), abs(v - v0), abs(z - z0)})")
    # print(f"valor do campo no ponto máximo: {cht.campo(alpha, u, v, z, sigN[-1], u0, v0, z0)}")
    # print(f"norma do vetor no ponto máximo: {cht.norm_sig(alpha, u, v, z, sigN[-1], u0, v0, z0)}")

    # print("==========================================")
    # idxmaxPoint = np.argmax(PointsNeg[:, 2])
    # maxPoint = PointsNeg[idxmaxPoint]
    # u, v, z = maxPoint

    # ax.scatter(maxPoint[0], maxPoint[1], maxPoint[2], color = 'darkgreen', marker = '.', s = 10)
    # print(f"Idx do ponto máximo: {idxmaxPoint}")
    # print(f"Ponto maximo: ({u, v, z})")
    # print(f"Valor de F e G no ponto: F = {F(u0, v0, z0, u, v, z)} | G = {G(alpha, u0, v0, z0, u, v, z)}")
    # print(f"Diferenças absolutas: ({abs(u - u0), abs(v - v0), abs(z - z0)})")
    # print(f"valor do campo no ponto máximo: {cht.campo(alpha, u, v, z, sigN[idxmaxPoint], u0, v0, z0)}")
    # print(f"norma do vetor no ponto máximo: {cht.norm_sig(alpha, u, v, z, sigN[idxmaxPoint], u0, v0, z0)}")
    # print(f"lambdaS = {fun.lbdas(u, v, z)}")
    # print(f"lambdaF = {fun.lbdaf(u, v, z)}")
    # print(f"lambdaZ = {af.lambdz(u, v, z, alpha)}")
    # print(f"Sigma = {sigN[idxmaxPoint]}")


    # print("==========================================")
    idxmaxPoint = np.argmax(PointsPos[:, 2])
    maxPoint = PointsPos[idxmaxPoint]
    u, v, z = maxPoint

    ax.scatter(maxPoint[0], maxPoint[1], maxPoint[2], color = 'yellow', marker = 'o', s = 10)
    print(f"Idx do ponto máximo: {idxmaxPoint}")
    print(f"Ponto maximo: ({u, v, z})")
    print(f"Valor de F e G no ponto: F = {F(u0, v0, z0, u, v, z)} | G = {G(alpha, u0, v0, z0, u, v, z)}")
    # print(f"Diferenças absolutas: ({abs(u - u0), abs(v - v0), abs(z - z0)})")
    print(f"valor do campo no ponto máximo: {cht.campo(alpha, u, v, z, sigP[- 1], u0, v0, z0)}")
    print(f"norma do vetor no ponto máximo: {cht.norm_sig(alpha, u, v, z, sigP[- 1], u0, v0, z0)}")
    print(f"LAMBDAS:")
    print(f"lambdaS = {fun.lbdas(u, v, z)}")
    print(f"lambdaF = {fun.lbdaf(u, v, z)}")
    print(f"lambdaZ = {af.lambdz(u, v, z, alpha)}")
    print(f"Sigma = {sigP[- 1]}")
    
    #__________INICIALIZANDO PLOTAGEM DE GRAFICO DAS VELOCIDADES__________#
    #ag.lamb_graph(alpha, [h, N])
    
    #__________INICIALIZANDO PLOTAGEM DA SUPERFICIE DE COINCIDENCIA__________#
    # base = cs.argmin_branch(branches)
    # auto_branches = b.Branches_auto([h, N], Num_z, alpha, base - 0.5, altura = 1)
    # auto_colors = b.Branches_auto_colors()

    #cs.plotCoicindenceCurves(ax, auto_branches, auto_colors)

    #Inicia plotagem
    plt.show()

main()