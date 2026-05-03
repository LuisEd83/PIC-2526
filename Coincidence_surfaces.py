import includes.Functions as fun
import includes._Branches as b
import includes.Auxiliar_Functions as af

import numpy as np
import matplotlib.pyplot as plt

#Função responsavel por determinar o menor Z (altura)
def argmin_branch(branches):
    Zm = 0 #Menor Z (teorico)
    for i in range(len(branches)):
        for j in range(len(branches[i])):
            if(branches[i][j][2] < Zm):
                Zm = branches[i][j][2] #Troca o Zm para o novo ponto com altura menor
    return Zm

def levelCurveWB(ax, auto_branches, auto_colors):
    if(af.branchSlow() and af.branchFast()):        #Plotagem dos dois ramos
        for i in range(len(auto_branches)):
            for j in range(len(auto_branches[i])):
                seg = np.array(auto_branches[i][j])
                ax.plot(seg[:, 0],
                        seg[:, 1],
                        seg[:, 2],
                        color = auto_colors[j],
                        linestyle = '-')
    elif(af.branchSlow() and (not af.branchFast())): #Plotagem do ramo lambda_f == lambda_z
        for i in range(len(auto_branches)):
            seg = np.array(auto_branches[i][0])
            ax.plot(seg[:, 0],
                    seg[:, 1],
                    seg[:, 2],
                    color = auto_colors[0],
                    linestyle = '-')
    elif((not af.branchSlow()) and af.branchFast()): #Plotagem do ramo lambda_s == lambda_z
        for i in range(len(auto_branches)):
            seg = np.array(auto_branches[i][1])
            ax.plot(seg[:, 0],
                    seg[:, 1],
                    seg[:, 2],
                    color = auto_colors[1],
                    linestyle = '-')

def levelCurve(alpha, curveConfig, numCortes):
    #Importando biblioteca para nao haver erro durante a lida do arquivo Functions.py
    import includes.Inicia as ini

    h, N = curveConfig

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
    
    #__________INICIALIZANDO PLOTAGEM DA CURVA RELACIONADA AOS AUTOVALORES__________#
    #branches = b.Branches_point(alpha, Point, [h, N])
    #base = argmin_branch(branches)
    auto_branches = b.Branches_auto([h, N], numCortes, alpha, base = 0, altura = 1)
    auto_colors = b.Branches_auto_colors()

    if(af.branchSlow() and af.branchFast()):        #Plotagem dos dois ramos
        for i in range(len(auto_branches)):
            for j in range(len(auto_branches[i])):
                seg = np.array(auto_branches[i][j])
                ax.plot(seg[:, 0],
                        seg[:, 1],
                        seg[:, 2],
                        color = auto_colors[j],
                        linestyle = '-')
    elif(af.branchSlow() and (not af.branchFast())): #Plotagem do ramo lambda_f == lambda_z
        for i in range(len(auto_branches)):
            seg = np.array(auto_branches[i][0])
            ax.plot(seg[:, 0],
                    seg[:, 1],
                    seg[:, 2],
                    color = auto_colors[0],
                    linestyle = '-')
    elif((not af.branchSlow()) and af.branchFast()): #Plotagem do ramo lambda_s == lambda_z
        for i in range(len(auto_branches)):
            seg = np.array(auto_branches[i][1])
            ax.plot(seg[:, 0],
                    seg[:, 1],
                    seg[:, 2],
                    color = auto_colors[1],
                    linestyle = '-')

    #Inicia plotagem
    plt.show()

levelCurve(0.1, [0.01, 500], 10)