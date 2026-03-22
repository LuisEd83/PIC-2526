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
import Autovalues_graph as ag

import numpy as np
import matplotlib.pyplot as plt

#Definindo constantes:
h = 0.01        #Passo do metodo de Euler
N = 500         #Numero de realizacao de passo
alpha = 0.04    #Variável de controle
Num_z = 30      #Numero de curvas de nivel

sqr3 = np.sqrt(3) 

#Definindo o ponto inicial para teste
U0 = [0.5, 0.3, 0.3]

U1 = [0.5, 0.1, 0.1] #Ponto importante

U2 = [0.2, 0.4, 0.1] #Ponto inicial fora do prisma

#Função responsavel por determinar o menor Z (altura)
def argmin_branch(branches):
    Zm = 0 #Menor Z (teorico)
    for i in range(len(branches)):
        for j in range(len(branches[i])):
            if(branches[i][j][2] < Zm):
                Zm = branches[i][j][2] #Troca o Zm para o novo ponto com altura menor
    return Zm

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
                linestyle = '-',   #Tipo de linha (no caso sera a linha continua)
                zorder = 4)        #Ordem de prioridade

    #__________INICIALIZANDO PLOTAGEM DOS PONTOS RELACIONADA AO AUTOVALOR LAMBDA-Z__________#
    #Plotando o ponto inicial
    ax.plot(*Point, 'k.', zorder = 5)

    marcados = b.andar_branches(Point, branches, passo = 10) #Criar uma lista de tuplas onde estas são (i,j) -> i = branch e j = posicao na branch
    cor_pontos = []                                          #Armazena a cor do ponto
    for i,j in marcados:
        af.color_point(cor_pontos, branches[i][j], alpha)
    
    print(cor_pontos)

    for l, (i, j) in enumerate(marcados):
        ax.plot(*branches[i][j], color=cor_pontos[l], marker='.', zorder=5) #Plota o ponto marcado


    #__________INICIALIZANDO PLOTAGEM DA CURVA RELACIONADA AOS AUTOVALORES__________#
    base = argmin_branch(branches)
    auto_branches = b.Branches_auto([h, N], Num_z, alpha, base, altura = 1)
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
    elif((not af.branchSlow()) and af.branchFast()): #Plotagem do ramo lambda_f == lambda_z
        for i in range(len(auto_branches)):
            seg = np.array(auto_branches[i][0])
            ax.plot(seg[:, 0],
                    seg[:, 1],
                    seg[:, 2],
                    color = auto_colors[0],
                    linestyle = '-')
    elif(af.branchSlow() and (not af.branchFast())): #Plotagem do ramo lambda_s == lambda_z
        for i in range(len(auto_branches)):
            seg = np.array(auto_branches[i][1])
            ax.plot(seg[:, 0],
                    seg[:, 1],
                    seg[:, 2],
                    color = auto_colors[1],
                    linestyle = '-')


    ag.lamb_graph(alpha, Point, [h, N])

    #Inicia plotagem
    plt.show()

Prism_plot(U2)