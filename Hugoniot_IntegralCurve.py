import matplotlib.pyplot as plt
import includes.Auxiliar_Functions as af
import includes.Functions as fun
import includes.Inicia as ini
import LambdaZ_vec as lv

import matplotlib.pyplot as plt
from numpy import array

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


def plotHugoniotCurve(
        ax,
        alpha,
        fixed_point,
        branches            : list,
        branchesImplicita,
        colors
):
    from includes.Inicia import baricentrica

    #Mapeia para as coordenadas baricentricas se necessario
    if(baricentrica()):
        for branch in branches:
            for p in branch:
                u, v, z = p[0]

                u, v = fun.map(u, v, baricentrica())

                p[0][0] = u
                p[0][1] = v

        for branch in branchesImplicita:
            for pointBranchImplicita in branch:
                u, v = pointBranchImplicita

                u, v = fun.map(u, v, baricentrica())

                pointBranchImplicita[0] = u
                pointBranchImplicita[1] = v

    """for i in range(1, len(branches)):
        ponto_anterior = branches[i-1][-1]
        branches[i].insert(0, ponto_anterior)"""

    #Plotagem
    for i, branch in enumerate(branches):
        points = array([p[0] for p in branch])
        ax.plot(points[:, 0],
                points[:, 1],
                points[:, 2],
                color = colors[i],
                linestyle = '-')

    for branch in branchesImplicita:
        points = array(branch)
        ax.plot(points[:, 0],
                points[:, 1],
                fixed_point[2],
                color = 'k',
                linestyle = '-')

    #Abre espaco para os cones (vetores 3D)
    # if(alpha != 0):
        # lv.plotVecsLambdaZWB(ax, alpha, branches, 25, True) #Plotagem dos vetores

    #Plot do ponto fixo
    x_bar, y_bar = fun.map(fixed_point[0], fixed_point[1], ini.baricentrica())
    fixed_point_3d = [x_bar, y_bar, fixed_point[2]]

    ax.scatter(*fixed_point_3d, color='black', s=40, depthshade=False)


def plotHugoniotProjecoes(
        fixed_point,
        branches,
        branchesImplicitas,
        colors
):
    #Plot do ponto fixo no ambiente 3D
    x_bar, y_bar = fun.map(fixed_point[0], fixed_point[1], ini.baricentrica())
    fixed_point_3d = [x_bar, y_bar, fixed_point[2]]

    u0, v0, z0 = fixed_point_3d

    fig1, ax_uz = plt.subplots(figsize = (5,5)) 
    fig2, ax_vz = plt.subplots(figsize = (5,5))
    fig3, ax_uv = plt.subplots(figsize = (5,5))

    #Plot do Ponto fixo nas projecoes
    ax_uz.plot(u0, z0, marker='o', color='k')
    ax_vz.plot(v0, z0, marker='o', color='k')
    ax_uv.plot(u0, v0, marker='o', color='k')

    for i, branch in enumerate(branches):
        if len(branch) == 0:
            continue

        u_vals = [point[0][0] for point in branch]
        v_vals = [point[0][1] for point in branch]
        z_vals = [point[0][2] for point in branch]

        #============
        #Projecao u-z
        #============
        ax_uz.plot(
            u_vals,
            z_vals,
            linestyle='-',
            color = colors[i],
            markersize=2
        )

        #============
        #Projecao v-z
        #============
        ax_vz.plot(
            v_vals,
            z_vals,
            linestyle='-',
            color = colors[i],
            markersize=2
        )

        
        #============
        #Projecao u-v
        #============
        #Pontos da branch Hugoniot
        ax_uv.plot(
            u_vals,
            v_vals,
            linestyle = '-',
            color = colors[i],
            markersize = 2
        )

    #==========================
    #Plot da Hugoniot implicita
    #==========================
    for branch in branchesImplicitas:
        #Pontos da branch Hugoniot implicita
        u_imp = [p[0] for p in branch]
        v_imp = [p[1] for p in branch]

        ax_uv.plot(
            u_imp,
            v_imp,
            linestyle = '--',
            color = 'k',
            markersize = 2
        )

    #Projecao u vs z
    ax_uz.set_xlabel("u")
    ax_uz.set_ylabel("z")
    ax_uz.set_title("Projection u vs z")
    ax_uz.set_xlim(0, 1)
    ax_uz.set_ylim(0, 1)
    ax_uz.grid(True)

    #Projecao v vs z
    ax_vz.set_xlabel("v")
    ax_vz.set_ylabel("z")
    ax_vz.set_title("Projection v vs z")
    ax_vz.set_xlim(0, 1)
    ax_vz.set_ylim(0, 1)
    ax_vz.grid(True)

    #Projecao u vs v
    ax_uv.set_xlabel("u")
    ax_uv.set_ylabel("v")
    ax_uv.set_title("Projection u vs v")
    triangulo(ax_uv, ini.baricentrica())
    ax_uv.set_xlim(0, 1)
    ax_uv.set_ylim(0, 1)
    ax_uv.grid(True)


    plt.tight_layout()