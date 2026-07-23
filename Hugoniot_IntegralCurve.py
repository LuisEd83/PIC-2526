import includes.Numericals_methods as nm
import matplotlib.pyplot as plt
import includes.Inicia as ini
import includes.Auxiliar_Functions as af
import includes.Functions as fun
import includes.Campo_Hugoniot as ch
import includes._Branches as b
import numpy as np

#Variavel do grafico                         #
iValue = 0.01                                #
fValue = 0.99                                #
eMask = True                                 #
factor = 2                                   #
Resol = 500 * factor                         #

#Para curva
alpha = 0.0
fixed_point = [0.425, 0.5235, 0.2 ]
#0.1 0.3 0.2

#######################################################
#Para as funções feitas à mão
h, N = 0.001, 1000
integ_config = [h, N]
#######################################################

#mapeia para baricentrica ou não
mp = ini.baricentrica()

#Define as concentracoes minima e maxima
zmin, zmax = ini.concentrations() 

ax = ini.ambiente3d()
ax.view_init(elev=30., azim=-130.) #Initial Camera Position
fig, (ax_uz, ax_vz) = plt.subplots(1, 2, figsize=(12, 5)) #Para as projecoes

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


branches = b.Branches_Hugoniot(iValue, fValue, Resol, eMask, alpha, fixed_point, integ_config)

#Plot do ponto fixo no ambiente 3D
ax.plot(fixed_point[0], fixed_point[1], fixed_point[2], marker='o', color='blue')

#Plot do Ponto fixo nas projeções
ax_uz.plot(fixed_point[0], fixed_point[2], marker='o', color='blue')
ax_vz.plot(fixed_point[1], fixed_point[2], marker='o', color='blue')

print("**********************************************************")
print("Valor das funções F e G no ponto fixo:")
print(f"Valor de F(fixed_P): {ch.F(fixed_point[0], fixed_point[1], fixed_point[2], fixed_point[0], fixed_point[1], fixed_point[2])};")
print(f"Valor de G(fixed_P): {ch.G(alpha, fixed_point[0], fixed_point[1], fixed_point[2], fixed_point[0], fixed_point[1], fixed_point[2])}")
print("**********************************************************")

print(f"Tamanho das branches = {len(branches)}")

for i, branch in enumerate(branches):
    if len(branch) == 0:
        continue

    branch = np.array(branch)

    print(f"[Branch {i}] {len(branch)} pontos")

    ax.plot(
        branch[:, 0],
        branch[:, 1],
        branch[:, 2],
        marker='o',
        linestyle='-',
        markersize=2
    )

    #Projecao u-z (coluna 0 vs coluna 2)
    ax_uz.plot(
        branch[:, 0],
        branch[:, 2],
        marker='o',
        linestyle='-',
        markersize=2
    )

    #Projecao v-z (coluna 1 vs coluna 2)
    ax_vz.plot(
        branch[:, 1],
        branch[:, 2],
        marker='o',
        linestyle='-',
        markersize=2
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

plt.tight_layout()
plt.show()