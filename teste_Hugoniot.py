import includes.Campo_Hugoniot_Teste as cht
import includes.Inicia as ini
import includes.Auxiliar_Functions as af
import includes.Functions as fun


import numpy as np
import matplotlib.pyplot as plt

h, N = 0.01, 500
alpha = 0.0
fixedPoint = [0.5, 0.3, 0.2]

arrayPos = cht.HugonioutEuler_method(alpha, fixedPoint, [h, N])
arrayNeg = cht.HugonioutEuler_method(alpha, fixedPoint, [-h, N])

arrayNeg = arrayNeg[1:]
arrayPos = arrayPos[1:]

# 1. Extrai apenas a parte [u, v, z] de cada elemento da lista de pontos
PointsPos = np.array([item[0] for item in arrayPos])
PointsNeg = np.array([item[0] for item in arrayNeg])
print(PointsPos[-1])
print(PointsNeg[-1])


# 2. Se precisar dos sigmas isolados em algum lugar:
# sigP = np.array([item[1] for item in arrayPos])
# sigN = np.array([item[1] for item in arrayNeg])

# 3. Agora a concatenação vai funcionar perfeitamente, pois ambos são matrizes Nx3 puras
points = np.concatenate((PointsPos, PointsNeg))
print(points.size)

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


ax.plot(points[:, 0], points[:, 1], points[:, 2], color='black', linestyle='-')
ax.scatter(fixedPoint[0], fixedPoint[1], fixedPoint[2], color = 'red', marker = 'o', s = 3)

plt.show()

