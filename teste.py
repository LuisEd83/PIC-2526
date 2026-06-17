import includes.Numericals_methods as nm
import matplotlib.pyplot as plt
import includes.Inicia as ini
import includes.Auxiliar_Functions as af
import includes.Functions as fun
import includes.Campo_Hugoniot as ch
import includes._Branches as b

import Euler_hug as eh

#Variavel do grafico                         #
iValue = 0.01                                #
fValue = 0.99                                #
eMask = True                                 #
factor = 2                                   #
Resol = 500 * factor                         #

#Para curva
alpha = 0.1
fixed_point = [0.1, 0.5, 0.2]
dz = 0.09

def pointInP(Point):
    if(not af.transparencia()):
        return 1
    if(ini.baricentrica()):
        return (af.if_PointInEq(Point) and (0.0 <= Point[2] <= 1.0))
    else:
        return (af.if_PointInRet(Point) and (0.0 <= Point[2] <= 1.0))

#_______Fixed_Points______#
#fixed_point = [0.2, 0.6, 0.3]  #Fixed_point importante

#_______Point_______#
#Point = [0.1999999999999947, 0.5606192233373988, 0.0] #ESTE PONTO EH ESPECIAL

#######################################################
#Para as funções feitas à mão
h, N = 0.05, 200
integ_config = [h, N]
#######################################################

#mapeia para baricentrica ou não
mp = ini.baricentrica()

#Define as concentracoes minima e maxima
zmin, zmax = ini.concentrations() 

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


P = b.Branches_Hugoniot(iValue, fValue, Resol, eMask, alpha, fixed_point, integ_config, dz)

#PLot do ponto inicial
#ax.plot(Point[0], Point[1], Point[2], marker = 'o', color = 'red')
ax.plot(fixed_point[0], fixed_point[1], fixed_point[2], marker = 'o', color = 'blue')

Points_filtered = []
#Points_filtered1 = []

print("**********************************************************")
print("Valor das funções F e G no ponto fixo:")
print(f"Valor de F(fixed_P): {ch.F(fixed_point[0], fixed_point[1], fixed_point[2], fixed_point[0], fixed_point[1], fixed_point[2])};")
print(f"Valor de G(fixed_P): {ch.G(alpha, fixed_point[0], fixed_point[1], fixed_point[2], fixed_point[0], fixed_point[1], fixed_point[2])}")
print("**********************************************************")

for i in range(len(P)):
    for j in range(len(P[i])):

        point = P[i][j]

        if pointInP(point) or (not af.transparencia()):

            Points_filtered.append(point)

            print(f"Ponto [{i}][{j}]: {point}")

            print(
                f"Valor de F: "
                f"{ch.F(fixed_point[0], fixed_point[1], fixed_point[2],
                        point[0], point[1], point[2])}"
            )

            print(
                f"Valor de G: "
                f"{ch.G(alpha,
                        fixed_point[0], fixed_point[1], fixed_point[2],
                        point[0], point[1], point[2])}"
            )

            print("--------------------------------------------")


for point in Points_filtered:
    ax.plot(
        point[0],
        point[1],
        point[2],
        marker='.',
        color='k',
        linestyle='None'
    )

"""
for i in range(len(P1)):
    if(pointInP(P1[i]) or (not af.transparencia())):
        Points_filtered1.append(P1[i])

for i in range(len(Points_filtered1)):
    ax.plot(
        P1[i][0],
        P1[i][1],
        P1[i][2],
        marker = 'o',
        color = 'k',
        linestyle = '-'
    )
        
"""

plt.show()



"""
#######################################################################################
fixed_point = [0.4, 0.2, 0.3]
Point = [0.400000000000119, 0.154386641840665, 0.0] 
Solucao candidata 1: (np.float64(0.400000000000119), np.float64(0.154386641840665), 0)
Solucao candidata 2: (np.float64(0.3765256472230938), np.float64(0.19279246676909106), 0)
Solucao candidata 3: (np.float64(0.24279522599290781), np.float64(0.43675449656584436), 0)
#######################################################################################
"""

"""
## Runge-Kutta feito à mão
#P = nm.RungeKutta4(alpha, fixed_point, Point, integ_config)
#P1 = nm.RungeKutta4(alpha, fixed_point, Point, [-integ_config[0], integ_config[1]])

## Euler feito à mão
#P = eh.Euler_method(alpha, fixed_point, Point, integ_config)
#P1 = eh.Euler_method(alpha, fixed_point, Point, [-integ_config[0], integ_config[1]])

## Scipy
#P = nm.runge_Kutta_Scipy(alpha, fixed_point, Point, integ_config1)
#P1 = nm.runge_Kutta_Scipy(alpha, fixed_point, Point, integ_config2)
"""