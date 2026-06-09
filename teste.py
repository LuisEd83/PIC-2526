import includes.Numericals_methods as nm
import matplotlib.pyplot as plt
import includes.Inicia as ini
import includes.Auxiliar_Functions as af
import includes.Functions as fun
import includes.Campo_Hugoniot as ch

import Euler_hug as eh

def pointInP(Point):
    if(not af.transparencia()):
        return 1
    if(ini.baricentrica()):
        return (af.if_PointInEq(Point) and (0.0 <= Point[2] <= 1.0))
    else:
        return (af.if_PointInRet(Point) and (0.0 <= Point[2] <= 1.0))

###fixed_point = [0.55, 0.3, 0.2]###
fixed_point = [0.2, 0.6, 0.3]

alpha = 0.0
#Teste 1
#Point = [0.5503759, 0.284621, 0.0]
#Point = [0.3253, 0.5597, 0.0]
#Point = [0.5413, 0.2924, 0.0]

#Teste 2
#Point = [0.1764, 0.6192, 0.0]
#Point = [0.169, 0.6265, 0.0]
Point = [0.1999999999999947, 0.5606192233373988, 0.0]

#######################################################
#Para o scipy
integ_config1 = { #h > 0
    's_inicial': 0.0,
    's_final': 20.0,
    'n_pontos': 20
}

integ_config2 = { #h < 0
    's_inicial': 20.0,
    's_final': 0.0,
    'n_pontos': 20
}
#######################################################

#######################################################
#Para as funções feitas à mão
h, N = 0.1, 200
integ_config = [h, N]
#######################################################

#mapeia para baricentrica ou não
mp = ini.baricentrica()

#Define as concentracoes minima e maxima
zmin, zmax = ini.concentrations() 

ax = ini.ambiente3d()
ax.view_init(elev=30., azim=-130.) #Initial Camera Position

if(af.transparencia() or 1):
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

## Runge-Kutta feito à mão
#P = nm.RungeKutta4(alpha, fixed_point, Point, integ_config)
#P1 = nm.RungeKutta4(alpha, fixed_point, Point, [-integ_config[0], integ_config[1]])

## Euler feito à mão
#P = eh.Euler_method(alpha, fixed_point, Point, integ_config)
P = eh.Euler_method(alpha, fixed_point, Point, [-integ_config[0], integ_config[1]])

## Scipy
#P = nm.runge_Kutta_Scipy(alpha, fixed_point, Point, integ_config1)
#P1 = nm.runge_Kutta_Scipy(alpha, fixed_point, Point, integ_config2)

#PLot do ponto inicial
ax.plot(Point[0], Point[1], Point[2], marker = 'o', color = 'red')
ax.plot(fixed_point[0], fixed_point[1], fixed_point[2], marker = 'o', color = 'blue')

Points_filtered = []
Points_filtered1 = []

print("**********************************************************")
print("Valor das funções F e G no ponto fixo:")
print(f"Valor de F(fixed_P): {ch.F(fixed_point[0], fixed_point[1], fixed_point[2], fixed_point[0], fixed_point[1], fixed_point[2])};")
print(f"Valor de G(fixed_P): {ch.G(alpha, fixed_point[0], fixed_point[1], fixed_point[2], fixed_point[0], fixed_point[1], fixed_point[2])}")
print("**********************************************************")

for i in range(len(P)):
    if(pointInP(P[i]) or (not af.transparencia())):
        Points_filtered.append(P[i])
        print(f"Valor de F(P[{i}]): {ch.F(fixed_point[0], fixed_point[1], fixed_point[2], P[i][0], P[i][1], P[i][2])};")
        print(f"Valor de G(P[{i}]): {ch.G(alpha, fixed_point[0], fixed_point[1], fixed_point[2], P[i][0], P[i][1], P[i][2])}")
        print("--------------------------------------------")


for i in range(len(Points_filtered)):
    ax.plot(
        P[i][0],
        P[i][1],
        P[i][2],
        marker = '.',
        color = 'k',
        linestyle = '-'
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