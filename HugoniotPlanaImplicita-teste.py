# -*- coding: utf-8 -*-
"""
Spyder Editor

Curva de Hugoniot plana pelo ponto base (u0, v0, z0).

Funciona bem no triangulo sem mapeamento
"""
# 
 
import numpy as np  #funcoes matemáticas
#import sympy as sym #Funcoes simbolicas
import matplotlib.pyplot as plt #Para plotagem de graficos
#from mpl_toolkits.mplot3d import axes3d

import includes.Functions as fun
import includes.Inicia as ini
#import DadosIniciais as di

umin, umax, vmin, vmax, triang = ini.dominio()


#mapeia para baricentrica ou não
mp = 0

# Fixa as viscosidades
#u0, v0, z0  = di.Ul()
muw0, muo, mug = 1.0, 2.0, 0.5


# Concentracoes iniciais
cmin, cmax = 0.0, 1.0


#The fixed point U0
"""u0 = 0.439
v0 = 0.5274
z0 = 0.685"""

"""
Possui dois pontos singulares
u0 = 0.6
v0 = 0.3
z0 = 0.5
"""
u0 = 0.418
v0 = 0.555
z0 = 0.0


# Define um grid no quadrado [0, 1] \times [0, 1]
x = np.linspace(umin, umax, 500)
y = np.linspace(vmin, vmax, 500)
X, Y = np.meshgrid(x, y)
#P, Q = fun.map(X,Y, mp) #fun.mapeia o grid para o triangulo equilatero


#Plotagem de Hugoniots planas

#Hugoniot por L no plano z =  z0, como curva de nivel zero
def H(u,v):
    return (fun.fw(u,v,z0) - fun.fw(u0,v0,z0))*(v-v0) - (fun.fo(u,v,z0) - fun.fo(u0,v0,z0))*(u-u0)

plt.figure('saturation triangle')

#Curva de Hugoniot no nivel c fornecido
Hl = H(X, Y) 

mask = (X + Y < 1)
Hl = np.where(X + Y < 1, Hl, np.nan)
# Guarde a referência da curva na variável 'cs'
cs = plt.contour(X, Y, Hl, [0], linestyles="dashed", colors='k')

pontos_hugoniot = []

# Forma atualizada (Matplotlib 3.8+): acesse get_paths() direto do cs
for path in cs.get_paths():
    v = path.vertices
    u_pts = v[:, 0]  # Coordenadas X (u)
    v_pts = v[:, 1]  # Coordenadas Y (v)
    
    # Organiza em matriz de pares (u, v)
    pts = np.column_stack((u_pts, v_pts))
    pontos_hugoniot.append(pts)

# Se houver ao menos um segmento de curva encontrado:
if pontos_hugoniot:
    pontos = pontos_hugoniot[0]
    # print("Primeiros 5 pontos da curva (u, v):")
    # print(pontos[:5]))
    for i in range(2000):
        plt.scatter(pontos[i][0], pontos[i][1], color = 'red', marker = '.')
    # plt.plot(pontos[:, 0], pontos[:, 1], linestyle = '-', color = 'blue')
    # Para salvar em arquivo:
    np.savetxt('pontos_hugoniot.txt', pontos, header='u v', comments='')

u_value = [1, 0]
v_value = [0, 1]
plt.plot(u_value, v_value, linestyle = '-', color = 'red') 

# Plota o ponto base 
# Esta plotagem deve vir após a plotagem implicita do contour
u0, v0 = fun.map(u0,v0,mp)
plt.plot(u0, v0, 'ko')

plt.grid(True)

plt.show()


