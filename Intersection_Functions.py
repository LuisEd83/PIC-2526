import includes.Numericals_methods as nm
import includes.Inicia as ini
import includes.Auxiliar_Functions as af
import includes.Functions as fun

import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

def a_sympy(z):
    return sp.atan(z)

def F(u0, v0, z0, u, v, z):
    #Define-se os termos (para simplicidade):
    firstT  = (fun.fw(u, v, z) - fun.fw(u0, v0, z0)) * (v - v0)
    secondT = (fun.fo(u, v, z) - fun.fo(u0, v0, z0)) * (u - u0)

    return firstT - secondT

def G(alpha, u0, v0, z0, u, v, z):
    #Define-se os termos (para simplicidade):
    firstT  = (fun.fw(u, v, z) - fun.fw(u0, v0, z0)) * (u0 * (z - z0) + alpha * (af.a(z) - af.a(z0)))
    secondT = fun.fw(u0, v0, z0) * (z - z0) * (u - u0)

    return firstT - secondT

"""
#######################################################################################################
Codigo simbolico (para caso seja utilizado posteriormente)
#######################################################################################################
#Define-se variaveis simbolicas:
alpha = sp.symbols("alpha")         #Variavel de controle
u, v, z    = sp.symbols("u v z")    #Variaveis
u0, v0, z0 = sp.symbols("u0 v0 z0") #Ponto fixo

#Define-se as funcoes simbolicas 
F_sym = sp.Lambda((u0, v0, z0, u, v, z), F(u0, v0, z0, u, v, z))                #Funcao F simbolica
G_sym = sp.Lambda((alpha, u0, v0, z0, u, v, z), G(alpha, u0, v0, z0, u, v, z))  #Funcao G simbolica

#Expande e simplifique as funcoes
#F_sym = sp.expand(F_sym)
F_sym = sp.simplify(F_sym)
G_sym = sp.simplify(G_sym)
#######################################################################################################
"""

#Define-se as variaveis
alpha = 0.1
u0, v0, z0 = 0.55, 0.3, 0.2
z = 0

#Variavel para resolucao
Resol = 500

#Grade de valores (respeitando u + v <= 1)
u_vals = np.linspace(0.01, 0.98, Resol)   #Cria a reta u com 200 pontos
v_vals = np.linspace(0.01, 0.98, Resol)   #Cria a reta v com 200 pontos
U, V   = np.meshgrid(u_vals, v_vals)    #Mescla u e v e retorna uma tupla de valores

#Mascara regioes invalidas (u + v > 1), i.e, armazena a regiao que nao eh de interesse nosso
mask = (U + V >= 1)

ZF = F(u0, v0, z0, U, V, z).astype(float)           #Cria uma 'lista' de pontos que foram varridos dados U e V acima
ZG = G(alpha, u0, v0, z0, U, V, z).astype(float)    #Cria uma 'lista' de pontos que foram varridos dados U e V acima

#Limpa os pontos invalidos das listas acima criadas
ZF[mask] = np.nan
ZG[mask] = np.nan

#Plot
fig, ax = plt.subplots(figsize=(7, 6))

#Curva F = 0
ax.contour(U, V, ZF, levels=[z], colors='blue',  linewidths=2) #O level (curva de nivel) eh determinada pelo valor de z

#Curva G = 0
ax.contour(U, V, ZG, levels=[z], colors='red', linewidths=2) #O level (curva de nivel) eh determinada pelo valor de z

#Legenda manual
ax.plot([], [], color = 'blue', label=f'F(u, v, {z}) = 0')
ax.plot([], [], color = 'red',  label=f'G(u, v, {z}) = 0')
ax.plot([], [], color = 'black', label = 'u + v = 1')

#Reta u + v = 1
u_value = [1, 0]
v_value = [0, 1]
plt.plot(u_value, v_value, linestyle = '-', color = 'k') 

ax.set_title(f'Level curve z = {z}')
ax.set_xlabel('u')
ax.set_ylabel('v')
ax.legend()
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)

plt.grid(True)
plt.tight_layout()
plt.show()
