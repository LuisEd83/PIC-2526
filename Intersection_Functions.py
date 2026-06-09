import includes.Campo_Hugoniot as ch

import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import fsolve

##############################################
#            Variaveis globais               #
##############################################
#Define-se as variaveis                      #
alpha = 0.1                                  #
u0, v0, z0 = 0.2, 0.6, 0.3                   #
z = 0                                        #
num_N = 10                                   #
dz = 1/(num_N * 100)                         #
#Variavel para resolucao                     #
factor = 2                                   #
Resol = 500 * factor                         #
#Variavel de tolerancia                      #
TOL_residual = 1e-8                          #
TOL = 1e-3                                   #
##############################################

#Sistema contendo as funcoes acima
def hugoniotSystem(vars):
    u_s, v_s = vars #Extrai as variaveis do sistema
    return [
        ch.F(u0, v0, z0, u_s, v_s, z),         #Funcoes definidas acima
        ch.G(alpha, u0, v0, z0, u_s, v_s, z)   #Funcoes definidas acima        
    ]


#Grade de valores (respeitando u + v <= 1)
u_vals = np.linspace(0.01, 0.98, Resol)   #Cria a reta u com 200 pontos
v_vals = np.linspace(0.01, 0.98, Resol)   #Cria a reta v com 200 pontos
U, V   = np.meshgrid(u_vals, v_vals)      #Mescla u e v e retorna uma tupla de valores

#Mascara regioes invalidas (u + v > 1), i.e, armazena a regiao que nao eh de interesse nosso
mask = (U + V > 1)

ZF = ch.F(u0, v0, z0, U, V, z).astype(float)           #Cria uma 'lista' de pontos que foram varridos dados U e V acima
ZG = ch.G(alpha, u0, v0, z0, U, V, z).astype(float)    #Cria uma 'lista' de pontos que foram varridos dados U e V acima


#filtra da grade os pontos (u, v) onde simultaneamente |F| ≈ 0 e |G| ≈ 0
"""
[Explicacao] -> Esses pontos jah estao perto da solucao real,
                o que aumenta a chance de convergência do fsolve
                (que serah utilizada posteriormente).
"""
proximity = (np.abs(ZF) < TOL) & (np.abs(ZG) < TOL) & ~mask
u_guesses = U[proximity]    #Variaveis de chute (guesses)
v_guesses = V[proximity]    #Variaveis de chute (guesses)

#Lista para armazenar as possiveis solucoes
candidatas = []

for u_g, v_g in zip(u_guesses, v_guesses):
    try:
        """
        [Explicacao] -> o fsolve eh projetado para encontrar raizes
                        de equacoes ou sistemas de equacoes nao li-
                        neares.
                    ->  u_g e v_g sao variaveis de chute (guess) jah
                        filtradas anteriormente.
        """
        sol = fsolve(hugoniotSystem, [u_g, v_g], full_output=True)
        u_s, v_s = sol[0]                                                            #Possivel solucao do sistema

        residual = np.linalg.norm(hugoniotSystem([u_s, v_s]))                        #Calcula a norma dos valores obtidos do sistema
        inside  = (u_s > 0) and (v_s > 0) and (u_s + v_s < 1)                        #Verifica se u_s e v_s estao dentro do dominio desejado
        not_dup = all(
            np.linalg.norm([u_s - us, v_s - vs]) > TOL for us, vs, zs in candidatas  #Verifica se u_s e v_s nao eh duplicata de uma solucao jah encontrada
            )
        
        if ((residual < TOL_residual) and (inside) and (not_dup)):
            candidatas.append((u_s, v_s, z))

    except:
        continue

print("--------------------------------------------------")
print(f"Número de soluções encontradas: {len(candidatas)}")
print("--------------------------------------------------")
print("########################################################################")
for i in range(len(candidatas)):
    #Print das solucoes aproximadas:
    print(f"Solucao candidata {i+1}: {candidatas[i]}")
print("########################################################################")

################_______________________________PLOTAGEM_______________________________################
#Limpa os pontos invalidos das listas acima criadas
ZF[mask] = np.nan
ZG[mask] = np.nan

#Criacao da figura
fig, ax = plt.subplots(figsize=(7, 6))

#######################################################################################################
#Reta u + v = 1
u_value = [1, 0]
v_value = [0, 1]
plt.plot(u_value, v_value, linestyle = '-', color = 'k') 
#######################################################################################################

#Curva F = 0 e Curva G = 0
ax.contour(U, V, ZF, levels=[z], colors='blue',  linewidths=2)  #O level (curva de nivel) eh determinada pelo valor de z
ax.contour(U, V, ZG, levels=[z], colors='red', linewidths=2)    #O level (curva de nivel) eh determinada pelo valor de z

#Legenda manual
ax.plot([], [], color = 'blue', label=f'F(u, v, {z}) = 0')      #Legenda -> Curva da funcao F
ax.plot([], [], color = 'red',  label=f'G(u, v, {z}) = 0')      #Legenda -> Curva da funcao G
ax.plot([], [], color = 'black', label = 'u + v = 1')           #Legenda -> Curva u + v = 1

#Configuracoes da figure
ax.set_title(f'Level curve z = {z}')    #Titulo da figura
ax.set_xlabel('u')                      #Legenda do eixo X
ax.set_ylabel('v')                      #Legenda do eixo Y
ax.legend()                             #Ativa a legenda
ax.set_xlim(0, 1)                       #Limite do "plano carteziano" no eixo X     
ax.set_ylim(0, 1)                       #Limite do "plano carteziano" no eixo Y

plt.grid(True)                          #Habilita a grade
plt.tight_layout()                      #Ajusta automaticamente o espaçamento
plt.show()                              #Mostra a figura
#######################################################################################################

"""
#######################################################################################################
#Codigo simbolico (para caso seja utilizado posteriormente)
#######################################################################################################
import sympy as sp
def a_sympy(z):
    return sp.atan(z)

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

print(G_sym)

#######################################################################################################
"""
