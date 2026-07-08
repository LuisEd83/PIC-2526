import includes.Auxiliar_Functions as af
import includes.Campo_Hugoniot as ch

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

##############################################
#            Variaveis globais               #
##############################################
#Define-se as variaveis                      #
alpha = 0.5                                  #
u0, v0, z0 = 0.5, 0.3, 0.2                   #
z = 1.0                                     #
num_N = 10                                   #
dz = 1/(num_N * 100)                         #
#Variavel do grafico                         #
iValue = 0.01                                #
fValue = 0.99                                #
eMask = True                                 #
factor = 2                                   #
Resol = 500 * factor                         #
#Variavel de tolerancia                      #
TOL_residual = 1e-8                          #
TOL = 1e-3                                   #
##############################################

matplotlib.use('TkAgg')  #backend mais fluido que o padrão

U, V, ZF, ZG, c = af.hugoniotSystemSolver(iValue, fValue, Resol, eMask, u0, v0, z0, z, alpha) #Descarto os candidatos usando _

print("--------------------------")
print("Valor campo Hugoniot:")
for i in range(len(c)):
    print(f"Ponto [{c[i][0], c[i][1], c[i][2]}]: ({ch.CampoHug(alpha, u0, v0, z0, c[i][0], c[i][1], c[i][2])})")
print("--------------------------")

print("--------------------------")
print("Valor do lambdaZ no ponto de interseccao:")
for i in range(len(c)):
    print(f"Ponto {i+1}: {af.lambdz(*c[i], alpha)}")
print("Valor do lambdaZ no ponto fixo:")
print(f"Ponto (fixo): {af.lambdz(u0, v0, z0, alpha)}")
print("--------------------------")

#Mascara regioes invalidas (u + v > 1), i.e, armazena a regiao que nao eh de interesse nosso
mask = (U + V > 1)

################_______________________________PLOTAGEM_______________________________################
#Limpa os pontos invalidos das listas acima criadas
if(eMask):
    ZF[mask] = np.nan
    ZG[mask] = np.nan

#Criacao da figura
fig = plt.figure(figsize=(14,6))

#___________SUBFIGURA 1___________#
#Variavel da primeira sub figura
ax1 = fig.add_subplot(121)

#######################################################################################################
#Reta u + v = 1
u_value = [1, 0]
v_value = [0, 1]
plt.plot(u_value, v_value, linestyle = '-', color = 'k') 
#######################################################################################################

#Curva F = 0 e Curva G = 0
ax1.contour(U, V, ZF, levels=[0], colors='blue', linewidths=2)   #O level (curva de nivel) eh determinada pelo valor de z
ax1.contour(U, V, ZG, levels=[0], colors='red', linewidths=2)    #O level (curva de nivel) eh determinada pelo valor de z

#Regiao onde estao localizadas as raizes do sistema
"""
[Explicacao] -> Esses pontos jah estao perto da solucao real,
                o que aumenta a chance de convergência do fsolve
                (que serah utilizada posteriormente).
"""
proximity = (np.abs(ZF) < TOL) & (np.abs(ZG) < TOL) & (~mask if eMask else True) #Calcula os indices onde ZF e ZG sao proximos de 0
u_guesses = U[proximity]    #Variaveis de chute (guesses)
v_guesses = V[proximity]    #Variaveis de chute (guesses)
ax1.scatter(u_guesses, v_guesses, color = "gray")

#Legenda manual
ax1.plot([], [], color = 'blue', label=f'F(u, v, {z}) = 0')      #Legenda -> Curva da funcao F
ax1.plot([], [], color = 'red',  label=f'G(u, v, {z}) = 0')      #Legenda -> Curva da funcao G
ax1.plot([], [], color = 'black', label = 'u + v = 1')           #Legenda -> Curva u + v = 1
ax1.plot([], [], color = 'gray', marker = 'o', linestyle = '', label = 'Root region')

#Configuracoes da primeira subfigure
ax1.set_title(f'Level curve z = {z}')    #Titulo da figura
ax1.set_xlabel('u')                      #Legenda do eixo X
ax1.set_ylabel('v')                      #Legenda do eixo Y
ax1.legend()                             #Ativa a legenda
ax1.set_xlim(0, 1)                       #Limite do "plano carteziano" no eixo X     
ax1.set_ylim(0, 1)                       #Limite do "plano carteziano" no eixo Y
ax1.grid(True)

#___________SUBFIGURA 2 (superficies)___________#
#Variavel da segunda sub figura
ax2 = fig.add_subplot(122, projection='3d')

#Eh o "pulo" entre os pontos que serao selecionados
step = 12  #Isto eh para deixar o plot mais fluido

#Escolha dos pontos para plot de superficie
Us  = U[::step, ::step]
Vs  = V[::step, ::step]
ZFs = ZF[::step, ::step]
ZGs = ZG[::step, ::step]

#Superficie com F(u,v,z) = 0 e G(u, v, z) = 0
ax2.plot_surface(Us, Vs, ZFs, cmap='Blues', alpha=0.4, #alpha diz respeito a opacidade
                 linewidth=0, antialiased=False)  #antialiased=False ajuda bastante no que diz respeito ao plot mais fluido

ax2.plot_surface(Us, Vs, ZGs, cmap='Reds', alpha=0.4,
                 linewidth=0, antialiased=False)

#Legenda manual
ax2.plot([], [], color = 'blue', label=f'F(u, v, {z})')      #Legenda -> Superficie da funcao F
ax2.plot([], [], color = 'red',  label=f'G(u, v, {z})')      #Legenda -> Superficie da funcao G

#Plot do plano values = 0
#Encontra os pontos de raiz (onde F = 0 e G = 0)
F_root = np.abs(ZFs) < TOL #Seleciona os indices onde ZFs eh menor que TOL
G_root = np.abs(ZGs) < TOL #Seleciona os indices onde ZGs eh menor que TOL

U_F_root, V_F_root, Z_F_root = Us[F_root], Vs[F_root], ZFs[F_root] #O ultimo eh altura F = 0
U_G_root, V_G_root, Z_G_root = Us[G_root], Vs[G_root], ZGs[G_root] #O ultimo eh altura F = 0

#Plot dos pontos da interseccao
ax2.scatter(U_F_root, V_F_root, Z_F_root, color='blue', marker = 'o', s=5, zorder = 4)
ax2.scatter(U_G_root, V_G_root, Z_G_root, color='red', marker = 'o', s=5, zorder = 4)

#Legenda manual
ax2.plot([],[], color = 'blue', marker = 'o', linestyle = '', label = 'F = 0')
ax2.plot([],[], color = 'red', marker = 'o', linestyle = '', label = 'G = 0')

ax2.set_title('F and G surfaces')
ax2.set_xlabel('u')
ax2.set_ylabel('v')
ax2.set_zlabel('value')
ax2.legend()

ax2.grid(True)                          #Habilita a grade
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
