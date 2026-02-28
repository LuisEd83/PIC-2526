"""
Modulo: Integracao pelo metodo de Euler

Objetivo:
- Utilizar o metodo de Euler para armazenar o conjunto de pontos que informam o comportamento
dos campos a partir de um ponto escolhido.
"""

import Determinants_Functions as df
from numpy import array 

#Definindo uma funcao que implementa o metodo de Euler para integracao
def Euler_method(point : list, integ_config : list):
    #Extraindo ponto inicial
    u0, v0, z0 = point

    #Extraindo configuracao de integracao
    h, N = integ_config

    #Definindo arrays (listas) que vao armazenar os resultados e inicializando com os valores iniciais
    Array_u = [u0]
    Array_v = [v0]
    Array_z = [z0]
    
    #Inicializando variaveis para o metodo de Euler
    uk = u0
    vk = v0
    zk = z0

    #Calculando e armazenando os resultados do metodo de Euler:
    for _ in range(N): #O laco vai repetir N vezes
        #Realizo o passo (obs: kp1 = k + 1)
        ukp1 = uk + h * df.P(uk, vk, zk)
        vkp1 = vk + h * df.Q(uk, vk, zk)
        zkp1 = zk + h * df.R(uk, vk, zk)

        #Armazeno os valores em seus respectivos arrays
        Array_u.append(ukp1)
        Array_v.append(vkp1)
        Array_z.append(zkp1)

        #Atualizo os valores das variaives uk, vk e zk de acorddo com o passo h 
        uk += h
        vk += h
        zk += h

    #Definindo um array para armazenar pontos (com tamanho N + 1):
    Array = array([Array_u, Array_v, Array_z])

    return Array