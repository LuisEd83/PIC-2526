#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 11:17:14 2026

@author: cido
"""

"""
Modulo: Integracao pelo metodo de Euler

Objetivo:
- Utilizar o metodo de Euler para armazenar o conjunto de pontos que informam o comportamento
dos compoentes do campo a partir de um ponto inicial escolhido.
"""

import includes.Campo_Ez as df
import includes.Campo_Hugoniot as ch

from includes.Functions import lbdc
from includes.Inicia import baricentrica

from numpy import array, ceil, abs, log10, sqrt, linspace

#Definindo uma funcao que implementa o metodo de Euler para integracao
def Euler_method(alpha, fixed_point : list, P_chute : list, integ_config : list):
    #Extraindo ponto fixado
    u0, v0, z0 = fixed_point 

    #Extraindo ponto inicial para integracao
    u1, v1, z1 = P_chute

    #Extraindo configuracao de integracao
    h, N = integ_config
    
    #Inicializando variaveis para o metodo de Euler
    uk = u1
    vk = v1
    zk = z1

    #Inicializando uma lista de pontos:
    Points = [[u1, v1, z1]]

    #Calculando e armazenando os resultados do metodo de Euler:
    for _ in range(N): #O laco vai repetir N vezes
        #Realizo o passo (obs: kp1 = k + 1)
        #(alpha, u0, v0, z0, u, v, z):
        ukp1 = uk + h * ch.Hug1(alpha, u0, v0, z0, uk, vk, zk)
        vkp1 = vk + h * ch.Hug2(alpha, u0, v0, z0, uk, vk, zk)
        zkp1 = zk + h * ch.Hug3(alpha, u0, v0, z0, uk, vk, zk)

        #Armazenando os valores na lista de pontos
        Points.append([ukp1, vkp1, zkp1])

        #Atualizo os valores das variaives uk, vk e zk 
        uk = ukp1
        vk = vkp1
        zk = zkp1

    #Definindo um array para armazenar pontos (com tamanho N + 1):
    Points = array(Points, float)
    return Points
