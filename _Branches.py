"""
Modulo: Branches

Objetivo:
- Criar uma função que retorna ramos (e sua respectiva coloracao)

Os ramos devem:
- Ter uma coloracao de acordo com o autovalor lambdaS, lambdaF e lambdaZ  

"""

import Euler_Integration as ei
from Inicia import baricentrica

import numpy as np

def Branches(alpha, Point : list, integ_config : list):
    #Extraindo os pontos:
    u, v, z = Point

    #Extraindo configuracao:
    h, N = integ_config

    #Definindo escolha (baricentrica):
    mp = baricentrica

    #Constante:
    sqr3 = np.sqrt(3)
    def pointInT(Point, bar):
        #Extraindo pontos
        u, v, = Point

        if(bar):
            #Teorema de Viviani
            h1 = np.abs(v)
            h2 = (np.abs(-sqr3 * u + v))/2
            h3 = (np.abs(sqr3 * u + v - sqr3))/2

            h = sqr3/2

            if(h1 + h2 + h3 - h == 0):
                return 1
        else:
            #Definicao de triangulo retangulo
            if((0 <= u <= 1) and (0 <= v <= 1) and (0 <= u + v <= 1)):
                return 1
            
        return 0
    
    """
    
    #ORGANIZAÇÃO DOS PONTOS
    
    """

    array_ph = ei.Euler_method(alpha, Point, integ_config) #Array com h > 0
    array_mh = ei.Euler_method(alpha, Point, [-integ_config[0], integ_config[1]]) #Array com h < 0

    #Retirando o ponto inicial (Point):
    array_ph = array_ph[~np.all(array_ph == Point, axis = 1)]
    array_mh = array_mh[~np.all(array_mh == Point, axis = 1)]

    #Extraindo as colunas destes arrays:
    coluna_ph = array_ph[:, 2]
    coluna_mh = array_mh[:, 2]

    #Criando uma variável para armazenar os pontos na ordem correta
    org_points = np.array([])

    #Variavel de controle de laco
    i = 0

    while(1):
        if(coluna_ph[0] - coluna_mh[i] > 0):
            org_points = np.concatenate([array_mh, np.array(Point).reshape(1, -1), array_ph])
            break
        elif(coluna_ph[0] - coluna_mh[i] < 0):
            org_points = np.concatenate([array_ph, np.array(Point).reshape(1, -1), array_mh])
            break
        else:
            i += 1

        if(i == len(coluna_mh)):
            print("[ERROR] - Impossibilidade de determinar ordem")
            exit()

    """
    
    #CRIAÇÃO DAS BRANCHES
    
    """

    #Criando a variavel para armazenar as branches:
    branches_list = []

    #Primeiro branch -> armazena os pontos do org_points que vem antes do Point
    index_point = np.where(np.all(org_points == Point, axis=1))[0][0]

    first_branch = org_points[:index_point + 1] #Copia parte do org_points

    #Selecionando o index que contem maior valor de z
    #indexHZ = points_concatened_C[:, 2].argmax()
