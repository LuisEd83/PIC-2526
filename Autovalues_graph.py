"""
Modulo: Grafico de autovalores


Objetivo:
- Plotar um único grafico;
- No grafico haverao 2 eixos: x-axis = tempo, y-axis = valores dos autovalores.
- Ainda no grafico haverao 3 curvas

- Extrair os pontos da curva dentro do prisma para encontrar todos os valores dos autovalores.

"""
import Numericals_methods as nm
from Auxiliar_Functions import lambdz

import numpy as np
import matplotlib.pyplot as plot

def lamb_graph(alpha, Point : list, integ_config : list):
    #Importando biblioteca para nao haver erro durante a lida do arquivo Functions.py
    import Functions as fun
    
    #Y-AXIS#
    #Colecao de pontos dentro do prisma:
    array_ph = nm.Euler_method(alpha, Point, integ_config) #Array com h > 0
    array_mh = nm.Euler_method(alpha, Point, [-integ_config[0], integ_config[1]]) #Array com h < 0

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
            #Inverte o array se o primeiro elemento do array_mh estiver mais proximo do Point
            if(np.linalg.norm(array_mh[-1] - Point) > np.linalg.norm(array_mh[0] - Point)): 
                array_mh = np.flip(array_mh, axis = 0) #Inverte apenas os elementos do array

            org_points = np.concatenate([array_mh, np.array(Point).reshape(1, -1), array_ph])
            break
        elif(coluna_ph[0] - coluna_mh[i] < 0):
            #Inverte o array se o primeiro elemento do array_ph estiver mais proximo do Point
            if(np.linalg.norm(array_ph[-1] - Point) > np.linalg.norm(array_ph[0] - Point)):
                array_ph = np.flip(array_ph, axis = 0) #Inverte apenas os elementos do array

            org_points = np.concatenate([array_ph, np.array(Point).reshape(1, -1), array_mh])
            break
        else:
            i += 1

        if(i == len(coluna_mh)):
            print("[ERROR] - Impossibilidade de determinar ordem")
            exit()
    
    #Procuranto o Point
    index_Point = 0         #Inicia em zero
    for i in range(len(org_points)):
        if((org_points[i][0] == Point[0]) and (org_points[i][1] == Point[1]) and (org_points[i][2] == Point[2])):
            index_Point = i                                   #Armazena o index do Point no org_points
            index_Point = index_Point/(2*integ_config[1] + 1) #Corrige para a escala do grafico

    #Inicializando vetores dos lambdas
    array_LambS = []
    array_LambF = []
    array_LambZ = []

    #Armazenando todos os valores dos lambdas:
    for i in range(len(org_points)):
        array_LambS.append(fun.lbdas(*org_points[i]))
        array_LambF.append(fun.lbdaf(*org_points[i]))
        array_LambZ.append(lambdz(*org_points[i], alpha))
    
    #X-AXIS#
    #Inicializando lista
    k = []

    for i in range(2*integ_config[1] + 1):
        k.append(i/(2*integ_config[1] + 1))

    #Iniciando a figura
    plot.figure()

    plot.axvline(x = index_Point, color = 'k', linestyle = '--')

    #Plotando pontos
    plot.plot(k, array_LambS, color = 'b', linestyle = '-', label = '\u03BBs')
    plot.plot(k, array_LambF, color = 'r', linestyle = '-', label = '\u03BBf')
    plot.plot(k, array_LambZ, color = 'purple', linestyle = '-', label = '\u03BBz')

    plot.grid(True)
    plot.xlabel("Eixo X - Passos")
    plot.ylabel("Eixo Y - Valores numéricos")

    plot.show()

lamb_graph(0.04, [0.5, 0.45, 0.1], [0.01, 500])