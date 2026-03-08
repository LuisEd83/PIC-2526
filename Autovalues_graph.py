"""
Modulo: Grafico de autovalores


Objetivo:
- Plotar um único grafico;
- No grafico haverao 2 eixos: x-axis = tempo, y-axis = valores dos autovalores.
- Ainda no grafico haverao 3 curvas

- Extrair os pontos da curva dentro do prisma para encontrar todos os valores dos autovalores.

"""

from Auxiliar_Functions import Array_Concatenated, Points_Filter, lambdz
from Inicia import baricentrica

import matplotlib.pyplot as plot

def lamb_graph(alpha, Point : list, integ_config : list):
    #Importando biblioteca para nao haver erro durante a lida do arquivo Functions.py
    import Functions as fun
    
    #mapeia para baricentrica ou não
    mp = baricentrica()
    
    #Y-AXIS#
    #Colecao de pontos dentro do prisma:
    array_points_c = Array_Concatenated(alpha, Point, integ_config) #Extrai todos os pontos
    
    #Inicializando vetores dos lambdas
    array_LambS = []
    array_LambF = []
    array_LambZ = []

    #Armazenando todos os valores dos lambdas:
    for i in range(len(array_points_c)):
        array_LambS.append(fun.lbdas(*array_points_c[i]))
        array_LambF.append(fun.lbdaf(*array_points_c[i]))
        array_LambZ.append(lambdz(*array_points_c[i], alpha))
    
    #X-AXIS#
    #Inicializando lista
    k = []

    for i in range(2*integ_config[1] + 1):
        k.append(i/(2*integ_config[1] + 1))

    #Iniciando a figura
    plot.figure()

    #Plotando pontos
    plot.plot(k, array_LambS, color = 'b', linestyle = '-', label = '\u03BBs')
    plot.plot(k, array_LambF, color = 'r', linestyle = '-', label = '\u03BBf')
    plot.plot(k, array_LambZ, color = 'purple', linestyle = '-', label = '\u03BBz')

    plot.grid(True)
    plot.xlabel("Eixo X - Passos")
    plot.ylabel("Eixo Y - Valores numéricos")

    plot.show()

lamb_graph(0.3, [0.5, 0.15, 0.1], [0.01, 500])