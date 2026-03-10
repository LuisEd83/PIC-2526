"""
Modulo: Branches

Objetivo:
- Criar uma função que retorna ramos (e sua respectiva coloracao)

Os ramos devem:
- Ter uma coloracao de acordo com o autovalor lambdaS, lambdaF e lambdaZ  

"""

import Euler_Integration as ei
from Inicia import baricentrica
from Functions import fw, lbdas, lbdaf 

import numpy as np

#Redefinindo funcoes para nao haver erro circular
def az(z):
    return np.cos(z)

def lambdz(u, v, z, alpha):
    return(fw(u, v, z)/(u + alpha*az(z)))

def Branches(alpha, Point : list, integ_config : list):
    #Definindo escolha (baricentrica):
    mp = baricentrica

    #Constante:
    sqr3 = np.sqrt(3)
    def pointInT(Point):
        #Extraindo pontos
        u, v, z = Point

        if(mp):
            #Teorema de Viviani
            h1 = np.abs(v)
            h2 = (np.abs(-sqr3 * u + v))/2
            h3 = (np.abs(sqr3 * u + v - sqr3))/2

            h = sqr3/2

            if((h1 + h2 + h3 - h == 0) and (0.0 <= z <= 1.0)):
                return 1
        else:
            #Definicao de triangulo retangulo
            if((0.0 <= u <= 1) and (0.0 <= v <= 1) and (0.0 <= u + v <= 1) and (0.0 <= z <= 1.0)):
                return 1
            
        return 0
    
    def indexs(array_org, lista_index : list, bool_value : bool, current_size, correc_error):
        #Condicao de parada
        if(current_size <= 0):
            return -1

        #Criando array para armazenar parte da lista
        array_temp = np.array([])

        if(bool_value):
            for i in range(len(array_org)): #Para percorrer todo o array de pontos.
                if(pointInT(*array_org[i]) == 0):
                    bool_value = False                                                               #Inverte o valor do bool_value
                    lista_index.append(correc_error + i)                                             #Armazena o index que ocorreu o if e soma com o erro (advindo da recursividade)
                    correc_error += (i+1)                                                            #Adiciona o erro atual com o antigo
                    current_size -= (i+1)                                                            #Retiro parte do tamanho do array
                    array_temp = array_org[i+1:].copy()                                              #Cria uma copia do array a partir de um indice i
                    return indexs(array_temp, lista_index, bool_value, current_size, correc_error)   #Chama a funcao novamente
        else:
            for i in range(len(array_org)):
                if(pointInT(*array_org[i]) == 1):
                    bool_value = True
                    lista_index.append(correc_error + i)
                    correc_error += (i+1)
                    current_size -= (i+1)
                    array_temp = array_org[i+1:].copy()
                    return indexs(array_temp, lista_index, bool_value, current_size, correc_error)


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

    #Criando a variaveis chamar funcao:
    points_indexs = []              #Armazenamento de indices
    bool_value = False              #Variavel booleana para controle
    size_array = len(org_points)    #Variavel para armazenar tamanho do array dos pontos organizados
    
    if(pointInT(*org_points[0])): #Determina valor booleano para funcao
        bool_value = True
    
    #Retorna os indices onde ocorre a transicao de "dentro do prisma" para "fora do prisma" e vice-versa
    indexs(org_points, points_indexs, bool_value, size_array, correc_error = 0)

    #Criando a variavel para armazenar as branches:
    branches_list = []

    #variavel para armazenar o inicio da posicao:
    p = 0
    for i in range(len(points_indexs)):
        branches_list.append([org_points[p : points_indexs[i]]])
        p = points_indexs[i]

    branches_list.append(org_points[p:])    #Adicionando o ultimo segmento
    branches = np.array(branches_list)

    return branches

#A funcao a seguir retorna um array em que cada elemento terah um valor do tipo
# 1 - lambdaS => Azul
# 2 - lambdaF => Vermelho
# 3 - LambdaZ => Roxo
# -1 - ERROR
def Branches_colors(alpha, branches):
    #Definindo lista que irah armazenar os numeros (1, 2 ou 3)
    colors = []

    for i in range(len(branches)):
        #Extrai o ponto cenntral (da metade do ramo)
        point = branches[i][len(branches[i])//2] #Divisão de inteiro

        #Definindo valores e comparando:
        lambdaS_value = lbdas(*point)
        lambdaF_value = lbdaf(*point)
        lambdaZ_value = lambdz(*point, alpha)

        if((lambdaZ_value < lambdaF_value < lambdaS_value) or (lambdaZ_value < lambdaS_value < lambdaF_value)):
            colors.append(1) #Adiciona a cor Azul
        elif((lambdaS_value < lambdaZ_value < lambdaF_value) or (lambdaF_value < lambdaZ_value < lambdaS_value)):
            colors.append(2) #Adiciona a cor Vermelha
        elif((lambdaF_value < lambdaS_value < lambdaZ_value) or (lambdaS_value < lambdaF_value < lambdaZ_value)):
            colors.append(3) #Adiciona a cor Roxa
        else:
            colors.append(-1)

    return np.array(colors)
