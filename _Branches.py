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
from Auxiliar_Functions import transparencia

import numpy as np

"""
1) Acertar a plotagem;
2) Alpha pequeno (máx 0.1)
3) Curva de Nível (lambdaS com lambdaZ e lambdaZ com lambdaF)
 -> lambdaS com lambdaZ = azul
 -> lambdaZ com lambdaF = vermelho
"""


#Redefinindo funcoes para nao haver erro circular
def lambdz(u, v, z, alpha):
    return(fw(u, v, z)/(u + alpha*np.cos(z)))

def Branches(alpha, Point : list, integ_config : list):
    #Determinando se deverao passar todos os pontos ou nao
    t = not transparencia()

    #Constante:
    sqr3 = np.sqrt(3)
    def pointInT(Point):
        #Extraindo pontos
        u, v, z = Point

        if(baricentrica):
            #Teorema de Viviani
            h1 = np.abs(v)
            h2 = (np.abs(-sqr3 * u + v))/2
            h3 = (np.abs(sqr3 * u + v - sqr3))/2

            h = sqr3/2

            if((np.isclose(h1 + h2 + h3, h)) and (0.0 <= z <= 1.0)):
                return 1
        else:
            #Definicao de triangulo retangulo
            if((0 <= u <= 1) and (0 <= v <= 1) and (0 <= u + v <= 1) and (0 <= z <= 1)):
                return 1
            
        return t #0 se transparencia == 1 e 1 se transparencia == 0
    
    def indexs(array_org, lista_index : list, bool_value : bool, current_size, correc_error):
        #Condicao de parada
        if(current_size <= 0):
            return

        #Criando array para armazenar parte da lista
        array_temp = np.array([])

        for i in range(len(array_org)):
            if(pointInT(array_org[i]) == int(not bool_value)):
                bool_value = not bool_value                                                      #Inverte o valor do bool_value
                lista_index.append(correc_error + i)                                             #Armazena o index que ocorreu o if e soma com o erro (advindo da recursividade)
                correc_error += (i+1)                                                            #Adiciona o erro atual com o antigo
                current_size -= (i+1)                                                            #Retiro parte do tamanho do array
                array_temp = array_org[i+1:].copy()                                              #Cria uma copia do array a partir de um indice i
                return indexs(array_temp, lista_index, bool_value, current_size, correc_error)


    """
    
    #ORGANIZAÇÃO DOS PONTOS
    
    """

    array_ph = ei.Euler_method(alpha, Point, integ_config)                          #Array com h > 0
    array_mh = ei.Euler_method(alpha, Point, [-integ_config[0], integ_config[1]])   #Array com h < 0

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
                array_ph = np.flip(array_ph, axis = 0) 

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
    
    if(transparencia()): #Se o usuario quiser com os pontos internos ao prisma
        #Criando a variaveis chamar funcao:
        points_indexs = []              #Armazenamento de indices
        bool_value = False              #Variavel booleana para controle
        size_array = len(org_points)    #Variavel para armazenar tamanho do array dos pontos organizados
        
        if(pointInT(org_points[0])):    #Determina valor booleano para funcao
            bool_value = not bool_value
        
        #Retorna os indices onde ocorre a transicao de "dentro do prisma" para "fora do prisma" e vice-versa
        indexs(org_points, points_indexs, bool_value, size_array, correc_error = 0)

        #Inicializando variaveis para loop:
        p = 0                                   #Variavel de posicao
        inside = bool(pointInT(org_points[0]))  #Variavel de controle

        #Isso alterna a escolha dos pontos
        for i in range(len(points_indexs)):
            if(inside): #Se os pontos estiverem dentro do prisma, ele os seleciona
                branches_list.append(org_points[p : points_indexs[i]])

            #Note que, se inside == False, ele apenas pula para o proximo indice, nao escolhendo os pontos que estao fora
            p = points_indexs[i]
            inside = not inside

        if(inside):
            branches_list.append(org_points[p:])    #Adicionando o ultimo segmento se a ultima iteracao fizer Inside = True

        #for i in range(len(branches_list)):
        #    print(f"Branch {i+1}: {branches_list[i]} \n")

        if(pointInT(Point)):
            #Localizando a branch e o indice do Point
            bp = 0       #Armazena a branch do Point
            index_bp = 0 #Armazena o indice do Point na branch 
            for i, branch in enumerate(branches_list):
                for j, p in enumerate(branch):
                     if np.all(p == Point):
                        bp = i
                        index_bp = j
            
            if(index_bp != 0): #Verifica se o Point é o primeiro elemento da branch do Point
                #Dividindo a primeira parte da branch do Point:
                branch_atual = branches_list[bp]

                branch_1 = branch_atual[:index_bp + 1] #Fatia até o Point (o incluindo)
                branch_2 = branch_atual[index_bp:]     #Fatia a partir do Point (o incluindo)

                branches_list.pop(bp)                  #Remove a branch do Point
                branches_list.insert(bp, branch_1)     #Substitui a branch que estava na posicao bp pela fatia que "vai" até o Point
                branches_list.insert(bp + 1, branch_2) #Cria uma nova branch que inicia a partir do Point

            #Agora deve-se tratar a branch bp + 1:
            branch_atual = branches_list[bp+1]
            index_max = np.argmax(branches_list[bp+1][:, 2]) #Extrai o indice do ponto maximo da branch bp + 1

            Point_max = branch_atual[index_max]
            if(not np.isclose(np.linalg.norm(Point_max - branch_atual[-1]), 0)): #Verifica se o ponto maximo nao eh o ultimo elemento da branch
                branch_1 = branch_atual[: index_max + 1]          #Fatia a branch a partir do Point ate o o ponto maximo da branch (o incluindo)
                branch_2 = branch_atual[index_max:]               #Fatia a branch a partir do ponto maximo (o incluinto)

                branches_list.pop(bp+1)                           #Remove a branch do ponto maximo
                branches_list.insert(bp+1, branch_1)              #Substitui a branch que estava na posicao bp+1 
                branches_list.insert(bp+2, branch_2)              #Cria uma nova branch que inicia a partir do ponto maximo

    else:
        #Transforma em um array numpy e recebe os pontos organizados
        branches_list = np.array([org_points])

    for i in range(len(branches_list)):
        print(f"Branch {i+1}: {branches_list[i]} \n")

    return branches_list #retorna como lista, para melhor eficiencia

#A funcao a seguir retorna um array em que cada elemento terah um valor do tipo
# 'b' - lambdaS => Azul
# 'r' - lambdaF => Vermelho
# 'purple' - LambdaZ => Roxo
# 'k' - ERROR
def Branches_colors(alpha, branches):
    #Definindo lista que irah armazenar os numeros ('b', 'r' ou 'purple')
    colors = []

    for i in range(len(branches)):
        #Extrai o ponto cenntral (da metade do ramo)
        point = branches[i][len(branches[i])//2] #Divisão de inteiro

        #Definindo valores e comparando:
        lambdaS_value = lbdas(*point)
        lambdaF_value = lbdaf(*point)
        lambdaZ_value = lambdz(*point, alpha)

        if((lambdaZ_value < lambdaF_value < lambdaS_value) or (lambdaZ_value < lambdaS_value < lambdaF_value)):
            colors.append('b') #Adiciona a cor Azul
        elif((lambdaS_value < lambdaZ_value < lambdaF_value) or (lambdaF_value < lambdaZ_value < lambdaS_value)):
            colors.append('r') #Adiciona a cor Vermelha
        elif((lambdaF_value < lambdaS_value < lambdaZ_value) or (lambdaS_value < lambdaF_value < lambdaZ_value)):
            colors.append('purple') #Adiciona a cor Roxa
        else:
            colors.append('k')

    return colors