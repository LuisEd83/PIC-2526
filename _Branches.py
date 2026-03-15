"""
Modulo: Branches

Objetivo:
- Criar uma função que retorna ramos (e sua respectiva coloracao)

Os ramos devem:
- Ter uma coloracao de acordo com o autovalor lambdaS, lambdaF e lambdaZ  

"""

import Numericals_methods as nm
from Inicia import baricentrica
from Functions import fw, lbdas, lbdaf 
from Auxiliar_Functions import transparencia, if_PointInEq, if_PointInRet

import numpy as np

""" 
1) Acertar a plotagem; (CONCLUIDO)
2) Alpha pequeno (máx 0.1)
3) Curva de Nível (lambdaS com lambdaZ e lambdaZ com lambdaF)
 -> lambdaS com lambdaZ = azul
 -> lambdaZ com lambdaF = vermelho
"""

#________________________BRANCHES DOS PONTOS________________________#

#Redefinindo funcoes para nao haver erro circular
def lambdz(u, v, z, alpha):
    return(fw(u, v, z)/(u + alpha*np.cos(z)))

def Branches_point(alpha, Point : list, integ_config : list):

    def pointInT(Point):
        if(not transparencia()):
            return 1

        if(baricentrica()):
            return (if_PointInEq(Point) and (0.0 <= Point[2] <= 1.0))
        else:
            return (if_PointInRet(Point) and (0.0 <= Point[2] <= 1.0))
    

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

    array_ph = nm.Euler_method(alpha, Point, integ_config)                          #Array com h > 0
    array_mh = nm.Euler_method(alpha, Point, [-integ_config[0], integ_config[1]])   #Array com h < 0

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
            bp = 0       #Armazena o indice da branch do Point
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
        
        else: #Se o Point nao estiver no prisma
            #Extraindo os primeiros pontos das branches
            distances = []         #lista contendo as distancias ate o Point em relacao a cada primeiro ponto de cada branch
            branch_index = 0       #indice da branch que possui o ponto mais proximo do Point
            for i in range(len(branches_list)):
                distances.append(np.linalg.norm(branches_list[i][0] - Point))

            branch_index = np.argmin(distances) #Armazena o indice da branch que possui o ponto mais proximo ao Point
            
            branch_atual = branches_list[branch_index]
            index_max = np.argmax(branches_list[branch_index][:, 2]) #Extrai o indice do ponto maximo da branch bp + 1

            Point_max = branch_atual[index_max]
            if(not np.isclose(np.linalg.norm(Point_max - branch_atual[-1]), 0)): #Verifica se o ponto maximo nao eh o ultimo elemento da branch
                branch_1 = branch_atual[: index_max + 1]                         #Fatia a branch a partir do Point ate o o ponto maximo da branch (o incluindo)
                branch_2 = branch_atual[index_max:]                              #Fatia a branch a partir do ponto maximo (o incluinto)

                branches_list.pop(branch_index)                                  #Remove a branch do ponto maximo
                branches_list.insert(branch_index, branch_1)                     #Substitui a branch que estava na posicao bp+1 
                branches_list.insert(branch_index + 1, branch_2)                 #Cria uma nova branch que inicia a partir do ponto maximo

            """
            Observacao: um caso que PODERIA quebrar esta logica seria se o Point estiver proximo a uma outra branch que, visualizando no grafico do Point_Prism.py,
            CLARAMENTE nao seria a branch que conteria o Point se ele estivesse no prisma.
            """

    else:
        #Transforma em um array numpy e recebe os pontos organizados
        branches_list = np.array([org_points])

    return branches_list #retorna como lista, para melhor eficiencia

#A funcao a seguir retorna um array em que cada elemento terah um valor do tipo
# 'b' - lambdaS => Azul
# 'r' - lambdaF => Vermelho
# 'purple' - LambdaZ => Roxo
# 'k' - ERROR
def Branches_point_colors(alpha, branches):
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

#________________________BRANCHES DOS AUTOVALORES________________________#

#A funcao a seguir retornarah um dois ramos:
#   -> Um ramo relacionado a LAMBDA_S == LAMBDA_Z
#   -> Um ramo relacionado a LAMBDA_F == LAMBDA_Z
def Branches_auto(bicetion_config : list, Num_z):

    #Definindo intervalos
    interval_1 = [0, 0.6] #Intervalo para LambdaS e LambdaZ
    interval_2 = [0.4, 1] #Intervalo para LambdaF e LambdaZ

    #Definindo o z inicial (inicio da curva de nivel):
    zk = 0
    dz = 1/Num_z #Diferenca entre duas curvas de nivel

    #Definindo a lista que vai armazenar as Branchs
    Branch_list = [] #Vai armazenar as branches LambdaS com LambdaZ e LambdaF com LambdaZ

    def pointInT(Point):
        if(not transparencia()):
            return 1

        if(baricentrica()):
            return if_PointInEq(Point)
        else:
            return if_PointInRet(Point)

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
                array_temp = array_org[i+1:]                                                     #Cria uma copia do array a partir de um indice i
                return indexs(array_temp, lista_index, bool_value, current_size, correc_error)


    for _ in range(Num_z):
        #Definindo uma lista que vai armazenar TODOS os pontos
        Points_list = []

        points_s = nm.Bisection_method(interval_1, bicetion_config, lbdas, zk) #Armazena os pontos onde LambdaS - LambdaZ == 0
        points_f = nm.Bisection_method(interval_2, bicetion_config, lbdaf, zk) #Armazena os pontos onde LambdaF - LambdaZ == 0
        Points_list = [points_s, points_f]

        #Definindo variavel que vai armazenar as branches da curva de nivel:
        branches_curva = []

        for i in range(2): 
            index_list = []                                     #Lista dos indices
            boolean_value = False                               #Inicia false, para termos controle
            c_size = len(Points_list[i])                        #Tamanho atual do Point_list

            if(pointInT(Points_list[i][0])): #Verifica se o primeiro ponto da lista estah ou nao no prisma (a depender da transparencia)
                boolean_value = not boolean_value

            print(f"\n--- i = {i} ---")
            print(f"Tamanho de Points_list[i]: {c_size}")
            print(f"Primeiros pontos: {Points_list[i][:3]}")
            print(f"boolean_value (dentro do prisma?): {boolean_value}")

            indexs(Points_list[i], index_list, boolean_value, c_size, correc_error = 0)

            print(f"index_list após indexs(): {index_list}")
            print(f"inside inicial: {bool(pointInT(Points_list[i][0]))}")


            #Inicializando variaveis para loop:
            p = 0                                       #Variavel de posicao
            inside = bool(pointInT(Points_list[i][0]))  #Variavel de controle

            #Definindo variavel que vai armazenar as branches da curva de nivel:
            branch_curva = []

            #Isso alterna a escolha dos pontos
            for j in range(len(index_list)):
                if(inside): #Se os pontos estiverem dentro do prisma, ele os seleciona
                    branch_curva.append(Points_list[i][p : index_list[j]])

                #Note que, se inside == False, ele apenas pula para o proximo indice, nao escolhendo os pontos que estao fora
                p = index_list[j]
                inside = not inside

            if(inside):
                branch_curva.append(Points_list[i][p:]) #Adicionando o ultimo segmento se a ultima iteracao fizer Inside = True

            branches_curva.append(branch_curva) #Armazena os lados LambdaS ou LambdaF (com LambdaZ)
        
        Branch_list.append(branches_curva) #Vai armazenar a curva de nivel inteira para z = zk
        
        #Atualizo o zk:
        zk += dz #Proxima iteracao -> nova curva de nivel

    
    #print(f"Branch {1}: {Branch_list[0]}\n")

    return Branch_list
            

#A funcao a seguir retornarah um duas cores relacionadas aos ramos da funcao anterior:
#   -> Se o ramo for do tipo LAMBDA_S == LAMBDA_Z ==> color = marrom escuro (color='#5C4033')
#   -> Se o ramo for do tipo LAMBDA_F == LAMBDA_Z ==> color = marrom claro (color='#C4A484')
def Branches_auto_colors():
    #Definindo variavel para armazenamento das cores:
    colors = []

    #Pela definicao da funcao anterior (via logica da programacao), dado o indice i da Branch,  
    #Branch[i][0] SEMPRE sera marrom escuro e Branch[i][1] SEMPRE sera marrom claro, logo
    colors.append('#5C4033')
    colors.append('#C4A484')

    return colors