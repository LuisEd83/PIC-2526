"""
Modulo: Branches

Objetivo:
- Criar uma função que retorna ramos (e sua respectiva coloracao)

Os ramos devem:
- Ter uma coloracao de acordo com o autovalor lambdaS, lambdaF e lambdaZ  

"""

import includes.Numericals_methods as nm
from includes.Inicia import baricentrica
from includes.Functions import fw, lbdas, lbdaf 
from includes.Auxiliar_Functions import transparencia, if_PointInEq, if_PointInRet

import numpy as np

#________________________BRANCHES DOS PONTOS________________________#

#Redefinindo funcoes para nao haver erro circular
def lambdz(u, v, z, alpha):
    return(fw(u, v, z)/(u + alpha*np.cos(z)))

def pointInP(Point):
    if(not transparencia()):
        return 1
    if(baricentrica()):
        return (if_PointInEq(Point) and (0.0 <= Point[2] <= 1.0))
    else:
        return (if_PointInRet(Point) and (0.0 <= Point[2] <= 1.0))

def Branches_point(alpha, Point : list, integ_config : list):

    def indexs(array_org, lista_index : list, bool_value : bool, current_size, correc_error):
        #Condicao de parada
        if(current_size <= 0):
            return

        #Criando array para armazenar parte da lista
        array_temp = np.array([])

        for i in range(len(array_org)):
            if(pointInP(array_org[i]) == int(not bool_value)):
                bool_value = not bool_value                                                      #Inverte o valor do bool_value
                lista_index.append(correc_error + i)                                             #Armazena o index que ocorreu o if e soma com o erro (advindo da recursividade)
                correc_error += (i+1)                                                            #Adiciona o erro atual com o antigo
                current_size -= (i+1)                                                            #Retiro parte do tamanho do array
                array_temp = array_org[i+1:].copy()                                              #Cria uma copia do array a partir de um indice i
                return indexs(array_temp, lista_index, bool_value, current_size, correc_error)

    def lambd_comp(Point : list, alpha):
        #Definindo valores e comparando:
        lambdaS_value = lbdas(*Point)
        lambdaF_value = lbdaf(*Point)
        lambdaZ_value = lambdz(*Point, alpha)

        if((lambdaZ_value < lambdaF_value < lambdaS_value) or (lambdaZ_value < lambdaS_value < lambdaF_value)):
            return 1
        elif((lambdaS_value < lambdaZ_value < lambdaF_value) or (lambdaF_value < lambdaZ_value < lambdaS_value)):
            return 2
        elif((lambdaF_value < lambdaS_value < lambdaZ_value) or (lambdaS_value < lambdaF_value < lambdaZ_value)):
            return 3
        else:
            return -1

    def indexs_ccurv(lista_index : list, array_org, color_value, alpha, current_size, correc_error):
        #condição de parada
        if(current_size <= 0):
            return
        
        #Criando array para armazenar parte da lista
        array_temp = np.array([])

        for i in range(len(org_points)):
            if(lambd_comp(org_points[i], alpha) != color_value):
                color_value = lambd_comp(org_points[i], alpha)
                lista_index.append(correc_error + i)
                correc_error += (i+1)
                current_size -= (i+1)
                array_temp = array_org[i+1:].copy()
                return indexs(array_temp, lista_index, color_value, current_size, correc_error)

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
        
        if(pointInP(org_points[0])):    #Determina valor booleano para funcao
            bool_value = not bool_value
        
        #Retorna os indices onde ocorre a transicao de "dentro do prisma" para "fora do prisma" e vice-versa
        indexs(org_points, points_indexs, bool_value, size_array, correc_error = 0)

        #Inicializando variaveis para loop:
        p = 0                                   #Variavel de posicao
        inside = bool(pointInP(org_points[0]))  #Variavel de controle

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

        if(pointInP(Point)):
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
        #Criando a variaveis chamar funcao:
        points_indexs = []                                              #Armazenamento de indices
        color_value = lambd_comp(org_points[0], alpha)                  #Variavel inteira para controle
        size_array = len(org_points)                                    #Variavel para armazenar tamanho do array dos pontos organizados
        
        #Retorna os indices onde ocorre a transicao de "dentro do prisma" para "fora do prisma" e vice-versa
        indexs_ccurv(points_indexs, org_points, color_value, alpha, size_array, correc_error = 0)

        #Inicializando variaveis para loop:
        p = 0                                   #Variavel de posicao

        for i in range(len(points_indexs)):
            branches_list.append(org_points[p : points_indexs[i] + 1])
            p = points_indexs[i]
        
        branches_list.append(org_points[p:])

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
            colors.append('magenta') #Adiciona a cor Vermelha
        elif((lambdaF_value < lambdaS_value < lambdaZ_value) or (lambdaS_value < lambdaF_value < lambdaZ_value)):
            colors.append('r') #Adiciona a cor Roxa
        else:
            colors.append('k')

    return colors

#________________________BRANCHES DOS AUTOVALORES________________________#

#A funcao a seguir retornarah um dois ramos:
#   -> Um ramo relacionado a LAMBDA_S == LAMBDA_Z
#   -> Um ramo relacionado a LAMBDA_F == LAMBDA_Z
def Branches_auto(bicetion_config : list, Num_z, alpha, base, altura):

    #Definindo intervalo
    interval = [0, 1.0]

    #Definindo o z inicial (inicio da curva de nivel):
    zk = base
    dz = (altura - base)/Num_z #Diferenca entre duas curvas de nivel

    #Definindo a lista que vai armazenar as Branchs
    Branch_list = [] #Vai armazenar as branches LambdaS com LambdaZ e LambdaF com LambdaZ

    sqr3 = np.sqrt(3) #Constante
    def if_PointT(Point : list):
        if(not transparencia()):
            return 1
        
        #Extrai as componentes do ponto:
        u, v, z = Point

        if(baricentrica()):
            #______USANDO O TEOREMA DE VIVIANI______#

            #Definindo as distancias entre o ponto P e as arestas do triangulo equilatero
            h1 = np.abs(v)
            h2 = (np.abs(-sqr3 * u + v))/2
            h3 = (np.abs(sqr3 * u + v - sqr3))/2

            #Altura total
            h = sqr3/2

            #Relizando comparacao
            if(np.abs(h1 + h2 + h3 - h) <= 1e-6):
                return 1
        else:
            if((0 <= u <= 1) and (0 <= v <= 1) and (0 <= u + v <= 1)):
                return 1
            
        return 0

    for _ in range(Num_z):
        #Definindo uma lista que vai armazenar TODOS os pontos
        Points_list = []

        points_s = nm.Bisection_method(interval, bicetion_config, lbdas, zk, alpha) #Armazena os pontos onde LambdaS - LambdaZ == 0
        points_f = nm.Bisection_method(interval, bicetion_config, lbdaf, zk, alpha) #Armazena os pontos onde LambdaF - LambdaZ == 0
        
        points_s.insert(0, [0.0, 0.0, zk])

        if(not baricentrica()):
            points_s.append([0.0, 1.0, zk])

        Points_list = [points_s, points_f]

        #Definindo variavel que vai armazenar as branches da curva de nivel:
        branches_curva = []

        for i in range(len(Points_list)): 
            filter_list = []                       #Lista que vai armazenar os pontos filtrados
            
            for j in range(len(Points_list[i])):
                if(if_PointT(Points_list[i][j])):  #Verifica se o ponto esta no triangulo
                    filter_list.append(Points_list[i][j])

            if(len(filter_list) >= 2): #So adiciona se tiver pontos suficientes
                branches_curva.append(filter_list) #Armazena os lados LambdaS ou LambdaF (com LambdaZ)
        
        if(branches_curva): #So adiciona a curva se tiver pelo menos uma branch valida
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

#/////////////////////////////#/////////////////////////////#
#Funcao que marca os pontos a cada x passos a partir do ponto inicial
def andar_branches(Point, branches, passo):
    def ponto_igual(p1, p2, tol=1e-12):
        return np.allclose(p1, p2, atol=tol)

    # Tenta localizar o ponto exato nas branches
    local_point = []
    for i, branch in enumerate(branches):
        for j, p in enumerate(branch):
            if ponto_igual(p, Point):
                local_point = [i, j]
                break
        if local_point:
            break

    # Se não encontrou, busca o mais próximo
    if not local_point:
        melhor_i, melhor_j = None, None
        melhor_dist = float('inf')
        for i, branch in enumerate(branches):
            for j, p in enumerate(branch):
                dist = np.linalg.norm(np.array(p) - np.array(Point))
                if dist < melhor_dist:
                    melhor_dist = dist
                    melhor_i, melhor_j = i, j
        local_point = [melhor_i, melhor_j]

    # Achata apenas os índices
    flat = [(i, j) for i, branch in enumerate(branches) for j in range(len(branch))]

    # Localiza o índice global
    inicio_global = next(
        (k for k, (i, j) in enumerate(flat) if i == local_point[0] and j == local_point[1]),
        None
    )

    # Seleciona índices a cada 'passo', excluindo o ponto de partida
    candidatos = flat[inicio_global % passo::passo]
    candidatos = [(i, j) for (i, j) in candidatos 
                  if not (i == local_point[0] and j == local_point[1])]

    # Filtra conforme transparencia
    point_indexes = []
    if transparencia:
        for (i, j) in candidatos:
            if pointInP(branches[i][j]):
                point_indexes.append((i, j))
    else:
        point_indexes = candidatos

    return point_indexes