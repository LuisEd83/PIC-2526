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
from includes.Auxiliar_Functions import transparencia, if_PointInEq, if_PointInRet, colorPoint, hugoniotSystemSolver
from includes.Campo_Hugoniot import Hug3, F, G

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

    def lambd_comp(Point : list, alpha):#Compara os valores dos lambdas
        #Definindo valores e comparando:
        lambdaS_value = lbdas(*Point)
        lambdaF_value = lbdaf(*Point)
        lambdaZ_value = lambdz(*Point, alpha)

        if(lambdaZ_value < lambdaS_value < lambdaF_value):
            return 1
        elif(lambdaS_value < lambdaZ_value < lambdaF_value):
            return 2
        elif(lambdaS_value < lambdaF_value < lambdaZ_value):
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
    array_ph = array_ph[1:]
    array_mh = array_mh[1:]

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


        if(pointInP(Point)):
            #Localizando a branch e o indice do Point
            bp = 0       #Armazena o indice da branch do Point
            index_bp = 0 #Armazena o indice do Point na branch 
            for i, branch in enumerate(branches_list):
                for j, p in enumerate(branch):
                     if np.all(p == Point):
                        bp = i
                        index_bp = j
            
            if(index_bp != 0): #Verifica se o Point eh o primeiro elemento da branch do Point
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
        #Criando a variaveis para chamar funcao:
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

    #for i in range(len(branches_list)):
    #    print(f"Branch {i+1}: {branches_list[i]} \n")

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
        colors.append(colorPoint(alpha, point)) #Armazena cor do ponto da forma ['LambZ', 'LambS', 'LambF']


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


######################################################3
#Branches relacionadas ao Campo Hugoniot
def Branches_Hugoniot(
        initialValue  : float, #Valor inicial do intervalo da grade
        finalValue    : float, #Valor final do intervalo da grande
        Resol         : int,   #Resolucao da grade
        enableMask    : bool,  #Variavel que habilita o filtro para o triangulo           
        alpha         : float, #Variavel de controle
        fixed_point   : list,  #Ponto fixo
        inteConfig    : list,  #Condiguracao da integracao
        TOL = 1e-3             #Tolerancia
) -> list:
    """
    [Explicacao]
    O que esta funcao realiza?
    -> Integracao a partir de um ponto inicial definido a partir da funcao hugoniotSystemSolve
       para obter arrays (branches) de pontos.
    A funcao retornarah as branches da seguinte forma: [[], [], [], ..., []]


    Possivel adicao futura: Um metodo para correcao da curva (provavelmente Newton-Rapshon)
    """
    import includes.Campo_Hugoniot as ch

    #-------------------------------------------------------------
    #| Extracao de variaveis convenientes para melhor explicacao |
    #-------------------------------------------------------------

    #Extraindo h e N:
    h, N = inteConfig

    #Extraindo o ponto fixo:
    u0, v0, z0 = fixed_point

    #----------------------------------------------------------
    #| Calculo do primeiro conjunto de pontos para integracao |
    #---------------------------------------------------------

    #Defini-se uma variavel para armazenar os pontos iniciais de integracao
    bestPoint1 = None
    bestPoint2 = None #Isso supondo que ha duas raizes no plano z = 0

    #Plano z = constante
    """
    [Explicacao] - NO caso limite, z0 pode ser igual a 0.0, logo deve-se integrar a partir de um
                   plano um pouco acima.
    """
    zinit = 0.0 if (z0 > TOL) else 0.05
    _, _, _, _, candidatas = hugoniotSystemSolver(
                                                initialValue, finalValue, Resol, enableMask,
                                                u0, v0, z0,
                                                zinit, alpha
                                                )
    
    #Raio de exclusão ao redor do primeiro ponto encontrado no plano
    exclusion_radius = 0.035  #ajuste conforme necessario

    if(len(candidatas) > 0):
        #Primeiro ponto: mais proximo do fixed_point

        """
        [Explicacao] - lambda eh uma funcao anonima (uma função sem nome), definida em uma linha.
                       min varre a lista de candidatas em busca de minimizar a distancia euclidiana 
        """
        bestPoint1 = min(
            candidatas,
            key=lambda c: np.linalg.norm([c[0] - u0,
                                          c[1] - v0,
                                          c[2] - z0])
        )

        #Filtra candidatos fora da região de exclusão ao redor do bestPoint1
        candidatas_validas = [
            c for c in candidatas
            if np.linalg.norm([c[0] - bestPoint1[0],
                               c[1] - bestPoint1[1],
                               c[2] - bestPoint1[2]]) > exclusion_radius
        ]

        #Segundo ponto: mais distante do bestPoint1 (fora da regiao de exclusão)
        if(len(candidatas_validas) > 0):
            bestPoint2 = max(
                candidatas_validas,
                key=lambda c: np.linalg.norm([c[0] - bestPoint1[0],
                                            c[1] - bestPoint1[1],
                                            c[2] - bestPoint1[2]])
            )
        else:
            print("Nenhum segundo ponto válido fora da região de exclusão.")
            bestPoint2 = None
    else:
        print("Nenhum candidato encontrado.")
        return []

    #--------------------------------------------------------
    #| Calculo do segundo conjunto de pontos para integracao|
    #-------------------------------------------------------

    """
    [Explicacao] - A ideia aqui eh encontrar outros pontos em um outro plano z = constante
                   tal que z0 < constante =< 1.0 usando uma estrategia semelhante a bis-
                   seccao.
                 - O caso limite aqui eh quando z0 = 1.0... Neste caso eh soh nao procurar
                   outros pontos-solucao do sistema e integrar a partir do plano z = 0.0.
    """
    
    #Busca binaria para encontrar z_const com 2 candidatas (ou confirmar que ha apenas 1 ramo)
    z_low   = z0
    z_high  = 1.0
    z_const = z_high

    MAX_ITER = 10  #Limite de iteracoes da busca binaria
    cand = []

    for _ in range(MAX_ITER):
        if(abs(z_const - z0) <= TOL):
            break

        _, _, _, _, cand = hugoniotSystemSolver(
                                initialValue, finalValue, Resol, enableMask,
                                u0, v0, z0,
                                z_const, alpha
                                )

        if(len(cand) > 1):
            break  #Encontrou 2 candidatas -> define z_const

        #Apenas 1 (ou nenhuma) candidata -> tenta subir z_const
        z_const = (z_low + z_high)/2.0
        z_high = z_const

    #Selecao dos pontos em z = z_const
    bestPoint3 = None
    bestPoint4 = None

    if(len(cand) > 0):
        #Primeiro ponto: mais proximo do fixed_point (Essa eh a estrategia atual)
        bestPoint3 = min(
            cand,
            key=lambda c: np.linalg.norm([c[0] - u0,
                                          c[1] - v0,
                                          c[2] - z0])
        )

        #Filtra candidatos fora da regiao de exclusao ao redor do bestPoint3
        cand_validas = [
            c for c in cand
            if np.linalg.norm([c[0] - bestPoint3[0],
                               c[1] - bestPoint3[1],
                               c[2] - bestPoint3[2]]) > exclusion_radius
        ]

        #Segundo ponto: mais distante do bestPoint3
        if(len(cand_validas) > 0):
            bestPoint4 = max(
                cand_validas,
                key=lambda c: np.linalg.norm([c[0] - bestPoint3[0],
                                              c[1] - bestPoint3[1],
                                              c[2] - bestPoint3[2]])
            )
        else:
            print("Nenhum segundo conjunto de pontos válidos em z = z_const.")

    #--------------------------------
    #Integracao dos ramos inferiores|
    #-------------------------------

    """
    [Explicacao] - Buscando uma solucao geral, pensei em integrar a partir dos pontos determinados anteriormen-
                   te ateh um valor proximo ao plano z = z0, quando o ponto estiver proximo ao plano z = z0, 
                   quebro a integracao.
                   Definido um outro plano z = z1 entre o plano z = z0 e z = 1 (calculado via "busca binaria) e
                   calculado o proximo ponto-solucao (vao surgir alguns, entao devo escolher de forma precisa um
                   deles) e, com este ponto, inicio uma nova integracao em duas direcoes (analogo ao procedi-
                   mento anterior envolvendo o outro campo).
    """

    #Definindo a variavel que armazenarah todas as Branches
    branches = []

    #Definindo a variavel arr para armazenar o array atual
    arr = []

    """
    [Explicacao] - O trecho abaixo define um ou dois ramos inferiores a partir dos pontos encontrados
                   anteriormente.

    [Refinamento] - No lugar de dar um "append" na variavel branches de forma aleatoria, vamos esco-
                    da seguinte forma (supondo que ha dois pontos-solucao):
                    -> Se os sentidos da terceira componente do campo nos pontos forem diferentes, 
                       a prioridade eh o sentido positivo. Ou seja, se Hug3 > 0 em um dos pontos, es-
                       te serah o primeiro ponto a ser integrado e dado "append" na variavel branches.
                    -> Se os sentidos da terceira componente do campo nos pontos forem iguais, escolhe
                       o primeiro ponto.
                    Resumindo: Eh selecionado o primeiro ponto para integracao aquele que apresenta
                    Hug3 > 0. 
                    Note que, desta forma, eh "soh integrar", uma vez que, ao realizar o "append", os
                    array jah estarao organizados.
    """

    #Organizacao dos pontos
    """
    [Explicacao] - Note que a troca soh serah necessaria se, e somente se, o sentido em bestPoint1 for
                   negativo e em bestPoint2 for positivo.
                   -> Se Hug3 > 0 em bestPoint1 e Hug3 > 0 em bestPoint2 => Pontos organizados
                   -> Se Hug3 > 0 em bestPoint1 e Hug3 < 0 em bestPoint2 => Pontos organizados
                   -> Se Hug3 < 0 em bestPoint1 e Hug3 > 0 em bestPoint2 => Pontos desorganizados => Troca
                   -> Se Hug3 < 0 em bestPoint1 e Hug3 < 0 em bestPoint2 => Impossivel (tem que ver na teoria)
    """
    if(len(list(filter(None, [bestPoint1, bestPoint2]))) > 1):
        if(Hug3(alpha, u0, v0, z0, bestPoint1[0], bestPoint1[1], bestPoint1[2]) < 0 and
           Hug3(alpha, u0, v0, z0, bestPoint2[0], bestPoint2[1], bestPoint2[2]) > 0):
            temp = bestPoint1
            bestPoint1 = bestPoint2
            bestPoint2 = temp

    for bestPoint in filter(None, [bestPoint1, bestPoint2]):
        #Define o sentido da integracao:
        """
        [Explicacao] - Fica mais facil definir o sentido da integracao calculando apenas a componente
                       z do campo Hugoniot.
                     - A linha abaixo atribui 1 a variavel sense se o terceiro componente do campo for
                       maior que 0, caso contrario atribui 0
        """
        sense = 1 if (Hug3(alpha, u0, v0, z0, bestPoint[0], bestPoint[1], bestPoint[2]) > 0) else 0

        print(f"F = {ch.F(u0, v0, z0, bestPoint[0], bestPoint[1], bestPoint[2])}")
        print(f"G = {ch.G(alpha, u0, v0, z0, bestPoint[0], bestPoint[1], bestPoint[2])}")
        print("----------------------------------------------------------------------------")

        #Com o sentido, faz-se a integracao
        if(sense):
            arr = nm.HugonioutEuler_method(alpha, fixed_point, bestPoint, [h, N])
        else:
            arr = nm.HugonioutEuler_method(alpha, fixed_point, bestPoint, [-h, N])
        
        #Armazeno os pontos
        branches.append(arr)

    """
    [Observacao] - Agora a variavel branches estah fa seguinte forma:

                    [[Ramo inicial], [Ramo final]]

                   Ao realizar a proxima integracao, deveremos por os ramos entre os ramos inicial e fi-
                   nal de tal forma que os pontos fiquem ordenados. 
    """

    #--------------------------------
    #Integracao dos ramos superiores|
    #-------------------------------

    """
    [Explicacao] - O trecho abaixo define um ou dois ramos a partir dos pontos calculados
                   em um plano z = constante.
                 - Claro, possa ser que nao haja esses ramos, pois, em um caso limite, z0
                   eh justamente 1.0, logo, integrando a partir de z = 0.0, teria a curva
                   integral inteira.
                   
    [Refinamento] - A ordem de integracao vai se dar da seguinte forma:
                    -> Primeiro: Identifico o bestpoint que estah mais proximo do ultimo 
                       ponto do ramo inicial. Integro com h > 0 e h < 0 a partir deste 
                       ponto, organizo em relacao ao valor de z dos pontos gerados. Com
                       os pontos organizados, realizo o append em 'branches' entre o Ramo
                       inicial e final.
                    -> Segundo: O segundo ponto (caso exista) estarah proximo do ultimo 
                       ponto do Ramo final (questao de continuidade), integro com h > 0 e 
                       h < 0 e o ordeno (de forma analoga ao primeiro ponto). Porem ha um
                       detalhe: a integracao dos dois pontos podem gerar a mesma curva, 
                       portanto, deve-se "matar" os pontos duplicados. Com estes pontos, 
                       realizo o append em "branches" entre o ramo da etapa anterior. 
    """
    #Ramos a partir de z = z_const (bestPoint3 e bestPoint4)
    upper_points = [bestPoint3, bestPoint4]

    #Realiza a troca entre os bestpoints se necessário
    valid_upper = list(filter(None, upper_points))
    if len(valid_upper) > 1:
        last = branches[0][-1]  #Ultimo ponto do ramo inicial

        diff3 = np.array([last[0] - bestPoint3[0], last[1] - bestPoint3[1], last[2] - bestPoint3[2]])
        diff4 = np.array([last[0] - bestPoint4[0], last[1] - bestPoint4[1], last[2] - bestPoint4[2]])

        #Se bestPoint4 estiver mais proximo do último ponto do ramo inicial -> troca
        if np.linalg.norm(diff3) > np.linalg.norm(diff4):
            bestPoint3, bestPoint4 = bestPoint4, bestPoint3
            upper_points = [bestPoint3, bestPoint4]

    #-----------------------------------------------------------------------
    #Referencia para checagem de duplicatas: ultimo ponto do ramo acumulado|
    #----------------------------------------------------------------------
    insert_idx = len(branches) - 1  #Posicao antes do Ramo final

    for idx, bestPoint in enumerate(filter(None, upper_points)):

        arr_pos = nm.HugonioutEuler_method(alpha, fixed_point, bestPoint, [ h, N])
        arr_neg = nm.HugonioutEuler_method(alpha, fixed_point, bestPoint, [-h, N])

        arr_pos = arr_pos[1:]
        arr_neg = arr_neg[1:]

        col_pos = arr_pos[:, 2]
        col_neg = arr_neg[:, 2]

        i = 0
        while(True):
            dz = col_pos[0] - col_neg[i]

            if(dz > 0):
                if np.linalg.norm(arr_neg[-1] - bestPoint) > np.linalg.norm(arr_neg[0] - bestPoint):
                    arr_neg = np.flip(arr_neg, axis=0)
                org_points = np.concatenate([arr_neg, np.array(bestPoint).reshape(1, -1), arr_pos])
                break

            elif(dz < 0):
                if np.linalg.norm(arr_pos[-1] - bestPoint) > np.linalg.norm(arr_pos[0] - bestPoint):
                    arr_pos = np.flip(arr_pos, axis=0)
                org_points = np.concatenate([arr_pos, np.array(bestPoint).reshape(1, -1), arr_neg])
                break

            else:
                i += 1

            if(i == len(col_neg)):
                print("[ERROR] - Impossibilidade de determinar ordem")
                exit()

        #Remocoo de "duplicatas": compara contra o org_points do bestPoint anterior
        TOL_dup = max(TOL, abs(h) * 2)  # Tolerancia especifica para deduplicacao

        if idx > 0:
            prev_branch = branches[insert_idx - 1]

            mask = np.ones(len(org_points), dtype=bool)
            for j, pt in enumerate(org_points):
                if np.any(np.linalg.norm(prev_branch - pt, axis=1) <= TOL_dup):
                    mask[j] = False

            org_points = org_points[mask]

        if(len(org_points) == 0):
            continue #Pula, nao ha nada a inserir (curva totalmente duplicada)

        #Insere antes do Ramo final
        if(len(branches) < 2):
            #Nao ha ramo final definido (No caso, ramo inicial = ramo final), "appenda" diretamente
            branches.append(org_points)
        else:
            branches.insert(insert_idx, org_points)
        insert_idx += 1  #Proximo insert empurra mais um

    #[Guard] Mata as branches vazias antes de prosseguir
    branches = [b for b in branches if len(b) > 0]

    #---------------------------------------------
    #|Filtro os pontos que estao dentro do Prisma|
    #---------------------------------------------
    
    for i in range(len(branches)):                                  #Varre todas as branches
        mask = np.ones(len(branches[i]), dtype=bool)                #Reinicializa para cada branch
        for j, pt in enumerate(branches[i]):                        #Varre a i-esima branch
            if((not pointInP(pt)) or                                #Verifica se nao esta no prisma
               ((abs(F(u0, v0, z0, pt[0], pt[1], pt[2])) > TOL) and
               (abs(G(alpha, u0, v0, z0, pt[0], pt[1], pt[2])) > TOL))):                     
                mask[j] = False                                     #"Mata" o j-esimo ponto do i-esimo branch

        branches[i] = branches[i][mask]

    #[Guard] Mata as branches vazias antes de prosseguir
    branches = [b for b in branches if len(b) > 0]

    #------------------------------------------------------------------------------
    #|Qubra de Branches onde pontos consecutivos que estao distantes (ateh demais)|
    #-----------------------------------------------------------------------------

    #Define-se uma distancia maxima entre os pontos consecutivos
    DIST_MAX = abs(h) * 5
    new_branches = []   #Update das branches


    """
    [Explicacao] - Este algorismo funciona da seguinte forma:
                   -> 'Identifica' a i-esima branch e armazena;
                   -> Branch_temp eh iniciada com o primeiro ponto da i-esima branch
                   -> Varre a i-esima branch em busca de dois pontos que estao suficientemente distantes.
                      Enquanto nao encontra, a branch_temp acumula os pontos e, ao encontrar o j-esimo
                      ponto que esta distante do seu sucessor, new_branches armazena os pontos acumulados
                      pela branch_temp e reinicia branch_temp com o sucessor do j-esimo ponto.
    """
    for i in range(len(branches)):
        branch_atual = list(branches[i])          #Trabalha como lista de pontos
        branch_temp  = [branch_atual[0]]          #Inicia com o primeiro ponto da i-esima branch

        """
        [Explicacao] - Este trecho de codigo, como foi dito anteriormente, varre a i-esima branch em busca
                       de dois pontos que estao suficientemente distantes e, sabendo a sua localizacao na 
                       branch, ha a quebra da branch.
        """
        for j in range(len(branch_atual) - 1):
            dist = np.linalg.norm(np.array(branch_atual[j+1]) - np.array(branch_atual[j]))

            if((dist > DIST_MAX)):
                #Quebra: salva o trecho atual e inicia um novo
                new_branches.append(np.array(branch_temp))  #Armazena os pontos que foram acumulados
                branch_temp = [branch_atual[j+1]]           #Reinicia com o sucessor j-esimo ponto que demarca a quebra
            else:
                branch_temp.append(branch_atual[j+1]) #Acumula o ponto na variavel

        if(len(branch_temp) > 0): #Se branch_temp nao for nulo, i.e, se houver pontos "sobrando", appenda em new_branches
            new_branches.append(np.array(branch_temp))

    #Atualizacao das branches
    branches = new_branches
    
    #--------------------------------------------------------------------------------------------------------
    #|Qubra de Branches onde pontos estao localizados em um maximo local (com concavidade virada para baixo)|
    #-------------------------------------------------------------------------------------------------------
    """
    [Explicacao] - De acordo com a teoria, no ponto de maximo local (onde a concavidade esta virada para baixo)
                   ha a troca dos valores dos autovalores, i.e, ocorre 

    [Refinamento] - Como a quebra da i-esima branch gera duas novas branches, ocorre o aumento no len de branches.
                    Para contornar isto o for deve percorrer branch a branch.
    """
    #Variavel para novas branches
    new_branches = [] #Para reaproveitar variavel

    for branch_atual in branches:
        if len(branch_atual) == 0:
            continue  #Pula branches vazias

        max_index = np.argmax(branch_atual[:, 2]) #Localiza o ponto com o maior z na branch
        """
        [Temporário] - teste de lambdas
        """
        print("######################################")
        print(f"Coordenadas do ponto maximo: {branch_atual[max_index][0], branch_atual[max_index][1], branch_atual[max_index][2]}")
        print(colorPoint(alpha, branch_atual[max_index]))
        print(f"Valores dos lambdas: {lbdas(*branch_atual[max_index]), lambdz(*branch_atual[max_index], alpha), lbdaf(*branch_atual[max_index])}")

        
        #Verifica-se se este ponto nao eh o ultimo na branch e se ele nao esta muito proximo do plano z = 1.0 
        if((not np.isclose(branch_atual[-1], branch_atual[max_index]).all()) and (abs(branch_atual[max_index][2] - 1.0) > TOL)):
            #Quebra da branch no index do ponto maximo
            branch1 = branch_atual[:max_index + 1]
            branch2 = branch_atual[max_index:]

            new_branches.append(branch1)
            new_branches.append(branch2)
        else:
            #Appenda o branch inteiro, sem realizar a quebra
            new_branches.append(branch_atual)
    
    branches = new_branches

    return branches