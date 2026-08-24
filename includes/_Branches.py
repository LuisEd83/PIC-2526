"""
Modulo: Branches

Objetivo:
- Criar uma função que retorna ramos (e sua respectiva coloracao)

Os ramos devem:
- Ter uma coloracao de acordo com o autovalor lambdaS, lambdaF e lambdaZ  

"""

import includes.Numericals_methods as nm
from includes.Inicia import baricentrica
from includes.Functions import fw, lbdas, lbdaf, map
from includes.Auxiliar_Functions import transparencia, if_PointInRet, colorPoint, az
from includes.Campo_Ez import R

import numpy as np

#________________________BRANCHES DA CURVA DE RAREFACAO________________________#

#Redefinindo funcoes para nao haver erro circular
def lambdz(u, v, z, alpha):
    return(fw(u, v, z)/(u + alpha*az(z)))

def pointInP(Point):
    if(not transparencia()):
        return 1
    
    return (if_PointInRet(Point) and (0.0 <= Point[2] <= 1.0))

def Branches_point(
        alpha           : float,   #Variavel de controle
        Point           : list,    #Ponto inicial
        integ_config    : list,    #Lista de configuracao para integracao
):
    #_________DEFINICAO DE FUNCOES AUXILIARES PARA ESTA FUNCAO________#
    #A funcao abaixo define a cor de cada ramo em funcao do lambdaZ
    def Branches_point_colors(alpha, branches):
        colors = []
        for i in range(len(branches)):
            point = branches[i][len(branches[i])//2]
            colors.append(colorPoint(alpha, point))
        return colors

    def indexs(array_org, bool_value_inicial):
        """
        Retorna a lista de indices onde ocorre transicao dentro/fora do prisma.
        bool_value_inicial deve ser pointInP(array_org[0]).
        """
        lista_index = []
        bool_value = bool_value_inicial
        for i in range(len(array_org)):
            if(pointInP(array_org[i]) == (not bool_value)):
                bool_value = not bool_value
                lista_index.append(i)
        return lista_index

    def filtraBranchesPeloPrisma(branches):
        novos_branches = [] #Variavel para armazenar as novas branches filtradas no prisma

        for branch in branches:                 #Itera sobre cada branch
            if(len(branch) == 0):
                continue

            estado_inicial = pointInP(branch[0])        #Estado (dentro/fora) do primeiro ponto
            transicoes = indexs(branch, estado_inicial)

            #Monta os limites dos segmentos: inicio em 0, fim em len(branch)
            limites = [0] + transicoes + [len(branch)]

            #Cada segmento k comeca com estado = estado_inicial, alternando a cada transicao
            estado_segmento = estado_inicial
            for k in range(len(limites) - 1):
                a, b = limites[k], limites[k+1]

                if(estado_segmento):              #Mantem apenas segmentos DENTRO do prisma
                    segmento = branch[a:b]
                    if(len(segmento) > 1):
                        novos_branches.append(segmento)

                estado_segmento = not estado_segmento  #Alterna para o proximo segmento

        return novos_branches

    #===========================
    #EXTRACAO DOS DADOS INICIAIS
    #===========================

    u0, v0, z0 = Point          #Coordenadas do ponto inicial
    h, N       = integ_config   #Configuracao para integracao

    #================================
    #CRIACAO E ORGANIZACAO DOS PONTOS
    #================================

    points_ph = nm.Euler_method(alpha, Point, [h, N])   #Ramo de integracao para h > 0
    points_mh = nm.Euler_method(alpha, Point, [-h, N])  #Ramo de integracao para h < 0

    #Retirada do ponto inicial (ele eh recolocado depois)
    points_ph = points_ph[1:]   
    points_mh = points_mh[1:]   

    #Para nao necessitar de recalculo futuramente, define-se uma flag de sentido de campo
    #OBS: isso NAO FUNCIONA se R == 0 (ao menos atualmente)
    fieldSense = True if(R(u0, v0, z0, alpha) > 0) else False

    #Para armazenar os ramos, define-se uma lista de armazenamento
    branches = []   #lista vazia

    #Realizo a concatenacao dos pontos
    if(fieldSense):
        points_mh = np.flip(points_mh, axis = 0)          #Inverte o array inferior
        points_mh = np.append(points_mh, [Point], axis = 0) #Insiro o ponto inicial no array inferior
        branches.extend([points_mh, points_ph]) #Armazeno os ramos
    else:
        points_ph = np.flip(points_ph, axis = 0)          #Inverte o array inferior
        points_ph = np.append(points_ph, [Point], axis = 0) #Insiro o ponto inicial no array inferior
        branches.extend([points_ph, points_mh]) #Armazeno os ramos


    #=====================
    #FILTRAGEM PELO PRISMA
    #=====================
    branches = filtraBranchesPeloPrisma(branches)

    #=================================================================
    #CORTES NAS BRANCHES EM FUNCAO NAS DIFERENCAS ENTRE OS AUTOVALORES
    #=================================================================

    #Defina uma nova branch para armazenar os novos ramos (que serao gerados pelos cortes)
    new_branches = []

    for branch in branches:
        #Encontrar todos os pontos de corte no branch atual
        cut_indices = []

        #Inicializamos as diferencas do primeiro ponto (sem abs)
        uk, vk, zk = branch[0]
        prev_diff_s = lbdas(uk, vk, zk) - lambdz(uk, vk, zk, alpha)
        prev_diff_f = lbdaf(uk, vk, zk) - lambdz(uk, vk, zk, alpha)

        #Iteramos a partir do segundo ponto (j = 1)
        for j in range(1, len(branch)):
            u, v, z = branch[j]
            
            diff_s = lbdas(u, v, z) - lambdz(u, v, z, alpha)
            diff_f = lbdaf(u, v, z) - lambdz(u, v, z, alpha)

            #Se o produto for < 0, houve mudanca de sinal entre o ponto j-1 e j
            if ((diff_s * prev_diff_s) < 0) or ((diff_f * prev_diff_f) < 0):
                cut_indices.append(j)
                
            #Atualiza as diferenças anteriores para a proxima iteracao
            prev_diff_s = diff_s
            prev_diff_f = diff_f
        
        #Se nao houver cortes, mantem a ramificacao intacta
        if(not cut_indices):
            new_branches.append(branch)
            continue
        
        #Fatiar o branch nos pontos encontrados
        start_idx = 0
        for cut_idx in cut_indices:
            #Inclui o trecho ate o ponto de corte (inclusive)
            sub_branch = branch[start_idx : cut_idx + 1]
            if(len(sub_branch) > 1):
                new_branches.append(sub_branch)
            start_idx = cut_idx + 1
        
        #Adiciona o restante da ramificação após o último corte
        remaining = branch[start_idx:]
        if(len(remaining) > 1):
            new_branches.append(remaining)

    #Atualiza a lista original apos terminar o processamento
    branches = new_branches

    #========================
    #COLORACAO DE CADA BRANCH
    #========================

    colors = Branches_point_colors(alpha, branches)

    return colors, branches

#________________________BRANCHES DOS AUTOVALORES________________________#

#A funcao a seguir retornarah um dois ramos:
#   -> Um ramo relacionado a LAMBDA_S == LAMBDA_Z
#   -> Um ramo relacionado a LAMBDA_F == LAMBDA_Z
def Branches_auto(bicetion_config : list, Num_z, alpha, base, altura):

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
        else:
            for branch in [points_s, points_f]:
                for point in branch:
                    u, v, z = point
                    u, v = map(u, v, baricentrica())

                    point[0] = u
                    point[1] = v
        
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

    #===========================
    #CORES DOS RAMOS SLOW E FAST
    #==========================
    colors = Branches_auto_colors()

    return colors, Branch_list
            

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
        alpha         : float, #Variavel de controle
        fixed_point   : list,  #Ponto fixo
        integ_config  : list,  #Condiguracao da integracao
        initialValue  : float, #Valor inicial do intervalo da grade
        finalValue    : float, #Valor final do intervalo da grande
        Resol         : int,   #Resolucao da grade
        enableMask    : bool,  #Variavel que habilita o filtro para o triangulo
        TOL = 1e-3               #Tolerancia
) -> list:
    from includes.Auxiliar_Functions import HugoniotPlanaImplicita, hugoniotSystemSolver
    from includes.Campo_Hugoniot import sig, sigm

    def indexs(array_org, bool_value_inicial):
        """
        Retorna a lista de indices onde ocorre transicao dentro/fora do prisma.
        bool_value_inicial deve ser pointInP(array_org[0]).
        """
        lista_index = []
        bool_value = bool_value_inicial
        for i in range(len(array_org)):
            if(pointInP(array_org[i]) == (not bool_value)):
                bool_value = not bool_value
                lista_index.append(i)
        return lista_index


    def filtraBranchesPeloPrisma(branchesHugoniot):
        novos_branches = [] #Variavel para armazenar as novas branches filtradas no prisma

        for branch in branchesHugoniot:                 #Itera sobre cada branch
            points = [p[0] for p in branch]             #Extrai apenas [u, v, z] de cada ponto do branch

            if(len(points) == 0):
                continue

            estado_inicial = pointInP(points[0])        #Estado (dentro/fora) do primeiro ponto
            transicoes = indexs(points, estado_inicial)

            #Monta os limites dos segmentos: inicio em 0, fim em len(branch)
            limites = [0] + transicoes + [len(branch)]

            #Cada segmento k comeca com estado = estado_inicial, alternando a cada transicao
            estado_segmento = estado_inicial
            for k in range(len(limites) - 1):
                a, b = limites[k], limites[k+1]

                if(estado_segmento):              #Mantem apenas segmentos DENTRO do prisma
                    segmento = branch[a:b]
                    if(len(segmento) > 0):
                        novos_branches.append(segmento)

                estado_segmento = not estado_segmento  #Alterna para o proximo segmento

        return novos_branches

    def colorHugoniot(
        alpha,
        branchesHugoniot
    ):
        #Variavel que armazenara as cores de cada branch:
        colors = []

        for branch in branchesHugoniot:
            if not branch:
                continue

            point = branch[len(branch)//2][0]

            color_list = colorPoint(alpha, point)
            colors.append(color_list[0])

        return colors

    #Extraindo o ponto inicial
    u0, v0, z0 = fixed_point

    #Extraindo configuracoes de integracao
    h, N = integ_config

    #conjunto de pontos para h>0
    pointPH = nm.HugonioutEuler_method(alpha, fixed_point, integ_config)
    pointMH = nm.HugonioutEuler_method(alpha, fixed_point, [-h, N])

    #Capta o sentido do campo no ponto inicial
    fieldSense = True if R(u0, v0, z0, alpha) > 0 else False

    #Variavel para armazenar o ultimo ponto do conjunto de pontos:
    lastPoint = pointPH[-1] if(fieldSense) else pointMH[-1]

    #Variavel que guarda as branches:
    branchesHugoniot = []

    #Calcula os pontos no plano z = z0
    branchesImplicita, pointsImplicita = HugoniotPlanaImplicita(alpha, u0, v0, z0, TOL)

    #Supondo que a curva Hugoniot VOLTE
    if((fieldSense and (abs(pointPH[-1][0][2] - z0) < TOL)) or ((not fieldSense) and (abs(pointMH[-1][0][2] - z0) < TOL))):

        todos_pontos = [
            (idx, pt) for idx, raizes in pointsImplicita for pt in raizes
        ]

        _, best_pt = min(todos_pontos, key = lambda item: np.sqrt(
            (item[1][0] - lastPoint[0][0])**2 +  # item[1][0] -> u
            (item[1][1] - lastPoint[0][1])**2 +  # item[1][1] -> v
            (z0 - lastPoint[0][2])**2
        ))

        bestPoint = [best_pt[0], best_pt[1], z0]

        #Anexa ao conjunto correto dado o sentido do campo:
        if(fieldSense):
            pointPH.append([bestPoint, sig(alpha, z0, u0, v0, z0)])
        else:
            pointMH.append([bestPoint, sig(alpha, z0, u0, v0, z0)])

        #Variavel de altura (zmax) para aplicar o Solver da hugoniot (com z != z0)
        hmax = pointMH[1][0][2] if(fieldSense) else pointPH[1][0][2]

        #Aplica o solver da Hugoniot para encontrar o proximo ponto para integracao
        _, _, _, _, pointsSolver = hugoniotSystemSolver(initialValue, finalValue, Resol, enableMask, u0, v0, z0, hmax, alpha, TOL = TOL)

        #Ponto mais proximo do solver
        newPoint = min(pointsSolver, key = lambda p: np.sqrt(
                        (p[0] - bestPoint[0])**2 +
                        (p[1] - bestPoint[1])**2 +
                        (p[2] - bestPoint[2])**2
                        ))

        #Calculamos o valor de sigma deste ponto
        sigNewPoint = sigm(alpha, newPoint[0], newPoint[1], newPoint[2], u0, v0, z0)

        #A partir deste novo ponto, integramos e encontramos o novo ramo
        pointsNewPoint = []
        if(fieldSense):
            pointsNewPoint = nm.HugonioutEuler_method(alpha, fixed_point, [-h, N], [newPoint[0], newPoint[1], newPoint[2], sigNewPoint], False)
        else:
            pointsNewPoint = nm.HugonioutEuler_method(alpha, fixed_point, integ_config, [newPoint[0], newPoint[1], newPoint[2], sigNewPoint], False)


        #======================
        # Motagem das branches
        #======================

        #Ordenacao de branches
        if(fieldSense):
            pointMH = pointMH[1:]
            branchesHugoniot.append(pointMH[::-1])
            pointPH = pointPH[:-35]
            branchesHugoniot.append(pointPH)
        else:    
            pointPH = pointPH[1:]
            branchesHugoniot.append(pointPH[::-1])
            pointMH = pointMH[:-35]
            branchesHugoniot.append(pointMH)

        #Apendda o ultimo ramo
        branchesHugoniot.append(pointsNewPoint)

    else: #Se a curva nao voltar para z = z0

        #======================
        # Motagem das branches
        #======================
        
        if(fieldSense):
            pointMH = pointMH[1:]
            branchesHugoniot.append(pointMH[::-1])
            branchesHugoniot.append(pointPH)
        else:    
            pointPH = pointPH[1:]
            branchesHugoniot.append(pointPH[::-1])
            branchesHugoniot.append(pointMH)

    #Cortes
    new_branches = []

    #Tolerancia de corte
    cutTOL = 5e-4

    for branch in branchesHugoniot:
        #Encontrar todos os pontos de corte no branch atual
        cut_indices = []

        uk, vk, zk = branch[0][0]
        prev_diff_s = lbdas(uk, vk, zk) - lambdz(uk, vk, zk, alpha)
        prev_diff_f = lbdaf(uk, vk, zk) - lambdz(uk, vk, zk, alpha)

        #Iteramos a partir do segundo ponto (j = 1)
        for j in range(1, len(branch)):
            u, v, z = branch[j][0]
            
            diff_s = lbdas(u, v, z) - lambdz(u, v, z, alpha)
            diff_f = lbdaf(u, v, z) - lambdz(u, v, z, alpha)

            #Se o produto for < 0, houve mudanca de sinal entre o ponto j-1 e j
            if ((diff_s * prev_diff_s) < 0) or ((diff_f * prev_diff_f) < 0):
                cut_indices.append(j)
                
            #Atualiza as diferencas anteriores para a proxima iteracao
            prev_diff_s = diff_s
            prev_diff_f = diff_f
        
        #Se não houver cortes, mantém a ramificação intacta
        if(not cut_indices):
            new_branches.append(branch)
            continue
        
        #Fatiar o branch nos pontos encontrados
        start_idx = 0
        for cut_idx in cut_indices:
            #Inclui o trecho ate o ponto de corte (inclusive)
            sub_branch = branch[start_idx : cut_idx + 1]
            if sub_branch:
                new_branches.append(sub_branch)
            start_idx = cut_idx + 1
        
        #Adiciona o restante da ramificacao apos o ultimo corte
        remaining = branch[start_idx:]
        if(remaining):
            new_branches.append(remaining)

    #Atualiza a lista original apos terminar o processamento
    branchesHugoniot = new_branches

    #Filtra as branches para que tenham um numero minimo de pontos
    num_min = 5
    branchesHugoniot = [branch for branch in branchesHugoniot if(len(branch) >= num_min)]

    #Variavel que armazena as cores:
    colors = colorHugoniot(alpha, branchesHugoniot)

    #Filtra para os pontos DENTRO do prisma
    branchesHugoniot = filtraBranchesPeloPrisma(branchesHugoniot)

    return branchesImplicita, colors, branchesHugoniot
