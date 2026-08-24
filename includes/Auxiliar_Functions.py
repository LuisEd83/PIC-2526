"""
Modulo: Funcoes auxiliares

Objetivo:
- Armazenar as funcoes que serao utilizadas no codigo e selecionar as necessarias

"""

import includes.Numericals_methods as nm
import includes.Functions as fun

import numpy as np

#Definindo funcao para determinar se o usuário quer ver ou nao os pontos dentro do prisma
def transparencia():
    x = 1 #x = 1 para ver apenas os pontos dentro do prisma; x = 0 e para ver todos os pontos
    return x

#Definindo uma funcao responsavel por permitir a plotagem dos ramos relacionados a curva de
#nivel lambda_s = lambda_z
def branchSlow(): 
    s = 0 #s = 1 para permitir que seja plotado o ramo
    return s

#Definindo uma funcao responsavel por permitir a plotagem dos ramos relacionados a curva de
#nivel lambda_f = lambda_z
def branchFast(): 
    f = 1 #s = 1 para permitir que seja plotado o ramo
    return f

def points_curv(): #Esta funcao habilita os pontos na curva integral.
    x = 1   #x = 1 para habilitar
    return x

##############################################################

"""
def a(z):
    return np.sin(z)

def az(z):
    return np.cos(z)

def a(z):
    return np.arctan(z)
def az(z): #Da/dz
    return 1/(1+z**2)
"""

#Funcao sigmoide
#Supondo k pequeno:
k = 0.01
def a(z):
    return 1/(1+np.exp(-k*z)) - 0.5

def az(z):
    return k*a(z)*(1-a(z))

def lambdz(u, v, z, alpha):
    return(fun.fw(u, v, z)/(u + alpha*az(z)))

#Derivada direcional de lambdZ com direcao dos campos
# def Dlambdz(u, v, z, sigma, u0, v0, z0, alpha):
#     from includes.Campo_Hugoniot import campo

#     #Calculando gradiente:
#     pc = (fun.Du(u, v, z) * (u + alpha*az(z)) - fun.fw(u, v, z))/((u+alpha*az(z))**2)                       #Primeira componente
#     sc = (fun.Dv(u, v, z))/(u + alpha*az(z))                                                                #Segunda componente
#     tc = (fun.Dc(u, v, z) * (u + alpha*az(z)) + fun.fw(u, v, z) * (alpha*np.sin(z)))/((u + alpha*az(z))**2) #Terceira componente

#     #Extraindo campos (componentes de direção)
#     camps = campo(alpha, u, v, z, sigma, u0, v0, z0)

#     #O return sera o produto intero entre o gradiente e a direcao da curva integral (que eh, no fim das contas, os campos)
#     return (pc*camps[0] + sc*camps[1] + tc*camps[2])

#def color_point(colors : list, Point, alpha): #Determina a cor do ponto a partir da derivada do Lambda_z
#    if(Dlambdz(*Point, alpha) > 0):
#        colors.append('r')
#    elif(Dlambdz(*Point, alpha) < 0):
#        colors.append('b')
#    else:
#        colors.append('white')

##############################################################
#Definindo uma constante para a funcao posterior
sqr3 = np.sqrt(3)

#Como está sendo utilizado o diagrama ternario
#Então:
def if_PointInEq(Point : list): #Retorna 1 se verdadeiro
    #Extraindo pontos
    u, v, z = Point

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
    
    return 0

def if_PointInRet(Point : list): #Retorna 1 se verdadeiro
    #Extraindo pontos
    u, v, z = Point

    if((0 <= u <= 1) and (0 <= v <= 1) and (0 <= u + v <= 1)):
        return 1

    return 0

#Definindo uma funcao para concatenar os pontos para h > 0 e h < 0
def Array_Concatenated(alpha, Point : list, integ_config : list):
    array_ph = nm.Euler_method(alpha, Point, integ_config) #Array dos pontos tal que h > 0
    array_mh = nm.Euler_method(alpha, Point, [-integ_config[0], integ_config[1]]) #Array dos pontos tal que h < 0

    #Retirando o ponto inicial de array_mp
    array_mh = array_mh[~np.all(array_mh == Point, axis = 1)]

    #Concatenando os arrays na ordem: h > 0 e h < 0
    return np.concatenate((array_ph, array_mh))

#Definindo um filtro para armazenar os pontos que pertencam ao dominio do prisma:
def Points_Filter(array, bar, N):
    #Definindo uma lista para armazenar os pontos ditos corretos
    Array_pc = []

    if(bar): #bar tem relacao com a funcao baricentrica 
        for i in range(2*N+1):
            if(if_PointInEq(array[i])):
                Array_pc.append(array[i])
    else:
        for i in range(2*N+1):
            if(if_PointInRet(array[i])):
                Array_pc.append(array[i])
    
    return np.array(Array_pc, float)

#####################################################################
#Defininco funcao que retorna a cor do ponto
def colorPoint(alpha, point):
        lambdaS_value = fun.lbdas(*point)
        lambdaF_value = fun.lbdaf(*point)
        lambdaZ_value = lambdz(*point, alpha)
        if(lambdaZ_value < lambdaS_value < lambdaF_value):
            return ['b', 'g', 'r']
        elif(lambdaS_value < lambdaZ_value < lambdaF_value):
            return ['g', 'b', 'r']
        elif(lambdaS_value < lambdaF_value < lambdaZ_value):
            return ['r', 'b', 'g']
        else:
            return ['k', 'k', 'k']
        
#Funcao que plota um cone - substtuindo a necessidade de plotar um vetor
def plotCone(ax, origin, direction, color, height=0.05, radius=0.015, opacity=1.0, n_theta=12):
    origin    = np.array(origin,    dtype=float)
    direction = np.array(direction, dtype=float)

    norm = np.linalg.norm(direction)
    if norm == 0:
        return
    d = direction / norm

    perp = np.array([1, 0, 0]) if abs(d[0]) < 0.9 else np.array([0, 1, 0])
    u = np.cross(d, perp);  u /= np.linalg.norm(u)
    v = np.cross(d, u);     v /= np.linalg.norm(v)

    apex        = origin + d * (height / 2)
    base_center = origin - d * (height / 2)

    theta = np.linspace(0, 2 * np.pi, n_theta)
    ring  = base_center[:, None] + radius * (np.outer(u, np.cos(theta))
                                           + np.outer(v, np.sin(theta)))

    # Superfície lateral + tampa em um único plot_surface
    # Empilha: [apex_row, ring] → plot_surface interpola entre eles
    apex_row = np.tile(apex[:, None], (1, n_theta))  # repete apex n_theta vezes
    X = np.vstack([apex_row[0], ring[0]])
    Y = np.vstack([apex_row[1], ring[1]])
    Z = np.vstack([apex_row[2], ring[2]])

    ax.plot_surface(X, Y, Z, color=color, alpha=opacity, linewidth=0,
                    antialiased=False, shade=True, zorder = 4)  # antialiased=False é mais rápido


#Funcao que seleciona os pontos de forma simetrica 
def stride_sample_symmetric(points, delta):
    """
    Amostra `points` simetricamente a partir do ponto de maior z,
    percorrendo para frente (+delta) e para tras (-delta).
    Retorna lista ordenada (de trás para frente).
    """
    if len(points) == 0:
        return []

    # Acha o índice do ponto com maior z (índice 2)
    peak_idx = max(range(len(points)), key=lambda i: points[i][2])

    # Percorre para trás (peak → 0)
    backward = list(range(peak_idx, -1, -delta))
    # Percorre para frente (peak → fim)
    forward  = list(range(peak_idx, len(points), delta))

    # Une, remove duplicata do peak e ordena
    indices = sorted(set(backward + forward))
    return [points[i] for i in indices]

#Funcao que remove pontos duplicados consecutivos de branches, ou seja, retira a conexao que ha entre duas branches
def remove_consecutive_duplicates(points, tol=1e-9):
    unique = [points[0]]
    for p in points[1:]:
        if not np.allclose(p, unique[-1], atol=tol):
            unique.append(p)
    return unique

################################################################
#Funcao que le os pontos iniciais de um arquivo txt:
def read_points():
    """Lê pontos de um arquivo txt, um ponto por linha: x y z"""
    points = []
    with open("Data/Initial_Points.txt", 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):  # ignora vazios e comentários
                continue
            coords = list(map(float, line.split()))
            points.append(coords)
    return points

#################################################################
def hugoniotSystemSolver(
                        initialValue : float, #Valor inicial do intervalo da grade
                        finalValue : float,   #Valor final do intervalo da grande
                        Resol : int,          #Resolucao da grade
                        enableMask : bool,    #Variavel que habilita o filtro para o triangulo
                        u0: float,            #Componente u do ponto fixo
                        v0: float,            #Componente v do ponto fixo
                        z0: float,            #Componente z do ponto fixo
                        z: float,             #Plano z constante
                        alpha: float,         #Variavel de controle
                        TOL = 1e-6,           #Tolerancia
                        TOL_residual = 1e-8   #Tolerancia 
                        ) -> list:
    
    """
    [Explicacao] - Dependencias para o codigo funcionar sem aparecer
                   o problema de "chamada circular"
    """
    import includes.Campo_Hugoniot as ch
    
    from scipy.optimize import fsolve

    def hugoniotSystem(vars, u0, v0, z0, z, alpha):
        u_s, v_s = vars #Extrai as variaveis para ser entradas do sistema
        return [
            ch.F(u0, v0, z0, u_s, v_s, z),         #Funcao F
            ch.G(alpha, u0, v0, z0, u_s, v_s, z)   #Funcao G        
        ]

    """
    [Explicacao] - Esta eh uma funcao que soluciona o sistema de Hugoniot acima e, como resultado,
                   retorna uma lista contendo as solucoes do sistema, os valores de F, G na malha
                   e a malha UxV.
    """

    #Grade de valores
    u_vals = np.linspace(initialValue, finalValue, Resol)   #Cria a reta u com 200 pontos
    v_vals = np.linspace(initialValue, finalValue, Resol)   #Cria a reta v com 200 pontos
    U, V   = np.meshgrid(u_vals, v_vals)      #Mescla u e v e retorna uma tupla de valores

    #Mascara regioes invalidas (u + v > 1), i.e, armazena a regiao que nao eh de interesse nosso
    mask = (U + V > 1)

    ZF = ch.F(u0, v0, z0, U, V, z).astype(float)           #Cria uma 'lista' de pontos que foram varridos dados U e V acima
    ZG = ch.G(alpha, u0, v0, z0, U, V, z).astype(float)    #Cria uma 'lista' de pontos que foram varridos dados U e V acima

    #filtra da grade os pontos (u, v) onde simultaneamente |F| ≈ 0 e |G| ≈ 0
    """
    [Explicacao] -> Esses pontos jah estao perto da solucao real,
                    o que aumenta a chance de convergência do fsolve
                    (que serah utilizada posteriormente).
    """

    proximity = (np.abs(ZF) < TOL) & (np.abs(ZG) < TOL) & (~mask if enableMask else True) #Calcula os indices onde ZF e ZG sao proximos de 0

    u_guesses = U[proximity]    #Variaveis de chute (guesses)
    v_guesses = V[proximity]    #Variaveis de chute (guesses)

    #Lista para armazenar as possiveis solucoes
    candidatas = []

    for u_g, v_g in zip(u_guesses, v_guesses):
        try:
            """
            [Explicacao] -> o fsolve eh projetado para encontrar raizes
                            de equacoes ou sistemas de equacoes nao li-
                            neares.
                        ->  u_g e v_g sao variaveis de chute (guess) jah
                            filtradas anteriormente.
            """
            sol = fsolve(hugoniotSystem, [u_g, v_g], args=(u0, v0, z0, z, alpha), full_output=True)
            u_s, v_s = sol[0]                                                            #Possivel solucao do sistema

            residual = np.linalg.norm(hugoniotSystem([u_s, v_s], u0, v0, z0, z, alpha))  #Calcula a norma dos valores obtidos do sistema, isto eh, a diferenca F-G dado u_s e v_s
            inside  = (u_s > 1e-9) and (v_s > 1e-9) and (u_s + v_s < 1)                        #Verifica se u_s e v_s estao dentro do dominio desejado
            not_dup = all(
                np.linalg.norm([u_s - us, v_s - vs]) > TOL for us, vs, zs in candidatas  #Verifica se u_s e v_s nao eh duplicata de uma solucao jah encontrada
                )
            
            if ((residual < TOL_residual) and (inside) and (not_dup)):
                candidatas.append((u_s, v_s, z))

        except:
            continue

    #__________LOG__________#
    print("--------------------------------------------------")
    print(f"Número de soluções encontradas: {len(candidatas)}")
    print("--------------------------------------------------")
    print("########################################################################")
    for i in range(len(candidatas)):
        #Print das solucoes aproximadas:
        print(f"Solucao candidata {i+1}: {candidatas[i]}")
    print("########################################################################")

    return [U, V, ZF, ZG, candidatas]

def sigma(
        u  : float,
        v  : float,
        z  : float,
        u0 : float,
        v0 : float,
        z0 : float,
        tol = 1e-2
):
    f = fun.fw(u, v, z)
    g = fun.fo(u, v, z)

    f0 = fun.fw(u0, v0, z0)
    g0 = fun.fo(u0, v0, z0)

    if(abs(u - u0) < tol):
        return (g - g0)/(v - v0)

    return (f - f0)/(u - u0)

def HugoniotPlanaImplicita(
        alpha : float, #Variavel de controle
        u0    : float, #Componente u do ponto fixo
        v0    : float, #Componente v do ponto fixo
        z0    : float, #Componente z do ponto fixo
        TOL   : float  #Tolerancia para a norma do Maximo
):  
    from matplotlib.pyplot import contour 
    
    def H(u,v):
        return (fun.fw(u,v,z0) - fun.fw(u0,v0,z0))*(v-v0) - (fun.fo(u,v,z0) - fun.fo(u0,v0,z0))*(u-u0)

    # def sig(u, v):
    #     FT = ((fun.fw(u,v,z0) * (v - v0))/(u0 + alpha * az(z0)) - (fun.fo(u,v,z0) - fun.fo(u0,v0,z0)))*(u - u0)
    #     ST = H(u, v)
    #     return FT - ST

    #Define um grid no quadrado [0, 1] X [0, 1]
    x = np.linspace(0.0, 1.0, 750)
    y = np.linspace(0.0, 1.0, 750)
    X, Y = np.meshgrid(x, y)

    Hl = H(X, Y)                    #Hl armazena os pontos em todo o grid [0, 1] X [0, 1]
    Hl = np.where(X + Y < 1, Hl, np.nan)  #Agora Hl armazena os pontos no triangulo X + Y < 1

    #Guarde a referencia da curva na variavel 'cs'
    cs = contour(X, Y, Hl, [0], alpha = 0)

    """
    =====================================================
    PASSO 1: separa cada Path em polilinhas desconectadas
    =====================================================
    [Explicacao] - cs.get_path retorna um segmento contendo varios pontos.
                   O que pode acontecer eh: um "segmento" contendo pontos
                   nao conectados, i.e, pontos que estao distantes, mas o
                   matplotlib os concatenam. Por exemplo, pode haver pon-
                   tos proximos a regiao u = 0 e pontos proximos a u = 1,
                   visualmente eh notorio que sao ramos distintos, mas, 
                   para o matplotlib, sao um unico segmento.
    """
    polilinhas = []
    for path in cs.get_paths():
        path.should_simplify = False
        for poly in path.to_polygons(closed_only = False): #Isto separa em poligonos
            if len(poly) > 1:                              #Se houver mais de um ponto, apeenda na polilinhas
                polilinhas.append(poly)

    """
    =====================================================
    PASSO 2: dentro de cada polilinha, separa em ramos
    nos pontos onde a curva toca a borda do triangulo
    (u=0, v=0, u+v=1)
    =====================================================
    [Explicacao] - Com os segmentos separados, podemos separar as branhes
                   usando os limites do triangulo (no caso, o retangulo).
    """

    dx = x[1] - x[0]     #dx varia dependendo da resolucao definida na linha 362
    tol_limite = 2 * dx  #tolerancia ~ 2 celulas de grid em volta do triangulo

    def bordaTriangulo(coordenada):
        u, v = coordenada
        return (u <= tol_limite) or (v <= tol_limite) or (u + v >= 1 - tol_limite)

    #Lista para armazenar as branches
    branches = []
    for poly in polilinhas:                             #para cada conjunto de pontos
        mask = np.array([bordaTriangulo(p) for p in poly]) #o varremos para encontrar os pontos-limites

        #Caso incomum (Professor, isso eh possivel?)
        #Se a polilinha nunca toca a borda, ela e um ramo fechado (loop) inteiro DENTRO do triangulo
        if(not mask.any()):
            branches.append(poly)
            continue

        #Caso normal: a polilinha vai de uma borda a outra sem tocar
        #a borda no meio -> e um unico ramo (comportamento tipico do contour)
        indices_borda = np.where(mask)[0] #Extrai os indices dos pontos na polilinha onde a condicao eh verdadeira

        #So ha necessidade real de split se houver toques de borda
        #nao-consecutivos no MEIO da polilinha (raro, mas possivel
        #perto de pontos de coincidencia/singularidade da H=0)
        cortes = [0]                            #Armazena o primeiro indice de corte (indice 0)
        for i in range(1, len(indices_borda)):
            if(indices_borda[i] - indices_borda[i-1] > 1): #ha ao menos 3 pontos que compoe o ramo
                #ha um trecho interior entre dois toques de borda:
                #fecha o ramo atual e comeca um novo a partir daqui
                cortes.append(indices_borda[i]) #Armazena os INDICES onde irao ocorrer a separacao em ramos
        cortes.append(len(poly))                #Armazena o ultimo indice de corte (indice len(poly))

        """
        [Explicacao] - Esse trecho (abaixo) define um intervalo
                       [a, b] (chamado trecho) e verifica se no
                       intervalo ha 2 pontos ou mais. Claro, ha
                       um truque do python em

                             zip(cortes[:-1], cortes[1:])

                       Para pegar pares de elementos consecuti-
                       vos de uma lista.
        """
        for a, b in zip(cortes[:-1], cortes[1:]): 
            trecho = poly[a:b]
            if len(trecho) > 10:
                branches.append(trecho)

    #=====================================================
    #PASSO 3: em cada ramo, encontra o(s) ponto(s) onde sig(u,v) = 0
    #via deteccao de mudanca de sinal + interpolacao/bissecao
    #=====================================================
    from scipy.optimize import brentq


    """
    [Explicacao] - Usando dois PONTOS consecutivos, pode-se 
                   interpolar usando a parametrizacao.
                   Com isto, podemos pegar um ponto qualquer
                   entre os dois pontos.
    """
    def sigSegmento(t, p_a, p_b):
        #t em [0,1] parametriza o segmento entre dois vertices consecutivos do ramo
        u = p_a[0] + t * (p_b[0] - p_a[0])
        v = p_a[1] + t * (p_b[1] - p_a[1])
        return H(u, v)

    resultados = []  #lista de (idx_do_ramo, [pontos onde sig(u,v) = 0])
    for idx_ramo, ramo in enumerate(branches):
        print(f"ramo {idx_ramo}: {len(ramo)} pontos")
        sig_vals = np.array([H(p[0], p[1]) for p in ramo]) #Armazena o valor dos sigmas de CADA PONTO do ramo 

        raizes_do_ramo = []                               #Armazena os pontos onde sig(u, v) = 0 na forma [idx, [pontos]]
        for i in range(len(ramo) - 1):
            s0, s1 = sig_vals[i], sig_vals[i+1]           #Extrai o valor do sigma do i-esimo ponto e seu sucessor

            if np.isnan(s0) or np.isnan(s1):                    #Se ui ou vi sao nan
                continue                                  #pula a atual iteracao

            if s0 == 0.0:                                 #Se a extremidade inferior tiver sig = 0
                raizes_do_ramo.append(tuple(ramo[i]))     #Armazena O PONTO.
                continue

            if s0 * s1 < 0:
                #ha mudanca de sinal nesse segmento -> raiz aqui dentro

                """
                [Explicacao] - o brenqt utiliza metodos para encontrar
                               raizes de funcao. Neste caso, eh notorio
                               que sigSegmento depende exclusivamente de
                               t, entao brentq ACHA, dado uma TOL, qual
                               valor de t para zerar o sig entre dois 
                               pontos.
                
                [Explicacao] - Portanto, como sabemos o valor de t que re-
                               sulta em sig = 0, eh so o utilizar para en-
                               contrar u e v que satisfacam sig(u,v) = 0.
                """
                t_raiz = brentq(
                    sigSegmento, 0.0, 1.0,
                    args=(ramo[i], ramo[i+1]), xtol=TOL
                )

                u_raiz = ramo[i][0] + t_raiz * (ramo[i+1][0] - ramo[i][0])
                v_raiz = ramo[i][1] + t_raiz * (ramo[i+1][1] - ramo[i][1])
                raizes_do_ramo.append((u_raiz, v_raiz))

        resultados.append((idx_ramo, raizes_do_ramo))

    return branches, resultados