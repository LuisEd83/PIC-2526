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
    f = 0 #s = 1 para permitir que seja plotado o ramo
    return f

def points_curv(): #Esta funcao habilita os pontos na curva integral.
    x = 1   #x = 1 para habilitar
    return x

##############################################################
#Derivada da funcao continua e crescente no intervalo [0,1]
def az(z): #Da/dz
    return np.cos(z) #Derivada de a(z) em relacao a z

def lambdz(u, v, z, alpha):
    return(fun.fw(u, v, z)/(u + alpha*az(z)))

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
        if((lambdaZ_value < lambdaF_value < lambdaS_value) or (lambdaZ_value < lambdaS_value < lambdaF_value)):
            return 'b'
        elif((lambdaS_value < lambdaZ_value < lambdaF_value) or (lambdaF_value < lambdaZ_value < lambdaS_value)):
            return 'g' #Adiciona a cor Vermelha
        elif((lambdaF_value < lambdaS_value < lambdaZ_value) or (lambdaS_value < lambdaF_value < lambdaZ_value)):
            return 'r' #Adiciona a cor Roxa
        else:
            return 'k'
        
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
