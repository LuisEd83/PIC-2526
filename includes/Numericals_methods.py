"""
Modulo: Integracao pelo metodo de Euler

Objetivo:
- Utilizar o metodo de Euler para armazenar o conjunto de pontos que informam o comportamento
dos compoentes do campo a partir de um ponto inicial escolhido.
"""

import includes.Campo_Ez as df
import includes.Campo_Hugoniot as ch

import includes.Functions as fun
from includes.Inicia import baricentrica

from numpy import array, ceil, abs, log10, sqrt, linspace

def az(z): #Da/dz
    return 1/(1+z**2)

def lambdz(u, v, z, alpha):
    return(fun.fw(u, v, z)/(u + alpha*az(z)))

#Definindo uma funcao que implementa o metodo de Euler para integracao
def Euler_method(alpha, point : list, integ_config : list):
    #Extraindo ponto inicial
    u0, v0, z0 = point

    #Extraindo configuracao de integracao
    h, N = integ_config
    
    #Inicializando variaveis para o metodo de Euler
    uk = u0
    vk = v0
    zk = z0

    #Inicializando uma lista de pontos:
    Points = [[u0, v0, z0]]

    #Calculando e armazenando os resultados do metodo de Euler:
    for _ in range(N): #O laco vai repetir N vezes
        #Realizo o passo (obs: kp1 = k + 1)
        ukp1 = uk + h * df.P(uk, vk, zk, alpha)
        vkp1 = vk + h * df.Q(uk, vk, zk, alpha)
        zkp1 = zk + h * df.R(uk, vk, zk, alpha)

        #Armazenando os valores na lista de pontos
        Points.append([ukp1, vkp1, zkp1])

        #Atualizo os valores das variaives uk, vk e zk 
        uk = ukp1
        vk = vkp1
        zk = zkp1

    #Definindo um array para armazenar pontos (com tamanho N + 1):
    Points = array(Points, float)

    return Points

#Definindo uma funcao que implementa o Metodo da Bissecao adaptado ao problema
#   -> Interval: Eh o intervalo [a, b]
#   -> bisection_config: Eh o passo e o numero de passos
#   -> function: eh a funcao lambdaS ou lambdaF
#   -> z: Eh a curva de nivel
#   -> tol: Eh a tolerancia
def Bisection_method(interval : list, bisection_config : list, function, z, alpha, e = 1e-10): 
    a0, b0 = interval       #Extraindo intervalo

    h, N = bisection_config #Extraindo as configuracoes para realizar o metodo

    if(h > 0.1):
        h = 0.01

    #Definindo o numero de iteracoes:
    Num_i = int(ceil((log10(b0 - a0) - log10(e))/log10(2)))

    #Definindo o valor inicial:
    vk = 0         #Reta horizontal v = constante (como eh o valor inicial, entao constante = 0)

    #Funcao que calcula a diferenca entre as funcoes
    def d(x, y, z): #x, y, z sao o uk, vk e z 
        if(alpha != 0):
            return function(x, y, z) - lambdz(x, y, z, alpha)
        return function(x, y, z) - fun.lbdc(x, y, z)
    
    Points = []                               #Variavel para armazenar os pontos
    altura = sqrt(3)/2 if baricentrica else 1 #Constante para controle de altura

    for _ in range(N):
        if(vk >= altura):                  #Quebra iteracao se a linha v = k estiver fora do triangulo
            break
        
        a, b = a0, b0                      #Armazena o intervalo original

        #Armazenando os valores da funcao d(x, y, z)
        da = d(a, vk, z)
        db = d(b, vk, z)

        if(abs(da) <= e):       #Se d(a, vk, z) for proximo de zero, isso significa que d(a, vk, z) eh uma zero espurio 
            a += h              #Adiciona um passo ao "a" 
            da = d(a, vk, z)    #Atualiza o da

        if(da*db > 0):          #Se nao houver mudanca de sinal no intervalo [a, b] para v = k, entao pulamos para o v = k + h
            vk += h
            continue

        for _ in range(Num_i):
            uk = (a+b)/2                      #Calculo o ponto medio do intervalo

            if(abs(uk - a) <= e * max(abs(uk), 1)):
                Points.append([uk, vk, z])    #Armazena o ponto em que a diferenca eh proximo a ZERO
                break
            else:
                product = da * d(uk, vk, z)

                if(product < 0):               #Se ha a troca de sinais
                    b = uk
                elif(product > 0):             #Se nao ha a troca de sinal
                    a = uk
                    da = d(a, vk, z)           #Atualiza o da para a proxima iteracao
                else:
                    Points.append([uk, vk, z]) #Raiz exata
                    break
        
        vk += h                                #Isso direciona o Metodo da bissecao para a proxima reta horizontal

    return Points


#Relacionado ao campo de vetores de Hugonoit
#   ->alpha:        variavel de controle
#   ->guess_point:  Ponto inicial de chute
#   ->Point:        Ponto inicial para a integração
#   ->integ_config: configuração para integração
#=> Retorno: um array numpy de pontos
def RungeKutta4(
        alpha              : float,
        fixed_point        : list,
        Point              : list,
        integ_config       : list
) -> array:
    #Extraindo ponto de chute
    u0, v0, z0 = fixed_point

    #Extraindo ponto inicial 
    u1, v1, z1 = Point

    #Extraindo configuracao de integracao
    h, N = integ_config
    
    #Inicializando variaveis para o metodo de Runge-Kutta
    uk = u1
    vk = v1
    zk = z1

    #Inicializando uma lista de pontos:
    Points = [[u1, v1, z1]]

    for _ in range(N):

        #Relacionado a variavel u:
        k1 = ch.Hug1(alpha, u0, v0, z0, uk, vk, zk)
        k2 = ch.Hug1(alpha, u0, v0, z0, uk + (h*k1)/2, vk + (h*k1)/2, zk + (h*k1)/2)
        k3 = ch.Hug1(alpha, u0, v0, z0, uk + (h*k2)/2, vk + (h*k2)/2, zk + (h*k2)/2)
        k4 = ch.Hug1(alpha, u0, v0, z0, uk + (h*k3), vk + (h*k3), zk + (h*k3))

        uk = uk + (h/6) * (k1 + 2*k2 + 2*k3 + k4) #Atualizo o valor de uk

        #Relacionado a variavel v:
        k1 = ch.Hug2(alpha, u0, v0, z0, uk, vk, zk)
        k2 = ch.Hug2(alpha, u0, v0, z0, uk + (h*k1)/2, vk + (h*k1)/2, zk + (h*k1)/2)
        k3 = ch.Hug2(alpha, u0, v0, z0, uk + (h*k2)/2, vk + (h*k2)/2, zk + (h*k2)/2)
        k4 = ch.Hug2(alpha, u0, v0, z0, uk + (h*k3), vk + (h*k3), zk + (h*k3))

        vk = vk + (h/6) * (k1 + 2*k2 + 2*k3 + k4) #Atualizo o valor de vk

        #Relacionado a variavel z:
        k1 = ch.Hug3(alpha, u0, v0, z0, uk, vk, zk)
        k2 = ch.Hug3(alpha, u0, v0, z0, uk + (h*k1)/2, vk + (h*k1)/2, zk + (h*k1)/2)
        k3 = ch.Hug3(alpha, u0, v0, z0, uk + (h*k2)/2, vk + (h*k2)/2, zk + (h*k2)/2)
        k4 = ch.Hug3(alpha, u0, v0, z0, uk + (h*k3), vk + (h*k3), zk + (h*k3))

        zk = zk + (h/6) * (k1 + 2*k2 + 2*k3 + k4) #Atualizo o valor de zk

        #Por conta da construcao do problema:
        if(abs(zk - z0) < 2e-1): #zk nao deve ficar proximo de z0
            break

        Points.append([uk, vk, zk])               #Armazeno o ponto

    #Definindo um array para armazenar pontos (com tamanho N + 1):
    Points = array(Points, float)

    return Points

def runge_Kutta_Scipy(
        alpha              : float,
        guess_point        : list,
        Point              : list,
        integ_config       : dict
):
    import scipy.integrate as solve_ip

    #Funcao que serah utilizada no runge-kutta 5 (do scipy)
    def function(s, state : list):
        u, v, z = state

        H = ch.CampoHug(alpha, guess_point[0], guess_point[1], guess_point[2], u, v, z) 

        return[H[0], H[1], H[2]]

    #Configuracao para integracao
    s_inicial = integ_config['s_inicial']
    s_final   = integ_config['s_final']
    n_pontos  = integ_config['n_pontos']


    s_span = (s_inicial, s_final)                   #Intervalo de integracao
    s_eval = linspace(s_inicial, s_final, n_pontos) #Particao do intervalo de integracao

    solucao = solve_ip.solve_ivp(
        function,
        s_span,
        Point,
        method = 'RK45',
        t_eval = s_eval
    )

    #solucao.y tem o formato (3, N) -> 3 variáveis ao longo de N pontos
    #solucao.y.T inverte para (N, 3) -> N pontos, cada um com 3 coordenadas [u, v, z]
    Points = array(solucao.y.T, dtype = float)

    return Points

def HugonioutEuler_method(alpha, fixed_point : list, P_chute : list, integ_config : list, TOL = 5e-3):
    #Extraindo ponto fixado
    u0, v0, z0 = fixed_point 

    #Extraindo ponto inicial para integracao
    u1, v1, z1 = P_chute #U1

    #Extraindo configuracao de integracao
    h, N = integ_config
    
    #Inicializando variaveis para o metodo de Euler
    uk = u1
    vk = v1
    zk = z1

    #Inicializando uma lista de pontos:
    Points = [[u1, v1, z1]]

    #Calculando e armazenando os resultados do metodo de Euler:
    for _ in range(N): #O laco vai repetir N vezes
        #Realizo o passo (obs: kp1 = k + 1)
        #(alpha, u0, v0, z0, u, v, z):
        ukp1 = uk + h * ch.Hug1(alpha, u0, v0, z0, uk, vk, zk)
        vkp1 = vk + h * ch.Hug2(alpha, u0, v0, z0, uk, vk, zk)
        zkp1 = zk + h * ch.Hug3(alpha, u0, v0, z0, uk, vk, zk)

        if((abs(zkp1 - z0) <= TOL)):
            break

        # if(ch.Norma(alpha, u0, v0, z0, ukp1, vkp1, zkp1) < TOL):
        #     break

        #Armazenando os valores na lista de pontos
        Points.append([ukp1, vkp1, zkp1])

        #Atualizo os valores das variaives uk, vk e zk 
        uk = ukp1
        vk = vkp1
        zk = zkp1

    #Definindo um array para armazenar pontos (com tamanho N + 1):
    Points = array(Points, float)
    return Points