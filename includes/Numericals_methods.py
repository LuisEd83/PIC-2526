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

from numpy import array, ceil, abs, log10, sqrt, linspace, exp

k = 0.01
def a(z):
    return 1/(1 + exp(-k*z)) - 0.5

def az(z):
    return k*a(z)*(1-a(z))

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
    altura = 1 #Constante para controle de altura

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

def smoller_Point(
                alpha : float, #Variavel de controle
                u0    : float, #Componente u do ponto fixo
                v0    : float, #Componente v do ponto fixo
                z0    : float, #Componente z do ponto fixo
                sig0  : float, #Componente sig do ponto fixo
                h     : float, #Passo
                TOL   : float  #Tolerancia para a norma do Maximo
                ) -> list:

    import includes.Campo_Ez as ce

    #Cria uma copia do passo
    h0 = h

    #Variavel para o criterio de parada 
    MAX_ITERATOR = 1000
    it = 0 

    while(True):
        uk = u0 + h0 * ce.P(u0, v0, z0, alpha)
        vk = v0 + h0 * ce.Q(u0, v0, z0, alpha)
        zk = z0 + h0 * ce.R(u0, v0, z0, alpha)
        sigk = lambdz(u0, v0, z0, alpha) + h0 * 0.5 #o 0.5 veio do Smoller

        #Norma do maximo
        #norm = abs(zk - z0) <= FUturamente
        norm = max(abs(uk - u0), abs(vk - v0), abs(zk - z0), abs(sigk - sig0))

        if(norm < TOL):
            h0 += h
            it += 1
            if(it == MAX_ITERATOR):
                print(f"[ERROR] - Ponto singular?")
                print(f"[ERROR] - Ponto fixo: {u0, v0, z0}")
                print(f"[ERROR] - Ponto calculado apos {it} iteracoes: {uk, vk, zk}")
                print(f"[ERROR] - Valor da norma do Maximo: {norm}")
                exit()

        else:
            print(f"[LOG] - Ponto calculado apos {it} iteracoes: {uk, vk, zk}")
            print(f"[LOG] - Valor da norma do Maximo: {norm}")
            print(f"[LOG] - valor do h0 : {h0} ")
            break
    
    return uk, vk, zk, sigk

def HugoniotNewtonRaphson_method(
        alpha           : float,
        current_point   : list,
        fixed_point     : list,
        tol      = 1e-11,
        max_iter = 100,
        gama     = 0.25             #Variavel de controle do passo do metodo de Newton-Raphson
):

    """
    [IDEIA] - Utilizar as funções F e G em conjunto com o método de Newton-Raphon para
              funções multivariáveis para localizar, no plano z = k, a sua raíz e cor-
              rigir, assim, a curva Hugoniot para cada ponto.

              Para o jacobiano que será necessário, já possuímos as derivadas parciais
              para cada variável (u, v e z), como estamos trabalhando em um plano z = k,
              a matriz jacobiana será 2x2, tornando o trabalho mais "fácil".
    """

    #Extraindo as componentes:
    #U0
    u0, v0, z0 = fixed_point

    #U1
    u1, v1, z1 = current_point

    def deltaUV(u, v):
        #Avaliacao das funcoes
        F = ch.F(u0, v0, z0, u, v, z1)
        G = ch.G(alpha, u0, v0, z0, u, v, z1)

        #Avaliacao das derivadas parciais(Jacobiano 2x2)
        Fu = ch.Fu(u0, v0, z0, u, v, z1)
        Fv = ch.Fv(u0, v0, z0, u, v, z1)
        Gu = ch.Gu(alpha, u0, v0, z0, u, v, z1)
        Gv = ch.Gv(alpha, u0, z0, u, v, z1)

        #Determinante do Jacobiano
        dJ = Fu * Gv - Fv * Gu
        # print(f"Valor do det|J(u,v)| = {dJ}")

        if(abs(dJ) < 1e-12):
            raise ValueError("Jacobiano singular ou muito próximo de zero.")

        #Solucao exata da inversao da matriz 2x2: J * delta = -F
        du = gama * ((Fv * G - Gv * F) / dJ)
        dv = gama * ((Gu * F - Fu * G) / dJ)

        return du, dv, F, G, dJ

    #Definindo o ponto de partida no plano z = z1
    uk, vk = u1, v1

    for it in range(max_iter):
        du, dv, F, G, dJ = deltaUV(uk, vk)
        #Critério de parada: norma dos residuos |F| e |G| abaixo da tolerancia
        if((abs(F) < tol) and (abs(G) < tol)):
            return uk, vk
        
        #Atualizacao dos pontos
        uk += du
        vk += dv


    raise RuntimeError("[ERROR] - O metodo de Newton-Raphson nao convergiu dentro do numero maximo de iteracoes")

def HugonioutEuler_method(
            alpha        : float,   #Variavel de controle   
            fixed_point  : list,    #Ponto fixo
            integ_config : list,    #Configuracao para integracao
            U4 = None,              #Ponto 
            flag = True,            #Flag para ativar ou nao o smoller
            tol = 1e-7,             #Tolerancia default
):
    #Extraindo ponto fixo:
    u0, v0, z0 = fixed_point
    sig0 = lambdz(u0, v0, z0, alpha) #sigma inicial

    #Extraindo configuracao da integracao
    h, N = integ_config

    #iniciando a lista de pontos
    Points = [[u0, v0, z0, sig0]] if(flag) else []

    #Iniciando os pontos iniciais para integracao
    uSmoller, vSmoller, zSmoller, sigSmoller = smoller_Point(alpha, u0, v0, z0, sig0, h, tol)

    uk, vk, zk, sigk = (uSmoller, vSmoller, zSmoller, sigSmoller) if(flag) else U4
    Points.append([uk, vk, zk, sigk])

    #nova tolerancia:
    zTol = abs(zSmoller - z0)

    if(ch.R_sig(alpha, uk, vk, zk, sigk, u0, v0, z0) > 0):
        if(((zk > z0) and (h < 0)) or ((zk < z0) and (h > 0))):
            h = -h
    else:
        if(((zk > z0) and (h > 0)) or ((zk < z0) and (h < 0))):
            h = -h

    for _ in range(N):
        if(abs(ch.F(u0, v0, z0, uk, vk, zk)) > tol or abs(ch.G(alpha, u0, v0, z0, uk, vk, zk)) > tol):
            try:
                uk, vk = HugoniotNewtonRaphson_method(alpha, [uk, vk, zk], fixed_point)
                sigk = ch.sigm(alpha, uk, vk, zk, u0, v0, z0)
            except Exception as e:
                print(f"[LOG] - {e}")

        ukp1    = uk   +  h * ch.P_sig(alpha, uk, vk, zk, sigk, u0, v0, z0)
        vkp1    = vk   +  h * ch.Q_sig(alpha, uk, vk, zk, sigk, u0, v0, z0)
        zkp1    = zk   +  h * ch.R_sig(alpha, uk, vk, zk, sigk, u0, v0, z0)
        sigkp1  = sigk +  h * ch.S(alpha, uk, vk, zk, sigk, u0, v0, z0)

        #Appendo o ponto
        Points.append([ukp1, vkp1, zkp1, sigkp1])

        #Se o ponto atual estiver proximo do plano z = z0 ou z = 0
        if((abs(zkp1 - z0) < zTol) or (abs(zkp1) < tol)):
            break;
        
        #Atualizo os componentes usando Euler
        uk, vk, zk, sigk = ukp1, vkp1, zkp1, sigkp1
            
    Points = array(Points, float)                       #Transformo em um array numpy
    return [[[p[0], p[1], p[2]], p[3]] for p in Points] #Retorno da forma [[u, v, z], sig]
