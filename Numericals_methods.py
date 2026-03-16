"""
Modulo: Integracao pelo metodo de Euler

Objetivo:
- Utilizar o metodo de Euler para armazenar o conjunto de pontos que informam o comportamento
dos compoentes do campo a partir de um ponto inicial escolhido.
"""

import Campo_Ez as df
from Functions import lbdc
from Inicia import baricentrica
from numpy import array, ceil, abs, log10, sqrt

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
def Bisection_method(interval : list, bisection_config : list, function, z, e = 1e-10): 
    a0, b0 = interval       #Extraindo intervalo
    h, N = bisection_config #Extraindo as configuracoes para realizar o metodo

    #Definindo o numero de iteracoes:
    Num_i = int(ceil((log10(b0 - a0) - log10(e))/log10(2)))

    #Definindo o valor inicial:
    vk = 0         #Reta horizontal v = constante (como eh o valor inicial, entao constante = 0)

    #Funcao que calcula a diferenca entre as funcoes
    def d(x, y, z): #x, y, z sao o uk, vk e z 
        return function(x, y, z) - lbdc(x, y, z)
    
    Points = []    #Variavel para armazenar os pontos
    altura = 0     #Constante para controle
    if(baricentrica):
        altura = sqrt(3)/2
    else:
        altura = 1

    for _ in range(N):
        a, b = a0, b0

        if(vk >= altura):                  #Quebra iteracao se a linha v = k estiver fora do triangulo
            break

        if(d(a, vk, z) * d(b, vk, z) > 0): #Pula iteracao se nao houver mudança de sinal no intervalo inicial
            vk += h
            continue

        for _ in range(Num_i):
            uk = (a+b)/2                   #Calculo o ponto medio do intervalo

            if(abs(uk - a) <= e * max(abs(uk), 1)):
                Points.append([uk, vk, z]) #Armazena o ponto em que a diferenca eh proximo a ZERO
                break
            else:
                if(d(a, vk, z) * d(uk, vk, z) < 0): 
                    b = uk
                elif(d(a, vk, z) * d(uk, vk, z) >= 0):
                    a = uk
        
        vk += h                            #Isso direciona o Metodo da bissecao para a proxima reta horizontal

    return Points
