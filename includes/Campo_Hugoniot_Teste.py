import includes.Functions as fun
import includes.Auxiliar_Functions as af
import includes.Campo_Ez as ce

from numpy import sqrt, array, isnan, linalg

def det_sig(
              alpha  : float,
              u     : float,
              v     : float,
              z     : float,
              sigma  : float
              ):
    """
    [Explicacao] - Este eh o determinante do sistema de equacoes em relacao ao sigma

                 - O calculo consiste em aplicar Laplace. Como na ultima linha da matriz
                 possui dois zeros, entao posso pegar o elemento de indice (3,3) e multi-
                 plicar pelo determinante do seu cofator.
    """

    #Primeiro termo
    firstT = (fun.fwu(u, v, z) - sigma) * (fun.fov(u,v,z) - sigma)

    #Segundo termo
    secondT = fun.fwv(u, v, z) * fun.fou(u, v, z)

    return (fun.fw(u, v, z) - sigma * (u + alpha * af.az(z))) * (firstT - secondT)
    
def det_usig(
              alpha : float,
              u     : float,
              v     : float,
              z     : float,
              sigma : float,
              u0    : float,
              v0    : float,
              z0    : float
              ):

    firstT1   = (v - v0) * (fun.fw(u, v, z) - sigma * (u + alpha * af.az(z)))
    secondT1 = fun.foc(u, v, z) * (alpha * (af.a(z) - af.a(z0)) - (z - z0) * u0)

    firstT2 = (u - u0) * (fun.fw(u, v, z) - sigma * (u + alpha * af.az(z)))
    secondT2 = fun.fwc(u, v, z) * (alpha * (af.a(z) - af.a(z0)) - (z - z0) * u0)

    return - fun.fwv(u, v, z) * (firstT1 - secondT1) + (fun.fov(u, v, z) - sigma) * (firstT2 - secondT2)


def det_vsig(
              alpha : float,
              u     : float,
              v     : float,
              z     : float,
              sigma : float,
              u0    : float,
              v0    : float,
              z0    : float
              ):

    firstT1   = (v - v0) * (fun.fw(u, v, z) - sigma * (u + alpha * af.az(z)))
    secondT1 = fun.foc(u, v, z) * (alpha * (af.a(z) - af.a(z0)) - (z - z0) * u0)

    firstT2 = (u - u0) * (fun.fw(u, v, z) - sigma * (u + alpha * af.az(z)))
    secondT2 = fun.fwc(u, v, z) * (alpha * (af.a(z) - af.a(z0)) - (z - z0) * u0)

    return (fun.fwu(u, v, z) - sigma) * (firstT1 - secondT1) - fun.fou(u, v, z) * (firstT2 - secondT2)

def det_zsig(
              alpha : float,
              u     : float,
              v     : float,
              z     : float,
              sigma : float,
              u0    : float,
              v0    : float,
              z0    : float
              ):
    
    """
    [Explicacao] - Este eh o determinante do sistema de equacoes em relacao [a z?]

                 - O calculo consiste em aplicar Laplace. Como na ultima linha da matriz
                 possui dois zeros, entao posso pegar o elemento de indice (3,3) e multi-
                 plicar pelo determinante do seu cofator.
    """
    #Primeiro termo
    firstT = (fun.fwu(u, v, z) - sigma) * (fun.fov(u,v,z) - sigma)

    #Segundo termo
    secondT = fun.fwv(u, v, z) * fun.fou(u, v, z)

    return (alpha * (af.a(z) - af.a(z0)) - (z - z0) * u0) * (firstT - secondT)

def norm_sig(
              alpha : float,
              u     : float,
              v     : float,
              z     : float,
              sigma : float,
              u0    : float,
              v0    : float,
              z0    : float
              ):
    return sqrt(
        det_sig(alpha, u, v, z, sigma)**2
        + det_usig(alpha, u, v, z, sigma, u0, v0, z0)**2
        + det_vsig(alpha, u, v, z, sigma, u0, v0, z0)**2
        + det_zsig(alpha, u, v, z, sigma, u0, v0, z0)**2
    )

def P_sig(
              alpha : float,
              u     : float,
              v     : float,
              z     : float,
              sigma : float,
              u0    : float,
              v0    : float,
              z0    : float
              ):
    return det_usig(alpha, u, v, z, sigma, u0, v0, z0)/norm_sig(alpha, u, v, z, sigma, u0, v0, z0)

def Q_sig(
            alpha : float,
            u     : float,
            v     : float,
            z     : float,
            sigma : float,
            u0    : float,
            v0    : float,
            z0    : float
        ):
    return det_vsig(alpha, u, v, z, sigma, u0, v0, z0)/norm_sig(alpha, u, v, z, sigma, u0, v0, z0)

def R_sig(
            alpha : float,
            u     : float,
            v     : float,
            z     : float,
            sigma : float,
            u0    : float,
            v0    : float,
            z0    : float
        ):
    return det_zsig(alpha, u, v, z, sigma, u0, v0, z0)/norm_sig(alpha, u, v, z, sigma, u0, v0, z0)

def S(
            alpha : float,
            u     : float,
            v     : float,
            z     : float,
            sigma : float,
            u0    : float,
            v0    : float,
            z0    : float
        ):
      return det_sig(alpha, u, v, z, sigma)/norm_sig(alpha, u, v, z, sigma, u0, v0, z0)

def campo(    
            alpha : float,
            u     : float,
            v     : float,
            z     : float,
            sigma : float,
            u0    : float,
            v0    : float,
            z0    : float
        ):
    return [P_sig(alpha, u, v, z, sigma, u0, v0, z0),
            Q_sig(alpha, u, v, z, sigma, u0, v0, z0),
            R_sig(alpha, u, v, z, sigma, u0, v0, z0),
            S(alpha, u, v, z, sigma, u0, v0, z0)]


def smoller_Point(
                alpha : float,
                u0    : float,
                v0    : float,
                z0    : float,
                h     : float,  
                TOL   : float
                ) -> list:
    
    h0 = h

    MAX_ITERATOR = 40
    it = 0
    while(True):
        uk = u0 + h0 * ce.P(u0, v0, z0, alpha)
        vk = v0 + h0 * ce.Q(u0, v0, z0, alpha)
        zk = z0 + h0 * ce.R(u0, v0, z0, alpha)
        sigk = af.lambdz(u0, v0, z0, alpha) + h0 * 0.5 #o 0.5 veio do Smoller
        norm = norm_sig(alpha, uk, vk, zk, sigk, u0, v0, z0)

        if any(isnan(comp) for comp in [uk, vk, zk, sigk]):
                print("[ERROR] - Ha componentes do tipo nan")
                exit()

        if(norm < TOL ):
            h0 += h
            it += 1
            if(it == MAX_ITERATOR):
                print(f"[ERROR] - Ponto singular?")
                print(f"[ERROR] - Ponto fixo: {u0, v0, z0}")
                print(f"[ERROR] - Ponto calculado apos {it} iteracoes: {uk, vk, zk}")
                print(f"[ERROR] - Valor da norma: {norm}")
                exit()

        else:
            print(f"[DEBUG] - Ponto calculado apos {it} iteracoes: {uk, vk, zk}")
            print(f"[DEBUG] - Valor da norma: {norm}")
            print(f"[DEBUG] - valor do h0 : {h0} ")
            break
    
    return [uk, vk, zk, sigk, h0]

#############################################################################################################
def HugonioutEuler_method(alpha, fixed_point : list, integ_config : list,  TOL = 1e-3):

    #Extraindo configuracao de integracao
    h, N = integ_config

    #Extraindo ponto fixado
    u0, v0, z0 = fixed_point 
    sig0 = af.lambdz(u0, v0, z0, alpha)

    #Inicializando uma lista de pontos:
    Points = [[u0, v0, z0, sig0]]

    #Inicializando variaveis para o metodo de Euler
    uk, vk, zk, sigk, h0 = smoller_Point(alpha, u0, v0, z0, h, TOL)
    _, _, _, _, candidatas = af.hugoniotSystemSolver(0.01, 1, 2000, True, u0, v0, z0, z0 + 200*z0*h0, alpha)

    if(len(candidatas) > 0):
        #Primeiro ponto: mais proximo do fixed_point

        """
        [Explicacao] - lambda eh uma funcao anonima (uma função sem nome), definida em uma linha.
                       min varre a lista de candidatas em busca de minimizar a distancia euclidiana 
        """
        newFixedPoint = min(
            candidatas,
            key=lambda c: linalg.norm([c[0] - u0,
                                          c[1] - v0,
                                          c[2] - z0])
        )

    uk, vk, zk = newFixedPoint

    Points.append([uk, vk, zk, sigk])

    #Calculando e armazenando os resultados do metodo de Euler:
    for i in range(N): #O laco vai repetir N vezes
        #Realizo o passo (obs: kp1 = k + 1)
        ukp1 = uk + h * P_sig(alpha, uk, vk, zk, sigk, *fixed_point)
        vkp1 = vk + h * Q_sig(alpha, uk, vk, zk, sigk, *fixed_point)
        zkp1 = zk + h * R_sig(alpha, uk, vk, zk, sigk, *fixed_point)
        sigkp1 = sigk + h * S(alpha, uk, vk, zk, sigk, *fixed_point)

        #Armazenando os valores na lista de pontos
        #print(f"Ponto {i} : {ukp1, vkp1, zkp1, sigkp1}")
        Points.append([ukp1, vkp1, zkp1, sigkp1])

        #Atualizo os valores das variaives uk, vk e zk 
        uk = ukp1
        vk = vkp1
        zk = zkp1
        sigk = sigkp1

    #Definindo um array para armazenar pontos (com tamanho N + 1):
    Points = array(Points, float)

    return [[[p[0], p[1], p[2]], p[3]] for p in Points] #[[u, v, z], sig]