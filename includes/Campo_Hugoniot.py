import includes.Functions as fun
import includes.Auxiliar_Functions as af
import includes.Campo_Ez as ce

from numpy import sqrt, array, linspace, meshgrid, where, nan, isnan, isclose
from matplotlib.pyplot import contour

#=============================================
# Antiga interpretacao do Hugoniot (3 dimensoes)
#=============================================

def h0(
        z,              #Concentracao variavel entre [0,1]
        z0             #Concentracao inicial
      )-> float:
    
    if(isclose(z, z0)):
        return af.az(z) #Retorna (d/dt)(a(z)) se z = z0
    
    return (af.a(z) - af.a(z0))/(z - z0)

def sig(
        alpha,          #Variavel de controle
        z,              #Variavel pertencente ao intervalo [0,1]
        u0,             #u inicial
        v0,             #v inicial
        z0              #z inicial
        )-> float:      #retorno esperado
    """
    [Explicacao] 
    """
    num = fun.fw(u0, v0, z0)        #Numerador
    den = u0 + alpha * h0(z, z0)    #Denominador

    return num/den

def sigm(
    alpha : float,
    u     : float,
    v     : float,
    z     : float,
    u0    : float,
    v0    : float,
    z0    : float
):
    if(isclose(u, u0)):
        return ((fun.fo(u, v, z0) - fun.fo(u0, v0, z0))/(v - v0))
    elif(isclose(v, v0)):
        return ((fun.fw(u, v, z0) - fun.fw(u0, v0, z0))/(u - v0))
    else:
        return sig(alpha, z, u0, v0, z0)

##############################
def F(
    u0, #Variavel do ponto fixo U0
    v0, #Variavel do ponto fixo U0
    z0, #Variavel do ponto fixo U0
    u,  #Variavel do ponto U1
    v,  #Variavel do ponto U1
    z   #Variavel do ponto U1
):
    #Define-se os termos (para simplicidade):
    firstT  = (fun.fw(u, v, z) - fun.fw(u0, v0, z0)) * (v - v0)
    secondT = (fun.fo(u, v, z) - fun.fo(u0, v0, z0)) * (u - u0)

    return firstT - secondT

def G(
    alpha,  #Varivavel de controle
    u0,     #Varivavel do ponto fixo U0
    v0,     #Varivavel do ponto fixo U0
    z0,     #Varivavel do ponto fixo U0
    u,      #Varivavel do ponto U1
    v,      #Varivavel do ponto U1
    z       #Varivavel do ponto U1
):
    #Define-se os termos (para simplicidade):
    firstT  = (fun.fw(u, v, z) - fun.fw(u0, v0, z0)) * (u0 * (z - z0) + alpha * (af.a(z) - af.a(z0)))
    secondT = fun.fw(u0, v0, z0) * (z - z0) * (u - u0)

    return firstT - secondT

##############################
def Fu(
       u0, #u inicial
       v0, #v inicial
       z0, #z inicial
       u, #Variaveis
       v, #Variaveis
       z  #Variaveis
       )-> float:
    
    """
    [Explicacao] - Esta eh a derivada pacial de F em relacao a variavel u
    """

    #Definindo o primeiro termo:
    firstT  = fun.fwu(u, v, z) * (v - v0)
    #Definindo o segundo termo:
    secondT = fun.fou(u, v, z) * (u - u0)
    #Definindo o terceiro termo:
    thirdT  = fun.fo(u, v, z) - fun.fo(u0, v0, z0)

    return firstT - secondT - thirdT

def Fv(
       u0, #u inicial
       v0, #v inicial
       z0, #z inicial
       u, #Variaveis
       v, #Variaveis
       z  #Variaveis
       )-> float:
    
    """
    [Explicacao] - Esta eh a derivada pacial de F em relacao a variavel v
    """

    #Definindo o primeiro termo:
    firstT  = fun.fwv(u, v, z) * (v - v0)
    #Definindo o segundo termo:
    secondT = fun.fov(u, v, z) * (u - u0)
    #Definindo o terceiro termo:
    thirdT  = fun.fw(u, v, z) - fun.fw(u0, v0, z0)

    return firstT - secondT + thirdT

def Fv(
       u0, #u inicial
       v0, #v inicial
       z0, #z inicial
       u, #Variaveis
       v, #Variaveis
       z  #Variaveis
       )-> float:
    
    """
    [Explicacao] - Esta eh a derivada pacial de F em relacao a variavel v
    """

    #Definindo o primeiro termo:
    firstT  = fun.fwv(u, v, z) * (v - v0)
    #Definindo o segundo termo:
    secondT = fun.fov(u, v, z) * (u - u0)
    #Definindo o terceiro termo:
    thirdT  = fun.fw(u, v, z) - fun.fw(u0, v0, z0)

    return firstT - secondT + thirdT

"""def Fz(
       u0, #u inicial
       v0, #v inicial
       u, #Variaveis
       v, #Variaveis
       z  #Variaveis
       )-> float:
    
    [Explicacao] - Esta eh a derivada pacial de F em relacao a variavel z

    #Definindo o primeiro termo:
    firstT  = fun.fwc(u, v, z) * (v - v0)
    #Definindo o segundo termo:
    secondT = fun.foc(u, v, z) * (u - u0)

    return firstT - secondT

def gradF(u0, v0, z0, u, v, z) -> list:
    return [Fu(u0, v0, z0, u, v, z), Fv(u0, v0, z0, u, v, z), Fz(u0, v0, u, v, z)]"""

##############################
def Gu(
       alpha,   #Variavel de controle
       u0,      #u inicial
       v0,      #v inicial
       z0,      #z inicial
       u,       #Variaveis
       v,       #Variaveis
       z        #Variaveis
       )-> float:
    
    """
    [Explicacao] - Esta eh a derivada pacial de G em relacao a variavel u
    """

    #Definindo o primeiro termo:
    firstT  = fun.fwu(u, v, z) * (u0 * (z - z0) + alpha * (af.a(z) - af.a(z0)))
    
    #Definindo o segundo termo:
    secondT = fun.fw(u0, v0, z0) * (z - z0)

    return firstT - secondT

def Gv(
       alpha,   #Variavel de controle
       u0,      #u inicial
       z0,      #z inicial
       u,       #Variaveis
       v,       #Variaveis
       z        #Variaveis
       )-> float:
    

    """
    [Explicacao] - Esta eh a derivada pacial de G em relacao a variavel v
    """

    #Definindo o primeiro termo:
    firstT  = fun.fwv(u, v, z) * (u0*(z - z0) + alpha * (af.a(z) - af.a(z0)))

    return firstT 

"""def Gz(
       alpha,   #Variavel de controle
       u0,      #u inicial
       v0,      #v inicial
       z0,      #z inicial
       u,       #Variaveis
       v,       #Variaveis
       z        #Variaveis
       )-> float:
    
    
    [Explicacao] - Esta eh a derivada pacial de G em relacao a variavel z
    z

    #Definindo o primeiro termo:
    firstT  = fun.fwc(u, v, z) * (u0*(z - z0) + alpha * (af.a(z) - af.a(z0)))

    #Definindo o segundo termo:
    secondT = (fun.fw(u, v, z) - fun.fw(u0, v0, z0)) * (u0 + alpha * af.az(z))

    #Defininfo o terceiro termo
    thirdT = fun.fw(u0, v0, z0) * (u - u0)

    return firstT + secondT - thirdT

def gradG(alpha, u0, v0, z0, u, v, z):
    return [Gu(alpha, u0, v0, z0, u, v, z), Gv(alpha, u0, z0, u, v, z), Gz(alpha, u0, v0, z0, u, v, z)]"""

"""def Det_u(alpha, u0, v0, z0, u, v, z):
    #definindo a diagonal principal
    diag_p = Fv(u0, v0, z0, u, v, z) * Gz(alpha, u0, v0, z0, u, v, z)

    #Definindo a diagonal secundaria
    diag_s = Fz(u0, v0, u, v, z) * Gv(alpha, u0, z0, u, v, z)

    return diag_p - diag_s

def Det_v(alpha, u0, v0, z0, u, v, z):
    #definindo a diagonal principal
    diag_p = Fu(u0, v0, z0, u, v, z) * Gz(alpha, u0, v0, z0, u, v, z)

    #Definindo a diagonal secundaria
    diag_s = Fz(u0, v0, u, v, z) * Gu(alpha, u0, v0, z0, u, v, z)

    return diag_p - diag_s

def Det_z(alpha, u0, v0, z0, u, v, z):
    #Definindo a diagonal principal
    diag_p = Fu(u0, v0, z0, u, v, z) * Gv(alpha, u0, z0, u, v, z)
    
    #Definindo a diagonal secundaria
    diag_s = Fv(u0, v0, z0, u, v, z) * Gu(alpha, u0, v0, z0, u, v, z)

    return diag_p - diag_s

#####################
def Norma(alpha, u0, v0, z0, u, v, z):
    return sqrt((Det_u(alpha, u0, v0, z0, u, v, z)**2) + 
                (Det_v(alpha, u0, v0, z0, u, v, z)**2) + 
                (Det_z(alpha, u0, v0, z0, u, v, z)**2))

def Hug1(alpha, u0, v0, z0, u, v, z):
    return Det_u(alpha, u0, v0, z0, u, v, z)/Norma(alpha, u0, v0, z0, u, v, z)

def Hug2(alpha, u0, v0, z0, u, v, z):
    return -Det_v(alpha, u0, v0, z0, u, v, z)/Norma(alpha, u0, v0, z0, u, v, z)

def Hug3(alpha, u0, v0, z0, u, v, z):
    return Det_z(alpha, u0, v0, z0, u, v, z)/Norma(alpha, u0, v0, z0, u, v, z)


def CampoHug(alpha, u0, v0, z0, u, v, z):
    return [Hug1(alpha, u0, v0, z0, u, v, z),
            Hug2(alpha, u0, v0, z0, u, v, z),
            Hug3(alpha, u0, v0, z0, u, v, z)]"""

#=============================================
# Nova interpretacao do Hugoniot (4 dimensoes)
#=============================================

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
    secondT1 = fun.foc(u, v, z) * (alpha * (af.a(z) - af.a(z0)) + (z - z0) * u0)

    firstT2 = (u - u0) * (fun.fw(u, v, z) - sigma * (u + alpha * af.az(z)))
    secondT2 = fun.fwc(u, v, z) * (alpha * (af.a(z) - af.a(z0)) + (z - z0) * u0)

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
    secondT1 = fun.foc(u, v, z) * (alpha * (af.a(z) - af.a(z0)) + (z - z0) * u0)

    firstT2 = (u - u0) * (fun.fw(u, v, z) - sigma * (u + alpha * af.az(z)))
    secondT2 = fun.fwc(u, v, z) * (alpha * (af.a(z) - af.a(z0)) + (z - z0) * u0)

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

    return (alpha * (af.a(z) - af.a(z0)) + (z - z0) * u0) * (firstT - secondT)

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

"""
def campo_orientado(alpha, uk, vk, zk, sigk, u0, v0, z0, h_sign):
    P = P_sig(alpha, uk, vk, zk, sigk, u0, v0, z0)
    Q = Q_sig(alpha, uk, vk, zk, sigk, u0, v0, z0)
    R = R_sig(alpha, uk, vk, zk, sigk, u0, v0, z0)
    S_val = S(alpha, uk, vk, zk, sigk, u0, v0, z0)

    # Forca a orientacao: se h > 0, queremos R (dz/ds) > 0;
    # se h < 0, queremos R < 0. Se nao bater, inverte o campo INTEIRO.
    if (h_sign > 0 and R < 0) or (h_sign < 0 and R > 0):
        P, Q, R, S_val = -P, -Q, -R, -S_val

    return P, Q, R, S_val

#############################################################################################################
def HugonioutEuler_method(alpha, fixed_point: list, integ_config: list, search_radius=0.01, TOL=1e-3):

    from numpy import array, dot, linalg, clip

    #Extraindo os valores das listas
    h, N = integ_config
    u0, v0, z0 = fixed_point
    sig0 = af.lambdz(u0, v0, z0, alpha)

    #iniciando a lista de pontos
    Points = [[u0, v0, z0, sig0]]

    # === Bootstrap: ponto inicial via solver implicito ===
    zh = z0 + (h/abs(h)) * search_radius
    _, _, _, _, candidatas = af.hugoniotSystemSolver(0.01, 1, 2000, True, u0, v0, z0, zh, alpha)

    if len(candidatas) > 0:
        u1, v1, z1 = min(
            candidatas,
            key=lambda c: linalg.norm([c[0]-u0, c[1]-v0, c[2]-z0])
        )
        sig1 = af.sigma(u1, v1, z0, u0, v0, alpha)
    else:
        print("Cheguei no smoller")
        u1, v1, z1, sig1, _ = smoller_Point(alpha, u0, v0, z0, h, TOL)

    uk, vk, zk, sigk = u1, v1, z1, sig1
    Points.append([uk, vk, zk, sigk])

    # === Inicializa a referencia de orientacao com o PRIMEIRO campo calculado ===
    P0 = P_sig(alpha, uk, vk, zk, sigk, u0, v0, z0)
    Q0 = Q_sig(alpha, uk, vk, zk, sigk, u0, v0, z0)
    R0 = R_sig(alpha, uk, vk, zk, sigk, u0, v0, z0)
    S0 = S(alpha, uk, vk, zk, sigk, u0, v0, z0)

    campo_anterior = array([P0, Q0, R0, S0], dtype=float)

    # Direcao real observada: do fixed_point ate o ponto do bootstrap
    disp = array([uk - u0, vk - v0, zk - z0], dtype=float)

    # Forca a orientacao do campo para concordar com o deslocamento JA CONHECIDO
    if dot(campo_anterior[:3], disp) < 0:
        campo_anterior = -campo_anterior

    #Flag para verificar 
    saiu_da_vizinhanca_z0 = False

    # === Parametros do passo adaptativo (regiao de coincidencia caracteristica) ===
    ESPACIAL_MIN = 1e-3   # limiar abaixo do qual consideramos "quase parado" em (u,v,z)
    H_MAX_FACTOR = 20      # limite de amplificacao do passo, em multiplos de |h|

    # === Loop de Euler com orientacao por continuidade + passo adaptativo ===
    for _ in range(N):
        P = P_sig(alpha, uk, vk, zk, sigk, u0, v0, z0)
        Q = Q_sig(alpha, uk, vk, zk, sigk, u0, v0, z0)
        R = R_sig(alpha, uk, vk, zk, sigk, u0, v0, z0)
        Sv = S(alpha, uk, vk, zk, sigk, u0, v0, z0)

        campo_atual = array([P, Q, R, Sv], dtype=float)

        # Continuidade: compara com o passo anterior, nao com um sinal fixo
        if dot(campo_atual, campo_anterior) < 0:
            campo_atual = -campo_atual

        campo_anterior = campo_atual
        P, Q, R, Sv = campo_atual

        # --- Passo adaptativo: amplifica h quando o avanco espacial (u,v,z) e pequeno ---
        espacial = (P**2 + Q**2 + R**2) ** 0.5

        if espacial < ESPACIAL_MIN:
            fator = ESPACIAL_MIN / max(espacial, 1e-8)
            fator = min(fator, H_MAX_FACTOR)
            h_efetivo = h * fator
        else:
            h_efetivo = h

        ukp1 = uk + h_efetivo * P
        vkp1 = vk + h_efetivo * Q
        zkp1 = zk + h_efetivo * R
        sigkp1 = sigk + h_efetivo * Sv

        Points.append([ukp1, vkp1, zkp1, sigkp1])

        if not saiu_da_vizinhanca_z0 and abs(zkp1 - z0) > TOL:
            saiu_da_vizinhanca_z0 = True

        if saiu_da_vizinhanca_z0 and abs(zkp1 - z0) < TOL:
            zh2 = 0
            if(h > 0):
                zh2 = z0 - (h/abs(h)) * search_radius
            else:
                zh2 = z0 + (h/abs(h)) * search_radius
            _, _, _, _, candidatas2 = af.hugoniotSystemSolver(0.01, 1, 5000, True, u0, v0, z0, zh2, alpha)

            if len(candidatas2) > 0:
                u_star, v_star, _ = min(
                    candidatas2,
                    key=lambda c: abs(af.sigma(c[0], c[1], z0, u0, v0, alpha) - sigkp1)
                )
                Points[-1] = [u_star, v_star, z0, sigkp1]

            continue

        uk, vk, zk, sigk = ukp1, vkp1, zkp1, sigkp1

    Points = array(Points, float)
    return [[[p[0], p[1], p[2]], p[3]] for p in Points]

"""

def smoller_Point(
                alpha : float, #Variavel de controle
                u0    : float, #Componente u do ponto fixo
                v0    : float, #Componente v do ponto fixo
                z0    : float, #Componente z do ponto fixo
                sig0  : float, #Componente sig do ponto fixo
                h     : float, #Passo
                TOL   : float  #Tolerancia para a norma do Maximo
                ) -> list:

    #Cria uma copia do passo
    h0 = h

    #Variavel para o criterio de parada 
    MAX_ITERATOR = 1000
    it = 0 

    while(True):
        uk = u0 + h0 * ce.P(u0, v0, z0, alpha)
        vk = v0 + h0 * ce.Q(u0, v0, z0, alpha)
        zk = z0 + h0 * ce.R(u0, v0, z0, alpha)
        sigk = af.lambdz(u0, v0, z0, alpha) + h0 * 0.5 #o 0.5 veio do Smoller

        #Norma do maximo
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

def HugoniotPlanaImplicita(
        alpha : float, #Variavel de controle
        u0    : float, #Componente u do ponto fixo
        v0    : float, #Componente v do ponto fixo
        z0    : float, #Componente z do ponto fixo
        TOL   : float  #Tolerancia para a norma do Maximo
):  
    def H(u,v):
        return (fun.fw(u,v,z0) - fun.fw(u0,v0,z0))*(v-v0) - (fun.fo(u,v,z0) - fun.fo(u0,v0,z0))*(u-u0)

    def sig(u, v):
        FT = ((fun.fw(u,v,z0) * (v - v0))/(u0 + alpha * af.az(z0)) - (fun.fo(u,v,z0) - fun.fo(u0,v0,z0)))*(u - u0)
        ST = H(u, v)
        return FT - ST

    #Define um grid no quadrado [0, 1] X [0, 1]
    x = linspace(0.0, 1.0, 750)
    y = linspace(0.0, 1.0, 750)
    X, Y = meshgrid(x, y)

    Hl = H(X, Y)                    #Hl armazena os pontos em todo o grid [0, 1] X [0, 1]
    Hl = where(X + Y < 1, Hl, nan)  #Agora Hl armazena os pontos no triangulo X + Y < 1

    #Guarde a referencia da curva na variavel 'cs'
    cs = contour(X, Y, Hl, [0], linestyles="dashed", colors='k')

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
        mask = array([bordaTriangulo(p) for p in poly]) #o varremos para encontrar os pontos-limites

        #Caso incomum (Professor, isso eh possivel?)
        #Se a polilinha nunca toca a borda, ela e um ramo fechado (loop) inteiro DENTRO do triangulo
        if(not mask.any()):
            branches.append(poly)
            continue

        #Caso normal: a polilinha vai de uma borda a outra sem tocar
        #a borda no meio -> e um unico ramo (comportamento tipico do contour)
        indices_borda = where(mask)[0] #Extrai os indices dos pontos na polilinha onde a condicao eh verdadeira

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
        return sig(u, v)

    resultados = []  #lista de (idx_do_ramo, [pontos onde sig(u,v) = 0])
    for idx_ramo, ramo in enumerate(branches):
        print(f"ramo {idx_ramo}: {len(ramo)} pontos")
        sig_vals = array([sig(p[0], p[1]) for p in ramo]) #Armazena o valor dos sigmas de CADA PONTO do ramo 

        raizes_do_ramo = []                               #Armazena os pontos onde sig(u, v) = 0 na forma [idx, [pontos]]
        for i in range(len(ramo) - 1):
            s0, s1 = sig_vals[i], sig_vals[i+1]           #Extrai o valor do sigma do i-esimo ponto e seu sucessor

            if isnan(s0) or isnan(s1):                    #Se ui ou vi sao nan
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

def HugoniotNewtonRaphson_method(
        alpha           : float,
        current_point   : list,
        fixed_point     : list,
        tol      = 1e-8,
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
        F = F(u0, v0, z0, u, v, z1)
        G = G(alpha, u0, v0, z0, u, v, z1)

        #Avaliacao das derivadas parciais(Jacobiano 2x2)
        Fu = Fu(u0, v0, z0, u, v, z1)
        Fv = Fv(u0, v0, z0, u, v, z1)
        Gu = Gu(alpha, u0, v0, z0, u, v, z1)
        Gv = Gv(alpha, u0, z0, u, v, z1)

        #Determinante do Jacobiano
        dJ = Fu * Gv - Fv * Gu
        # print(f"Valor do det|J(u,v)| = {dJ}")

        if(abs(dJ) < 1e-12):
            raise ValueError("Jacobiano singular ou muito próximo de zero.")

        #Solucao exata da inversao da matriz 2x2: J * delta = -F
        du = gama * ((Fv * G - Gv * F) / dJ)
        dv = gama * ((Gu * F - Fu * G) / dJ)

        return du, dv, F, G

    #Definindo o ponto de partida no plano z = z1
    uk, vk = u1, v1

    for _ in range(max_iter):
        du, dv, F, G = deltaUV(uk, vk)

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
            tol = 1e-5,             #Tolerancia default
            flag = True,            #Flag para ativar ou nao o smoller
):
    #Extraindo ponto fixo:
    u0, v0, z0 = fixed_point
    sig0 = af.lambdz(u0, v0, z0, alpha) #sigma inicial

    #Extraindo configuracao da integracao
    h, N = integ_config

    #iniciando a lista de pontos
    Points = [[u0, v0, z0, sig0]] if(flag) else []

    #Iniciando os pontos iniciais para integracao
    uSmoller, vSmoller, zSmoller, sigSmoller = smoller_Point(alpha, u0, v0, z0, sig0, h, tol)

    uk, vk, zk, sigk = (uSmoller, vSmoller, zSmoller, sigSmoller) if(flag) else (u0, v0, z0, sig0 )
    Points.append([uk, vk, zk, sigk])

    #nova tolerancia:
    zTol = abs(zSmoller - z0)
    print(f"Tolerância em relação a zTol: {zTol}")

    if(R_sig(alpha, uk, vk, zk, sigk, u0, v0, z0) > 0):
        if(((zk > z0) and (h < 0)) or ((zk < z0) and (h > 0))):
            h = -h
    else:
        if(((zk > z0) and (h > 0)) or ((zk < z0) and (h < 0))):
            h = -h

    for _ in range(N):
        ukp1    = uk   +  h * P_sig(alpha, uk, vk, zk, sigk, u0, v0, z0)
        vkp1    = vk   +  h * Q_sig(alpha, uk, vk, zk, sigk, u0, v0, z0)
        zkp1    = zk   +  h * R_sig(alpha, uk, vk, zk, sigk, u0, v0, z0)
        sigkp1  = sigk +  h * S(alpha, uk, vk, zk, sigk, u0, v0, z0)
        
        #Appendo o ponto
        Points.append([ukp1, vkp1, zkp1, sigkp1])

        #Se o ponto atual estiver proximo do plano z = z0 ou z = 0
        if((abs(zkp1 - z0) < zTol) or (abs(zkp1) < tol) or (abs(zkp1 - 1.0) < tol)):
            break;

        #Atualizo os componentes usando Euler
        uk, vk, zk, sigk = ukp1, vkp1, zkp1, sigkp1

    for i in reversed(range(len(Points))):
        u, v, z, sig = Points[i]

        if(abs(F(u0, v0, z0, u, v, z)) < tol and abs(G(alpha, u0, v0, z0, u, v, z)) < tol):
            break

        try:
            u, v = HugoniotNewtonRaphson_method(alpha, [u, v, z], fixed_point)
            sig  = S(alpha, u, v, z, sig, u0, v0, z0)
            Points[i] = [u, v, z, sig]
        except Exception as e:
            print(f"[LOG] - {e}")
            

    Points = array(Points, float)                       #Transformo em um array numpy
    return [[[p[0], p[1], p[2]], p[3]] for p in Points] #Retorno da forma [[u, v, z], sig]
