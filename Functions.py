# -*- coding: utf-8 -*-
"""
VsCode Editor

Definicao das principais funcoes dependentes do modelo

fronteira(u,v): testa se o ponto (u, v) estah fora da regiao de interesse

map(x, y, mp): mapeamento para o  triangulo equilatero
umap(x, y, mp): des-mapeamento do triangulo equilatero

muw(c): Viscosidade da agua

D(u,v,c): Denominador
Du(u,v,c): dD/du
Dv(u,v,c): dD/dv
Dc(u,v,c): dD/dc

fw(u,v,c): funcao fluxo agua
fwu(u,v,c): dfw/du
fwv(u,v,c): dfw/dv
fwc(u,v,c): dfw/dc

fo(u,v,c): funcao fluxo oleo
fou(u,v,c): dfo/du
fov(u,v,c): dfo/dv
foc(u,v,c): dfo/dc

fw(u,v,c): funcao fluxo gas
fou(u,v,c): dfg/du
fov(u,v,c): dfg/dv
foc(u,v,c): dfg/dc

lbdc(u,v,c): lambdac (autovalor contato)
lbdas(u,v,c): lambdas (autovalow slow)
lbdaf(u,v,c): lambdaf (autovalow fast)
eigvls(u,v,c): ambos autovalores (slow e fast)
"""

import numpy as np  #funcoes matemáticas
import Inicia as ini

# Concentracoes iniciais
cmin, cmax = ini.concentrations()

muw0, muo, mug = ini.viscosidades()
 
#Mapeamento
def map(x, y, mp): #Funcao que faz o mapeamento do triangulo equilatero
    if mp == 1:
        x = x + 0.5*y
        y = np.sqrt(3)/2*y
    return x,y

#Desmapeamento
def umap(x, y, mp): #Funcao que faz o des-mapeamento do triangulo equilatero
    if mp == 1:
        sqrt3 = np.sqrt(3)
        x = x - 1/sqrt3 * y
        y = 2/sqrt3*y
    return x,y

#Viscosidade da agua em funcao da concetracao do polimero   
def muw(c): #Viscosidade da agua
    #return muw0 + c
    return muw0*2**c

def muwc(c): #dmuw/dc
    #return 1.0
    ln2 = np.log(2)
    return ln2*muw0*2**c

#Denominador (mobilidade total)
def D(u,v,c): #Denominador
    w = 1 - u - v
    return u**2/muw(c) + v**2/muo + w**2/mug

def Du(u,v,c): #dD/du
    w = 1 - u - v
    return 2*u/muw(c) - 2*w/mug

def Dv(u,v,c): #dD/dv
    w = 1 - u - v
    return 2*v/muo - 2*w/mug

def Dc(u,v,c): ##dD/dc
    return -muwc(c)*u**2/muw(c)**2

#Funcoes de Fluxo

def fw(u,v,c): #fw
    return (u**2/muw(c))/D(u,v,c)

def fwu(u,v,c): #dfw/du
    denom = D(u, v, c)
    viscw = muw(c)
    return (2*u/viscw*denom - u**2/viscw*Du(u,v,c)) / (denom**2)

def fwv(u,v,c): #dfw/dv
    return (- u**2/muw(c)*Dv(u,v,c)) / (D(u,v,c)**2)

def fwc(u,v,c): #dfw/dc
    denom = D(u, v, c)
    viscw = muw(c)
    return u**2*(- (muwc(c)/viscw**2)*denom - Dc(u,v,c)/viscw) / (denom**2)
                        
def fo(u,v,c): #fo
    return (v**2/muo)/D(u,v,c)

def fou(u,v,c): #dfo/du
    return (-v**2/muo)*Du(u,v,c) / (D(u,v,c)**2)

def fov(u,v,c): #dfo/dv
    denom = D(u, v, c)
    return (2*v/muo*denom - v**2/muo*Dv(u,v,c)) / (denom**2)

def foc(u,v,c): #dfo/dc
    return -v**2/muo*Dc(u,v,c) / (D(u,v,c)**2)

def fg(u,v,c): #fg
    w = 1 - u - v
    return (w**2/mug)/D(u,v,c)

def fgu(u,v,c): #dfg/du
    denom = D(u, v, c)
    w = 1 - u - v
    return (-2*w/mug*denom-w**2/mug*Du(u,v,c))/(denom**2)

def fgv(u,v,c): #dfg/dv
    denom = D(u, v, c)
    w = 1 - u - v
    return (-2*w/mug*denom-w**2/mug*Dv(u,v,c))/(denom**2)

def fgc(u,v,c): #dfg/dv
    w = 1 - u - v
    return (-w**2/mug*Dc(u,v,c))/(D(u,v,c)**2)

################################################
# Autovalor contato lambdac 
#(Alpha == 0)
def lbdc(u,v,c): #lambdac = fw/sw
    return (u/muw(c))/D(u,v,c) # fw(u,v,c)/u

################################################
#Autovalores nos planos
    
 
# Os dois autovalores separados
def lbdas(u,v,c):
    #Jacobiana de A
    a11 = fwu(u,v,c)
    a12 = fwv(u,v,c)
    a21 = fou(u,v,c)
    a22 = fov(u,v,c)
    
    lambdas = 0.5*(a11 + a22 - np.sqrt((a22 - a11)**2 + 4*a21*a12))
          
    return lambdas

def lbdaf(u,v,c):
    #Jacobiana de A
    a11 = fwu(u,v,c)
    a12 = fwv(u,v,c)
    a21 = fou(u,v,c)
    a22 = fov(u,v,c)
    
    lambdaf = 0.5*(a11 + a22 + np.sqrt((a11 - a22)**2 + 4*a21*a12))
    
    return lambdaf

#Os dois autovalores juntos
def eigvls(u,v,c):
    #Jacobiana de A
    a11 = fwu(u,v,c)
    a12 = fwv(u,v,c)
    a21 = fou(u,v,c)
    a22 = fov(u,v,c)
    traco = a11 + a22
    discr = np.sqrt((a22 - a11)**2 + 4*a21*a12)
    
    lambdas = 0.5*(traco - discr)
    lambdaf = traco - lambdas
    
    return lambdas, lambdaf

