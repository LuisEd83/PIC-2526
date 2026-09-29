"""
Modulo: Curva Integral

Objetivos:
-> Enviar curvas integrais (metodo de Euler) para a janela do prisma via socket
"""

import socket

import includes._Branches as b
import includes.Auxiliar_Functions as af
import includes.Functions as fun
from includes.Inicia import baricentrica


def enviar_pontos(pontos, host="127.0.0.1", port=5050):
    """Envia uma lista de pontos (x, y, z) para a janela do prisma, um por linha."""
    msg = "".join(f"{x},{y},{z}\n" for x, y, z in pontos)
    with socket.create_connection((host, port)) as s:
        s.sendall(msg.encode("utf-8"))


def sendIntegralCurve(Point, branches, host="127.0.0.1", port=5050):
    """
    Mapeia os pontos das branches (e o ponto fixo) para a base baricentrica
    e envia tudo para a janela do prisma em uma unica conexao TCP.
    """
    mp = baricentrica()

    pontos = []

    for branch in branches:
        for p in branch:
            u, v, z = p
            u, v = fun.map(u, v, mp)
            p[0] = u
            p[1] = v
            pontos.append((u, v, z))

    #Ponto fixo
    x_bar, y_bar = fun.map(Point[0], Point[1], mp)
    pontos.append((x_bar, y_bar, Point[2]))

    enviar_pontos(pontos, host=host, port=port)


if __name__ == "__main__":
    Point = [0.3, 0.2, 0.2]
    alpha = 0.5
    h = 0.01
    N = 1500
    integ_config = [h, N]

    _, branches = b.Branches_point(alpha, Point, integ_config)
    sendIntegralCurve(Point, branches)