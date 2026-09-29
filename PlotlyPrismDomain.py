import sys

import plotly.graph_objects as go
import includes.Functions as fun
import includes.Inicia as ini

from PySide6.QtWidgets import QApplication, QMainWindow
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtNetwork import QTcpServer, QHostAddress
from PySide6.QtCore import QUrl

class Prisma:
    """Geometria do prisma (arestas), a partir do mapeamento baricentrico."""

    def __init__(self): #Inicializacao do prisma
        mp = ini.baricentrica()
        u_min, u_max = 0.0, 1.0
        v_min, v_max = 0.0, 1.0
        z_min, z_max = 0.0, 1.0

        self.I = [  #pontos inferiores do prisma
            (*fun.map(u_min, v_min, mp), z_min),
            (*fun.map(u_max, v_min, mp), z_min),
            (*fun.map(u_min, v_max, mp), z_min),
            (*fun.map(u_min, v_min, mp), z_min),
        ]
        self.S = [  #pontos superiores do prisma
            (*fun.map(u_min, v_min, mp), z_max),
            (*fun.map(u_max, v_min, mp), z_max),
            (*fun.map(u_min, v_max, mp), z_max),
            (*fun.map(u_min, v_min, mp), z_max),
        ]

    def build_figure(self, points=None):
        #Combino as coordenaas de cada ponto do respectivo conjunto
        #Cada variavel fica de forma analoga a Iu = (u_min, u_max, u_min, u_min) 
        Iu, Iv, Iz = zip(*self.I)
        Su, Sv, Sz = zip(*self.S)

        fig = go.Figure() #Crio a figura

        fig.add_trace(    #Adiciona novo traco, equivalente ao plt do matlpotlib
            go.Scatter3d( #A nova adicao eh um objeto do tipo scatter (no R3)
                x=Iu,     #Equivalente ao matplotlib
                y=Iv,
                z=Iz,
                mode="lines", #Configuracao do tipo de conexao entre os pontos
                line=dict(color="black", width=3) #Cor do ponto
            )
        )

        fig.add_trace(
            go.Scatter3d(
                x=Su,
                y=Sv,
                z=Sz,
                mode="lines",
                line=dict(color="black", width=3)
            )
        )


        for (iu, iv, iz), (su, sv, sz) in zip(self.I, self.S):
            fig.add_trace(
                go.Scatter3d(
                    x=[iu, su],
                    y=[iv, sv],
                    z=[iz, sz],
                    mode="lines",
                    line=dict(color="black", width=3)
                )
            )

        if points:
            xs, ys, zs = zip(*points)
            fig.add_trace(
                go.Scatter3d(
                    x=xs,
                    y=ys,
                    z=zs,
                    mode="markers",
                    marker=dict(size=3, color="red")
                )
            )

        fig.update_layout(
            title=dict(text="Prism domain", font=dict(size=22, family="Arial, sans-serif", color="#000000"), x=0.5, y=0.95),
            scene=dict(camera=dict(up=dict(x=0, y=0, z=1), eye=dict(x=-1.8, y=-1.5, z=1.2))),
        )
        return fig


class PrismaWindow(QMainWindow):
    """
    Janela do prisma. Escuta uma porta TCP; cada linha recebida eh um ponto
    no formato "u, v, z", ex: "0.3, 0.5, 0.2"
    """

    def __init__(self, host="127.0.0.1", port=5050):
        super().__init__()
        self.setWindowTitle("Prisma 3D")
        self.resize(1000, 800)

        self.prisma = Prisma()
        self.points = []

        self.view = QWebEngineView()
        self.setCentralWidget(self.view)
        self._render()

        self._buffers = {}
        self.server = QTcpServer(self)
        self.server.newConnection.connect(self._on_new_connection)
        self.server.listen(QHostAddress(host), port)
        print(f"[PrismaWindow] Escutando em {host}:{port}")

    def _on_new_connection(self):
        while self.server.hasPendingConnections():
            sock = self.server.nextPendingConnection()
            self._buffers[sock] = b""
            sock.readyRead.connect(lambda s=sock: self._on_ready_read(s))
            sock.disconnected.connect(lambda s=sock: self._buffers.pop(s, None))

    def _on_ready_read(self, sock):
        self._buffers[sock] += bytes(sock.readAll())
        while b"\n" in self._buffers[sock]:
            line, self._buffers[sock] = self._buffers[sock].split(b"\n", 1)
            line = line.strip()
            if line:
                self._add_point(line.decode("utf-8"))

    def _add_point(self, linha):
        try:
            x, y, z = (float(v) for v in linha.split(","))
        except ValueError:
            print(f"[PrismaWindow] linha invalida: {linha!r}")
            return
        self.points.append((x, y, z))
        self._render()

    def _render(self):
        fig = self.prisma.build_figure(points=self.points)
        html = fig.to_html(full_html=True, include_plotlyjs="cdn")
        self.view.setHtml(html, QUrl("https://localhost/"))


def main():
    app = QApplication(sys.argv)
    window = PrismaWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()