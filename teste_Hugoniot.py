import includes.Campo_Hugoniot_Teste as cht
import includes.Inicia as ini
import includes.Functions as fun

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

umin, umax, vmin, vmax, triang = ini.dominio()
mp = 0

alpha = 0.0
u0 = 0.5
v0 = 0.3
z0 = 0.5

branches, resultados = cht.HugoniotPlanaImplicita(alpha, u0, v0, z0, 1e-3)

x = np.linspace(umin, umax, 500)
y = np.linspace(vmin, vmax, 500)
X, Y = np.meshgrid(x, y)

# === Plot dos ramos, cada um com uma cor diferente ===
cores = cm.tab10.colors  # paleta com 10 cores distintas, cicla se houver mais ramos

for idx, ramo in enumerate(branches):
    ramo = np.array(ramo)
    cor = cores[idx % len(cores)]
    plt.plot(ramo[:, 0], ramo[:, 1], '-', color=cor, linewidth=1.5,
              label=f'ramo {idx}')

# === Plot dos pontos onde sig(u,v) = 0, por ramo ===
for idx_ramo, raizes in resultados:
    if len(raizes) == 0:
        continue
    raizes = np.array(raizes)
    plt.plot(raizes[:, 0], raizes[:, 1], 'x', color='red', markersize=10,
              markeredgewidth=2, zorder=5)

# === Ponto fixo ===
plt.plot(u0, v0, 'ko', markersize=8, label='ponto fixo')

plt.xlabel('u')
plt.ylabel('v')
plt.title(f'Hugoniot Plana Implícita — z0={z0}, alpha={alpha}')
plt.legend(fontsize=8, loc='best')
plt.grid(True)
plt.show()