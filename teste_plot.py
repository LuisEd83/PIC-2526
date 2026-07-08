import numpy as np
import matplotlib.pyplot as plt
from skimage import measure
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

import includes.Campo_Hugoniot as ch

# -------------------------------------------------------
# Flags de plotagem — habilite/desabilite conforme necessario
# -------------------------------------------------------

PLOT_SURFACES = False   # Subplot com F=0, G=0 e intersecao
PLOT_INTERSEC = True   # Subplot apenas com os pontos de intersecao

# -------------------------------------------------------
# Parametros
# -------------------------------------------------------

alpha       = 0.1
u0, v0, z0  = 0.5, 0.3, 0.2
fixed_point = [u0, v0, z0]

# -------------------------------------------------------
# Grade de avaliacao
# -------------------------------------------------------

N        = 575     # Resolucao por eixo
lim      = 0.8    # Raio da janela ao redor do ponto fixo
TOL_surf = 7e-5   # Tolerancia para coleta dos pontos de intersecao

uArr = np.linspace(0.0, 1.1, N)
vArr = np.linspace(0.0, 1.1, N)
zArr = np.linspace(0.0, 1.1, N)

U, V, Z = np.meshgrid(uArr, vArr, zArr, indexing='ij')

F_vals = ch.F(u0, v0, z0, U, V, Z)
G_vals = ch.G(alpha, u0, v0, z0, U, V, Z)

spacing = (uArr[1] - uArr[0],
           vArr[1] - vArr[0],
           zArr[1] - zArr[0])

# -------------------------------------------------------
# Marching Cubes (so se for plotar as superficies)
# -------------------------------------------------------

if PLOT_SURFACES:
    verts_F, faces_F, _, _ = measure.marching_cubes(F_vals, level=0, spacing=spacing)
    verts_G, faces_G, _, _ = measure.marching_cubes(G_vals, level=0, spacing=spacing)

    # Translada para o sistema de coordenadas real
    verts_F += np.array([uArr[0], vArr[0], zArr[0]])
    verts_G += np.array([uArr[0], vArr[0], zArr[0]])

    # ---------------------------------------------------
    # Filtra vertices e faces fora do prisma (u + v <= 1)
    # A face eh descartada se qualquer um dos seus tres
    # vertices estiver fora do prisma
    # ---------------------------------------------------

    mask_vF  = (verts_F[:, 0] + verts_F[:, 1]) <= 1.0
    verts_F  = verts_F[mask_vF]
    faces_F  = np.array([f for f in faces_F if all(mask_vF[f])])

    mask_vG  = (verts_G[:, 0] + verts_G[:, 1]) <= 1.0
    verts_G  = verts_G[mask_vG]
    faces_G  = np.array([f for f in faces_G if all(mask_vG[f])])

    # Reindexacao: apos remover vertices, os indices das
    # faces apontam para posicoes antigas — e necessario
    # remapear para os novos indices contiguos
    old_to_new_F          = np.full(mask_vF.shape, -1)
    old_to_new_F[mask_vF] = np.arange(mask_vF.sum())
    faces_F               = old_to_new_F[faces_F]

    old_to_new_G          = np.full(mask_vG.shape, -1)
    old_to_new_G[mask_vG] = np.arange(mask_vG.sum())
    faces_G               = old_to_new_G[faces_G]

# -------------------------------------------------------
# Coleta dos pontos de intersecao (so se necessario)
# Inclui o filtro do prisma: u + v <= 1
# -------------------------------------------------------

if PLOT_SURFACES or PLOT_INTERSEC:
    mask_intersec   = (
        (np.abs(F_vals) < TOL_surf) &
        (np.abs(G_vals) < TOL_surf) &
        (U + V <= 1.0)               # Limita ao prisma
    )
    intersec_points = np.column_stack([U[mask_intersec],
                                       V[mask_intersec],
                                       Z[mask_intersec]])
    print(f"Pontos na intersecao: {len(intersec_points)}")

# -------------------------------------------------------
# Monta a figura dinamicamente conforme as flags
# -------------------------------------------------------

n_plots = int(PLOT_SURFACES) + int(PLOT_INTERSEC)

if n_plots == 0:
    print("Nenhuma plotagem habilitada.")
    exit()

fig = plt.figure(figsize=(8 * n_plots, 7))
fig.patch.set_facecolor('#0d0d1a')

plot_idx = 1

# -------------------------------------------------------
# Funcao auxiliar: aplica estetica padrao em um eixo 3D
# -------------------------------------------------------

def style_ax(ax):
    ax.set_facecolor('#0d0d1a')
    ax.set_xlim(uArr[0], uArr[-1])
    ax.set_ylim(vArr[0], vArr[-1])
    ax.set_zlim(zArr[0], zArr[-1])
    ax.set_xlabel('u', color='white')
    ax.set_ylabel('v', color='white')
    ax.set_zlabel('z', color='white')
    ax.tick_params(colors='white')
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.pane.set_edgecolor('#333355')
    ax.yaxis.pane.set_edgecolor('#333355')
    ax.zaxis.pane.set_edgecolor('#333355')
    ax.grid(True, color='#333355', linewidth=0.4)

# -------------------------------------------------------
# Subplot 1: superficies F=0 e G=0 + intersecao
# -------------------------------------------------------

if PLOT_SURFACES:
    ax1 = fig.add_subplot(1, n_plots, plot_idx, projection='3d')
    style_ax(ax1)

    if len(faces_F) > 0:
        mesh_F = Poly3DCollection(verts_F[faces_F], alpha=0.20)
        mesh_F.set_facecolor('#4a90d9')
        mesh_F.set_edgecolor('none')
        ax1.add_collection3d(mesh_F)

    if len(faces_G) > 0:
        mesh_G = Poly3DCollection(verts_G[faces_G], alpha=0.15)
        mesh_G.set_facecolor('#e07b54')
        mesh_G.set_edgecolor('none')
        ax1.add_collection3d(mesh_G)

    if len(intersec_points) > 0:
        ax1.scatter(intersec_points[:, 0],
                    intersec_points[:, 1],
                    intersec_points[:, 2],
                    color='lime', s=4, zorder=10)

    ax1.scatter(*fixed_point, color='white', s=50, zorder=5)

    ax1.legend(handles=[
        Patch(facecolor='#4a90d9', alpha=0.5, label='F(u,v,z) = 0'),
        Patch(facecolor='#e07b54', alpha=0.5, label='G(u,v,z) = 0'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='lime',
               markersize=5, label='F=0 ∩ G=0', linestyle='None'),
    ], loc='upper left', facecolor='#0d0d1a', edgecolor='#333355',
       labelcolor='white', fontsize=8)

    ax1.set_title('Superficies F=0 e G=0', color='white', pad=10)
    plot_idx += 1

# -------------------------------------------------------
# Subplot 2: apenas os pontos de intersecao
# -------------------------------------------------------

if PLOT_INTERSEC:
    ax2 = fig.add_subplot(1, n_plots, plot_idx, projection='3d')
    style_ax(ax2)

    if len(intersec_points) > 0:
        ax2.scatter(intersec_points[:, 0],
                    intersec_points[:, 1],
                    intersec_points[:, 2],
                    color='lime', s=6, zorder=10)

    ax2.scatter(*fixed_point, color='white', s=50, zorder=5)

    ax2.legend(handles=[
        Line2D([0], [0], marker='o', color='w', markerfacecolor='lime',
               markersize=5, label='F=0 ∩ G=0', linestyle='None'),
    ], loc='upper left', facecolor='#0d0d1a', edgecolor='#333355',
       labelcolor='white', fontsize=8)

    ax2.set_title('Intersecao: F=0 ∩ G=0', color='white', pad=10)

# -------------------------------------------------------
# Para sobrepor as branches, descomente:
#
# for branch in branches:
#     branch = np.array(branch)
#     if PLOT_SURFACES:
#         ax1.plot(branch[:, 0], branch[:, 1], branch[:, 2],
#                  linewidth=1.5, color='yellow', zorder=15)
#     if PLOT_INTERSEC:
#         ax2.plot(branch[:, 0], branch[:, 1], branch[:, 2],
#                  linewidth=1.5, color='yellow', zorder=15)
# -------------------------------------------------------

plt.tight_layout()
plt.show()