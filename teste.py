"""
Seletor interativo de (u, v, z) em Flet.

- Triângulo retângulo isósceles de cateto 1: vértices em (u,v) = (0,0), (1,0), (0,1).
  Passe o mouse sobre o triângulo para pré-visualizar o ponto (marcador claro) e
  clique para confirmar a escolha (marcador sólido).
- Barra vertical ao lado: escolha a altura z em [0, 1] da mesma forma (hover + clique).
- Também é possível digitar a trinca (u, v, z) diretamente e aplicar.

Testado com flet==1.0.0
"""

import flet as ft
import flet.canvas as cv

# ----------------------------------------------------------------------
# Configuração geométrica
# ----------------------------------------------------------------------
SIZE = 320      # lado (em px) do quadrado que contém o triângulo
BAR_W = 60      # largura (em px) da barra de altura
BAR_H = 320     # altura (em px) da barra de altura


def clamp(x: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, x))


def project_to_triangle(u: float, v: float) -> tuple[float, float]:
    """Garante 0<=u, 0<=v, u+v<=1, projetando ortogonalmente na hipotenusa se preciso."""
    u = clamp(u, 0.0, 1.0)
    v = clamp(v, 0.0, 1.0)
    if u + v > 1.0:
        excesso = (u + v - 1.0) / 2.0
        u -= excesso
        v -= excesso
    return u, v


def main(page: ft.Page):
    page.title = "Seleção de (u, v, z)"
    page.window.width = 1000
    page.window.height = 560
    page.horizontal_alignment = ft.CrossAxisAlignment.CENTER
    page.padding = 24

    # -------------------- estado --------------------
    state = {"u": 0.25, "v": 0.25, "z": 0.5}

    # -------------------- marcadores (triângulo) --------------------
    marker = ft.Container(
        width=12, height=12, border_radius=6,
        bgcolor=ft.Colors.RED,
        border=ft.Border.all(1, ft.Colors.WHITE),
        left=0, top=0,
    )
    hover_marker = ft.Container(
        width=10, height=10, border_radius=5,
        bgcolor=ft.Colors.with_opacity(0.45, ft.Colors.RED),
        left=0, top=0, visible=False,
    )

    # -------------------- marcadores (barra de altura) --------------------
    z_marker = ft.Container(width=BAR_W, height=3, bgcolor=ft.Colors.BLUE, left=0, top=0)
    z_hover_marker = ft.Container(
        width=BAR_W, height=2,
        bgcolor=ft.Colors.with_opacity(0.45, ft.Colors.BLUE),
        left=0, top=0, visible=False,
    )

    info_text = ft.Text(size=18, weight=ft.FontWeight.BOLD)
    hover_uv_text = ft.Text(size=12, italic=True, color=ft.Colors.GREY_600, value=" ")
    hover_z_text = ft.Text(size=12, italic=True, color=ft.Colors.GREY_600, value=" ")
    error_text = ft.Text(size=12, color=ft.Colors.RED, value="")

    # -------------------- campos de entrada manual --------------------
    u_field = ft.TextField(label="u", width=90, value="0.250")
    v_field = ft.TextField(label="v", width=90, value="0.250")
    z_field = ft.TextField(label="z", width=90, value="0.500")

    # -------------------- conversões pixel <-> (u,v)/z --------------------
    def uv_to_pixel(u: float, v: float) -> tuple[float, float]:
        # origem (u=0,v=0) no canto inferior esquerdo do quadrado
        x = u * SIZE
        y = SIZE - v * SIZE
        return x, y

    def pixel_to_uv(x: float, y: float) -> tuple[float, float]:
        u = x / SIZE
        v = (SIZE - y) / SIZE
        return project_to_triangle(u, v)

    def z_to_pixel(z: float) -> float:
        return BAR_H - z * BAR_H

    def pixel_to_z(y: float) -> float:
        return clamp((BAR_H - y) / BAR_H, 0.0, 1.0)

    # -------------------- atualização de UI --------------------
    def update_info():
        info_text.value = f"(u, v, z) = ({state['u']:.3f}, {state['v']:.3f}, {state['z']:.3f})"
        u_field.value = f"{state['u']:.3f}"
        v_field.value = f"{state['v']:.3f}"
        z_field.value = f"{state['z']:.3f}"
        page.update()

    def set_uv(u: float, v: float, confirmar: bool):
        x, y = uv_to_pixel(u, v)
        if confirmar:
            marker.left = clamp(x - 6, -6, SIZE - 6)
            marker.top = clamp(y - 6, -6, SIZE - 6)
            state["u"], state["v"] = u, v
        else:
            hover_marker.visible = True
            hover_marker.left = clamp(x - 5, -5, SIZE - 5)
            hover_marker.top = clamp(y - 5, -5, SIZE - 5)

    def set_z(z: float, confirmar: bool):
        y = z_to_pixel(z)
        if confirmar:
            z_marker.top = clamp(y - 1.5, 0, BAR_H - 3)
            state["z"] = z
        else:
            z_hover_marker.visible = True
            z_hover_marker.top = clamp(y - 1, 0, BAR_H - 2)

    # -------------------- handlers: triângulo --------------------
    def tri_hover(e: ft.HoverEvent):
        p = e.local_position
        u, v = pixel_to_uv(p.x, p.y)
        set_uv(u, v, confirmar=False)
        hover_uv_text.value = f"pré-visualização: (u, v) = ({u:.3f}, {v:.3f})"
        page.update()

    def tri_exit(e: ft.HoverEvent):
        hover_marker.visible = False
        hover_uv_text.value = " "
        page.update()

    def tri_tap_down(e: ft.TapEvent):
        p = e.local_position
        u, v = pixel_to_uv(p.x, p.y)
        set_uv(u, v, confirmar=True)
        update_info()

    # -------------------- handlers: barra de altura --------------------
    def z_hover(e: ft.HoverEvent):
        p = e.local_position
        z = pixel_to_z(p.y)
        set_z(z, confirmar=False)
        hover_z_text.value = f"pré-visualização: z = {z:.3f}"
        page.update()

    def z_exit(e: ft.HoverEvent):
        z_hover_marker.visible = False
        hover_z_text.value = " "
        page.update()

    def z_tap_down(e: ft.TapEvent):
        p = e.local_position
        z = pixel_to_z(p.y)
        set_z(z, confirmar=True)
        update_info()

    # -------------------- entrada manual da trinca --------------------
    def apply_manual(e):
        try:
            u = float(u_field.value.replace(",", "."))
            v = float(v_field.value.replace(",", "."))
            z = float(z_field.value.replace(",", "."))
        except (ValueError, AttributeError):
            error_text.value = "Valores inválidos — use números (ex.: 0.25)."
            page.update()
            return
        error_text.value = ""
        u, v = project_to_triangle(u, v)
        z = clamp(z, 0.0, 1.0)
        set_uv(u, v, confirmar=True)
        set_z(z, confirmar=True)
        update_info()

    apply_btn = ft.Button(
        "Definir trinca (u, v, z)",
        icon=ft.Icons.CHECK,
        on_click=apply_manual,
    )

    # -------------------- desenho do triângulo --------------------
    triangle_paint_fill = ft.Paint(
        style=ft.PaintingStyle.FILL,
        color=ft.Colors.with_opacity(0.15, ft.Colors.BLUE),
    )
    triangle_paint_stroke = ft.Paint(
        style=ft.PaintingStyle.STROKE,
        stroke_width=2,
        color=ft.Colors.BLUE,
    )
    triangle_elements = [
        cv.Path.MoveTo(0, SIZE),   # (u=0, v=0)
        cv.Path.LineTo(SIZE, SIZE),  # (u=1, v=0)
        cv.Path.LineTo(0, 0),        # (u=0, v=1)
        cv.Path.Close(),
    ]
    triangle_canvas = cv.Canvas(
        [
            cv.Path(triangle_elements, paint=triangle_paint_fill),
            cv.Path(triangle_elements, paint=triangle_paint_stroke),
        ],
        width=SIZE,
        height=SIZE,
    )

    triangle_detector = ft.GestureDetector(
        content=ft.Stack(
            [triangle_canvas, hover_marker, marker],
            width=SIZE, height=SIZE,
        ),
        mouse_cursor=ft.MouseCursor.CLICK,
        on_hover=tri_hover,
        on_exit=tri_exit,
        on_tap_down=tri_tap_down,
    )

    # -------------------- desenho da barra de altura --------------------
    z_bar = ft.Container(
        width=BAR_W, height=BAR_H,
        gradient=ft.LinearGradient(
            begin=ft.Alignment.TOP_CENTER,
            end=ft.Alignment.BOTTOM_CENTER,
            colors=[
                ft.Colors.with_opacity(0.6, ft.Colors.BLUE),
                ft.Colors.with_opacity(0.05, ft.Colors.BLUE),
            ],
        ),
        border=ft.Border.all(2, ft.Colors.BLUE),
    )

    z_detector = ft.GestureDetector(
        content=ft.Stack(
            [z_bar, z_hover_marker, z_marker],
            width=BAR_W, height=BAR_H,
        ),
        mouse_cursor=ft.MouseCursor.CLICK,
        on_hover=z_hover,
        on_exit=z_exit,
        on_tap_down=z_tap_down,
    )

    # -------------------- layout --------------------
    page.add(
        ft.Row(
            [
                ft.Column(
                    [
                        ft.Text("Triângulo (cateto = 1) — passe o mouse e clique", weight=ft.FontWeight.BOLD),
                        triangle_detector,
                        hover_uv_text,
                    ],
                    horizontal_alignment=ft.CrossAxisAlignment.CENTER,
                ),
                ft.VerticalDivider(),
                ft.Column(
                    [
                        ft.Text("Altura z ∈ [0, 1]", weight=ft.FontWeight.BOLD),
                        z_detector,
                        hover_z_text,
                    ],
                    horizontal_alignment=ft.CrossAxisAlignment.CENTER,
                ),
                ft.VerticalDivider(),
                ft.Column(
                    [
                        ft.Text("Trinca selecionada", weight=ft.FontWeight.BOLD),
                        info_text,
                        ft.Divider(),
                        ft.Text("Entrada manual:"),
                        u_field,
                        v_field,
                        z_field,
                        apply_btn,
                        error_text,
                    ],
                    horizontal_alignment=ft.CrossAxisAlignment.START,
                ),
            ],
            alignment=ft.MainAxisAlignment.CENTER,
            vertical_alignment=ft.CrossAxisAlignment.START,
            spacing=30,
        )
    )

    set_uv(state["u"], state["v"], confirmar=True)
    set_z(state["z"], confirmar=True)
    update_info()


if __name__ == "__main__":
    ft.run(main)