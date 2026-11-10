"""GGH como CVP (y SVP) sobre la base pública.

Clave privada: una base buena R. Pública: B = U·R, la misma celosía con una base mala.
Cifrar es esconder un punto: c = m·B + e, con e pequeño. El atacante solo ve c y B.

    celosía:   L(B), la de la base pública
    objetivo:  c
    ataque:    la cadena de cadena.py (LLL → Babai → certificado → enumeración)
    lectura:   el punto es m·B; sus pasos en la base pública son el mensaje m      mensaje

Las tareas SVP piden el vector más corto de L(B): la misma cadena, sin objetivo.
Las tareas con «reducción módulo n» (Tier03 T01, T02 y T05) están en ggh_modulo.py.

Convención: las bases se dan por filas (cada fila es un vector).
"""

from cadena import resuelve_cvp, resuelve_svp
from celosias import coordenadas_en_base


def rompe_ggh_cvp(base, objetivo):
    """Punto de L(B) más cercano al objetivo, con su traza, y el mensaje que esconde."""
    traza = resuelve_cvp(base, objetivo)
    return traza, mensaje(base, traza.respuesta)


def rompe_ggh_svp(base):
    """Vector no nulo más corto de L(B), con su traza."""
    return resuelve_svp(base)


def mensaje(base_publica, punto):
    """El mensaje m tal que m·B = punto: los pasos del punto con la base pública (m = punto·B⁻¹)."""
    return coordenadas_en_base(base_publica, punto)
