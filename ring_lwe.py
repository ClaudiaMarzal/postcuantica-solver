"""Ring-LWE como CVP: b = a·s + e en Z_q[x]/(x^n + 1), con s y e pequeños.

Los polinomios se guardan como su lista de coeficientes, empezando por x^0. Se trabaja
módulo x^n + 1 (x^n se cambia por −1) y con coeficientes módulo q.

    1. multiplicar por a es multiplicar por una matriz:  a·s = rot(a)·s       matriz_negaciclica
       así que b = rot(a)·s + e es un LWE con m = k = n.
    2. Pero con m = k la celosía q-aria es todo Z^n (no sirve). Falta usar que s
       también es pequeño: s ≡ I·s son n ecuaciones más, de «error» −s.
       A_ext = [I; rot(a)] (2n ecuaciones, n incógnitas), objetivo (0, b)       sistema_ampliado
       (idea de Bai y Galbraith).
    3. la cadena de siempre sobre Λ_q(A_ext), de dimensión 2n, con objetivo (0, b)
                                                       cadena.resuelve_cvp  → punto (s, b − e)
    4. leer el secreto: las n primeras coordenadas del punto son s             rompe_ring_lwe
    5. descifrar la palabra: c2 − s·c1 = w + ⌊q/2⌋·m, con w pequeño            descifra_palabra

Convención: las bases se dan por filas (cada fila es un vector).
"""

from dataclasses import dataclass

from cadena import Traza, resuelve_cvp
from lwe import bit_mas_cercano, celosia_qaria, centrado


def matriz_negaciclica(a):
    """rot(a): la matriz n×n con rot(a)·s = coeficientes de a·s módulo x^n + 1.

    Su columna j son los coeficientes de x^j·a: cada columna es la anterior bajada un
    puesto, y el coeficiente que sale por abajo (x^n) entra arriba con signo menos,
    porque x^n ≡ −1. Es decir, rot(a)_ij = a_(i−j) si i ≥ j y −a_(n+i−j) si i < j.
    """
    n = len(a)
    return [[a[i - j] if i >= j else -a[n + i - j] for j in range(n)] for i in range(n)]


def producto_polinomios(a, s, q):
    """a·s en Z_q[x]/(x^n + 1), como rot(a)·s con los coeficientes módulo q."""
    return [sum(fila_j * sj for fila_j, sj in zip(fila, s)) % q for fila in matriz_negaciclica(a)]


def sistema_ampliado(a, b):
    """El LWE con las ecuaciones de «s pequeño»: A_ext = [I; rot(a)] y objetivo (0, …, 0, b).

    Las n primeras ecuaciones dicen s ≈ 0 (su error es −s); las n últimas, a·s ≈ b
    (su error es e). El punto buscado es (s, b − e), a distancia² ‖s‖² + ‖e‖².
    """
    n = len(a)
    identidad = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    return identidad + matriz_negaciclica(a), [0] * n + list(b)


@dataclass
class RoturaRingLWE:
    s: list                     # secreto, coeficientes centrados (pequeños)
    e: list                     # error centrado, b − a·s
    traza: Traza                # la cadena sobre Λ_q([I; rot(a)]) con objetivo (0, b)
    coherente: bool             # las n últimas coordenadas del punto son a·s (mod q)


def rompe_ring_lwe(a, b, q):
    """La cadena sobre Λ_q([I; rot(a)]), de dimensión 2n, con objetivo (0, b); s son las n primeras coordenadas."""
    n = len(a)
    A_ext, objetivo = sistema_ampliado(a, b)
    traza = resuelve_cvp(celosia_qaria(A_ext, q), objetivo)
    punto = traza.respuesta
    s = [centrado(x, q) for x in punto[:n]]
    a_por_s = producto_polinomios(a, s, q)
    e = [centrado(bi - x, q) for bi, x in zip(b, a_por_s)]
    coherente = all((x - y) % q == 0 for x, y in zip(punto[n:], a_por_s))
    return RoturaRingLWE(s, e, traza, coherente)


def descifra_palabra(s, criptograma, q):
    """Palabra de n bits de (c1, c2): cada coeficiente de c2 − s·c1 está cerca de 0 o de ⌊q/2⌋.

    Con c1 = a·r + e1 y c2 = b·r + e2 + ⌊q/2⌋·m, queda c2 − s·c1 = (e·r + e2 − s·e1) + ⌊q/2⌋·m,
    y el primer sumando es pequeño porque e, r, e2, s y e1 lo son.
    """
    c1, c2 = criptograma
    s_por_c1 = producto_polinomios(c1, s, q)
    return [bit_mas_cercano(x - y, q) for x, y in zip(c2, s_por_c1)]
