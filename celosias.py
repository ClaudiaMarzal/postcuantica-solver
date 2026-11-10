"""Algoritmos de celosías (con fpylll): las piezas con las que se construye la cadena.

    gram_schmidt            la «regla»: los b*_i y sus longitudes
    cociente_hadamard       lo ortogonal que es una base (1 = perpendicular, 0 = muy sesgada)
    reduce_lll              base mala → base buena de la misma celosía (polinómico)
    babai                   objetivo → punto cercano, capa a capa con dos redondeos (polinómico)
    cota_inferior_lambda1   el «suelo» min ‖b*_i‖: λ1 nunca está por debajo (polinómico)
    puntos_en_bola          enumeración: todos los puntos de una bola (exacta, exponencial)
    coordenadas_en_base     los «pasos» de un punto en una base (el mensaje de GGH)

No dependen de GGH ni de LWE: la cadena que las combina está en cadena.py.
Convención: las bases se dan por filas (cada fila es un vector).
"""

from fractions import Fraction
from math import prod, sqrt

from fpylll import GSO, LLL, Enumeration, IntegerMatrix


def norma2(v):
    """‖v‖², la norma al cuadrado: con enteros es exacta y evita raíces."""
    return sum(x * x for x in v)


def gram_schmidt(base):
    """Objeto GSO de fpylll con la ortogonalización de Gram-Schmidt ya calculada.

    M.get_r(i, i) es ‖b*_i‖²: lo que queda de b_i al quitarle su sombra sobre las
    flechas anteriores. Es la medida que usan Hadamard, Babai, el suelo y la enumeración.
    """
    M = GSO.Mat(IntegerMatrix.from_matrix(base))
    M.update_gso()
    return M


def cociente_hadamard(base):
    """Cociente de Hadamard H = (|det B| / ∏‖b_i‖)^(1/n), entre 0 y 1.

    Mide lo ortogonal que es la base: 1 si es ortogonal, cerca de 0 si es muy sesgada.
    Como |det B| = ∏‖b*_i‖ (vectores de Gram-Schmidt), H es la media geométrica de los
    cocientes ‖b*_i‖ / ‖b_i‖ (cuánto sobrevive de cada flecha), y así no hace falta
    calcular el determinante.
    """
    M = gram_schmidt(base)
    n = len(base)
    return prod(sqrt(M.get_r(i, i) / norma2(base[i])) for i in range(n)) ** (1 / n)


def reduce_lll(base):
    """Base LLL-reducida: genera la misma celosía, con vectores cortos y casi ortogonales.

    Solo usa los dos movimientos unimodulares (reducir: restar a una flecha el múltiplo
    entero de otra más cercano a su sombra; intercambiar: si no cumple la condición de
    Lovász), así que no mueve ni un punto. fpylll usa δ = 0.99 y η = 0.51, una
    reducción más fuerte que la δ = 3/4 de las explicaciones.
    """
    B = IntegerMatrix.from_matrix(base)
    LLL.reduction(B, delta=0.99, eta=0.51)
    return [list(fila) for fila in B]


def babai(base, objetivo):
    """Punto de la celosía cercano al objetivo (algoritmo del plano más cercano de Babai).

    Para j = n, …, 1: c_j = redondear(⟨y, b*_j⟩ / ‖b*_j‖²) (en qué capa cae y) y se resta
    c_j·b_j del objetivo. Una pasada, sin volver atrás. Solo da el más cercano de verdad
    si la base es «buena»; por eso se le pasa una base ya reducida con LLL.

    Ojo: CVP.babai de fpylll aplica LLL por dentro, así que aquí se usa MatGSO.babai,
    que trabaja con la base tal como se le da.
    """
    M = gram_schmidt(base)
    coeficientes = M.babai(objetivo)
    return list(M.B.multiply_left(coeficientes))


def cota_inferior_lambda1(base):
    """El «suelo»: cota inferior de λ1² (el cuadrado de la norma del vector no nulo más corto).

    Vale min ‖b*_i‖². Si v = c_1·b_1 + ... + c_j·b_j con c_j ≠ 0 (el último coeficiente no
    nulo), la parte de v en la dirección de b*_j es c_j·b*_j, así que
    ‖v‖ ≥ |c_j|·‖b*_j‖ ≥ ‖b*_j‖. Solo hace falta Gram-Schmidt: tiempo polinómico, sin
    calcular λ1 (que sería resolver el SVP).
    """
    M = gram_schmidt(base)
    return min(M.get_r(i, i) for i in range(len(base)))


def puntos_en_bola(base, radio2, centro=None):
    """Enumeración: todos los puntos de la celosía a distancia² ≤ radio2 del centro.

    Sin centro, la bola está centrada en el origen y se buscan vectores no nulos (SVP);
    de cada pareja ±v sale solo uno. Recorre un árbol con una coordenada por nivel (de la
    última a la primera): en cada nivel solo se miran las capas que cortan la bola, y se
    podan las ramas que ya se salen. Es exacta pero de coste exponencial en la dimensión.
    Va mucho más rápida si la base está reducida con LLL: el árbol depende de las
    longitudes de Gram-Schmidt.

    Devuelve (puntos, nodos): la lista de (d², punto) ordenada por distancia y el número
    de nodos del árbol que se han recorrido.
    """
    M = gram_schmidt(base)
    enumeracion = Enumeration(M, nr_solutions=10_000)
    objetivo = None if centro is None else M.from_canonical(tuple(map(float, centro)))
    # + 0.5: base y centro son enteros, así que las d² también; entran las que valen justo radio2
    soluciones = enumeracion.enumerate(0, len(base), radio2 + 0.5, 0, target=objetivo)
    puntos = []
    for d2, coeficientes in soluciones:
        c = [round(x) for x in coeficientes]
        puntos.append((round(d2), [sum(ci * fila[j] for ci, fila in zip(c, base)) for j in range(len(base[0]))]))
    return sorted(puntos), enumeracion.get_nodes()


def coordenadas_en_base(base, punto):
    """Los «pasos» x con x·B = punto: cuántas veces se usa cada flecha de la base.

    Es resolver el sistema lineal Bᵀ·xᵀ = puntoᵀ (Gauss con fracciones, exacto). Si el
    punto está en la celosía, x es entero. En GGH, con la base pública, x es el mensaje.
    """
    n = len(base)
    ampliada = [[Fraction(base[i][j]) for i in range(n)] + [Fraction(punto[j])] for j in range(n)]
    for col in range(n):
        p = next(i for i in range(col, n) if ampliada[i][col] != 0)
        ampliada[col], ampliada[p] = ampliada[p], ampliada[col]
        pivote = ampliada[col][col]
        ampliada[col] = [x / pivote for x in ampliada[col]]
        for i in range(n):
            if i != col and ampliada[i][col] != 0:
                f = ampliada[i][col]
                ampliada[i] = [x - f * y for x, y in zip(ampliada[i], ampliada[col])]
    x = [fila[n] for fila in ampliada]
    return [int(c) if c.denominator == 1 else c for c in x]
