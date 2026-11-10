"""LWE como CVP: del sistema de ecuaciones b = A·s + e (mod q) a una celosía.

Sin el error, los vectores posibles son A·s mod q. Al quitar el «mod» (cada ecuación
gana un entero t_i: b_i = <a_i, s> + e_i + q·t_i) y leer el sistema por columnas,
quedan las combinaciones enteras

    s_1·a_1 + ... + s_k·a_k + t_1·q·e_1 + ... + t_m·q·e_m

de las columnas a_j de A y de los vectores q·e_i. Eso es la celosía q-aria

    Λ_q(A) = A·Z^k + q·Z^m   (dimensión m)

y b es un punto de Λ_q(A) desplazado por el error e: b − e = A·s está en la celosía,
a distancia ‖e‖ de b. Así que encontrar s es un CVP, y el ataque es la misma cadena
que en GGH; solo cambian la celosía, el objetivo y el último paso:

    1. base de Λ_q(A)                         celosia_qaria
    2. LLL → Babai → certificado → (enumeración), con objetivo b
                                              cadena.resuelve_cvp   → v = A·s (mod q)
    3. leer el secreto: resolver A·s ≡ v      resuelve_mod          → s, y el error e = b − v
    4. (criptogramas) descifrar el bit        descifra_bit

Si hay tantas ecuaciones como incógnitas (m = k, A invertible), Λ_q(A) es todo Z^m:
cualquier b ya es un punto, el «error» sale 0 y el s no tiene sentido. Por eso
rompe_lwe avisa (celosia_es_todo_Zm); Ring-LWE lo evita añadiendo ecuaciones (ring_lwe.py).

Convención: las bases se dan por filas (cada fila es un vector), como en GGH.
"""

from dataclasses import dataclass

from cadena import Traza, resuelve_cvp


def gauss_jordan_mod(M, q):
    """Forma escalonada reducida de M módulo q (Gauss-Jordan).

    Devuelve (R, pivotes): R tiene una fila por cada pivote, y en las columnas
    pivote R vale la identidad. Como pivote solo se usan elementos invertibles mod q.
    """
    R = [[x % q for x in fila] for fila in M]
    pivotes = []
    fila = 0
    for col in range(len(R[0])):
        p = next((i for i in range(fila, len(R)) if _invertible(R[i][col], q)), None)
        if p is None:
            continue
        R[fila], R[p] = R[p], R[fila]
        inv = pow(R[fila][col], -1, q)
        R[fila] = [x * inv % q for x in R[fila]]
        for i in range(len(R)):
            if i != fila and R[i][col]:
                f = R[i][col]
                R[i] = [(x - f * y) % q for x, y in zip(R[i], R[fila])]
        pivotes.append(col)
        fila += 1
        if fila == len(R):
            break
    return R[:fila], pivotes


def _invertible(x, q):
    try:
        pow(x, -1, q)
        return True
    except ValueError:
        return False


def celosia_qaria(A, q):
    """Base (por filas) de la celosía q-aria Λ_q(A) = A·Z^k + q·Z^m.

    A es m×k (m ecuaciones, k incógnitas). Los generadores son las k columnas de A y
    los m vectores q·e_i; son k + m vectores en Z^m, así que sobran k. Para quedarse
    con una base:

      1. Se escalona A^T módulo q. Eso equivale a combinar las columnas de A (cambiar
         de incógnitas s' = A_1·s) hasta que en k ecuaciones «pivote» queden la
         identidad. Las k filas resultantes generan lo mismo que las columnas de A
         módulo q.
      2. Los q·e_i de las ecuaciones pivote sobran (se obtienen de esas filas y de los
         demás q·e_j), así que solo se añaden los q·e_j de las ecuaciones no pivote.

    Resultado: m vectores de Z^m. Si se ordenan las coordenadas poniendo primero las
    pivote, la base es

        [ I_k   A'_2^T  ]
        [  0   q·I_m-k  ]

    y su determinante es q^(m-k). Las coordenadas se dejan en el orden original de las
    ecuaciones, para que b se pueda usar tal cual.
    """
    m = len(A)
    At = [list(col) for col in zip(*A)]          # k × m
    R, pivotes = gauss_jordan_mod(At, q)
    base = [list(fila) for fila in R]
    for j in range(m):
        if j not in pivotes:
            base.append([q if i == j else 0 for i in range(m)])
    return base


def resuelve_mod(A, v, q):
    """Solución s de A·s ≡ v (mod q), con Gauss-Jordan sobre la matriz ampliada [A | v].

    Devuelve None si el sistema no tiene solución única módulo q (algún pivote cae en la
    columna de v, o falta un pivote en las de A).
    """
    k = len(A[0])
    ampliada = [fila + [vi] for fila, vi in zip(A, v)]
    R, pivotes = gauss_jordan_mod(ampliada, q)
    if pivotes != list(range(k)):
        return None
    return [R[j][k] for j in range(k)]


def centrado(x, q):
    """Representante de x mod q en (−q/2, q/2]: así un error −1 se ve como −1 y no como q − 1."""
    x %= q
    return x - q if x > q // 2 else x


def celosia_es_todo_Zm(A, q):
    """¿Es Λ_q(A) = Z^m? Pasa si A tiene rango m módulo q (p. ej. m = k con A invertible).

    Entonces todas las ecuaciones son pivote, la base q-aria es la identidad y no hay
    ningún q·e_j: cualquier objetivo es ya un punto de la celosía y el ataque no informa.
    """
    At = [list(col) for col in zip(*A)]
    _, pivotes = gauss_jordan_mod(At, q)
    return len(pivotes) == len(A)


@dataclass
class RoturaLWE:
    s: list | None              # secreto, en [0, q − 1]; None si A·s ≡ v no tiene solución única
    e: list | None              # error centrado, b − A·s
    traza: Traza                # la cadena sobre Λ_q(A) con objetivo b
    todo_Zm: bool               # Λ_q(A) = Z^m: el resultado no tiene sentido


def rompe_lwe(A, b, q):
    """Ataque primal: la cadena sobre Λ_q(A) con objetivo b y, después, Gauss módulo q."""
    traza = resuelve_cvp(celosia_qaria(A, q), b)
    v = traza.respuesta
    s = resuelve_mod(A, v, q)
    e = None if s is None else [centrado(bi - vi, q) for bi, vi in zip(b, v)]
    return RoturaLWE(s, e, traza, celosia_es_todo_Zm(A, q))


def bit_mas_cercano(x, q):
    """0 si x (mod q) está más cerca de 0, 1 si está más cerca de ⌊q/2⌋ (la lectura del reloj)."""
    x %= q
    return 0 if min(x, q - x) < abs(x - q // 2) else 1


def descifra_bit(s, criptograma, q):
    """Bit del criptograma de Regev (u, v): v − <u, s> = rᵀe + bit·⌊q/2⌋ está cerca de 0 o de ⌊q/2⌋."""
    u, v = criptograma
    return bit_mas_cercano(v - sum(ui * si for ui, si in zip(u, s)), q)
