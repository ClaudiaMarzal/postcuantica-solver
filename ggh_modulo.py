"""GGH con reducción módulo n (Tier03 T01, T02 y T05): hipótesis 4.

La historia que se supone: quien preparó la tarea tenía una base y un criptograma con
números pequeños, positivos y negativos, y antes de publicarlos pasó cada número al
rango de 0 a n − 1 (los redujo módulo n). En T02, con n = 113: −1 → 112, −21 → 92,
−45 → 68.

El ataque deshace ese paso y después hace lo de siempre:

    1. deshacer la reducción: a cada número mayor que n/2 se le resta n     representante_centrado
       (el representante centrado, entre −n/2 y n/2)                         base_centrada, centra
    2. la cadena de siempre con la base y el objetivo centrados            cadena.resuelve_cvp / _svp
    3. (CVP) el mensaje: los pasos del punto en la base centrada            ggh.mensaje

El supuesto: los números originales estaban entre −n/2 y n/2. Si alguno era mayor, el
paso 1 recupera otra base (otra celosía) y la respuesta cambia. Lo único que lo indica
es cuánto quedan los números tras centrar (tamano_tras_centrar): en T02 quedan −1, 16,
15 y −1; en T05, algunos siguen cerca de n/2.

Ojo: rehacer el cifrado con el resultado, (m·B_centrada + e) mod n, siempre reproduce el
objetivo del fichero, porque e se define como lo que falta para llegar a él. No sirve
de comprobación.

Uso:
    uv run ggh_modulo.py                       # las GGH con módulo de Csd/Material_asociado
    uv run ggh_modulo.py fichero.txt ...       # solo esos ficheros
"""

import sys

from cadena import resuelve_cvp, resuelve_svp
from celosias import cociente_hadamard
from ggh import mensaje


def representante_centrado(x, n):
    """El número congruente con x módulo n que está entre −n/2 y n/2: si x mod n > n/2, se le resta n."""
    x %= n
    return x - n if x > n // 2 else x


def centra(v, n):
    return [representante_centrado(x, n) for x in v]


def base_centrada(base, n):
    return [centra(fila, n) for fila in base]


def tamano_tras_centrar(base, n):
    """El mayor |número| de la base centrada, comparado con n/2 (1 = en el límite; cerca de 0 = muy pequeños)."""
    return max(abs(x) for fila in base_centrada(base, n) for x in fila) / (n / 2)


def rompe_ggh_modulo_cvp(base, objetivo, n):
    """Centra la base y el objetivo y aplica la cadena. Devuelve (traza, mensaje, base centrada, objetivo centrado)."""
    Bc, yc = base_centrada(base, n), centra(objetivo, n)
    traza = resuelve_cvp(Bc, yc)
    return traza, mensaje(Bc, traza.respuesta), Bc, yc


def rompe_ggh_modulo_svp(base, n):
    """Centra la base y busca el vector más corto con la cadena. Devuelve (traza, base centrada)."""
    Bc = base_centrada(base, n)
    return resuelve_svp(Bc), Bc


def v(x):
    return "(" + ", ".join(map(str, x)) + ")"


def resuelve_tarea_modulo(tarea, describe_traza, describe_tiempos):
    """Imprime la solución de una GGH con módulo según la hipótesis 4 (lo usa también solver.py)."""
    n = tarea.modulo
    if tarea.problema == "SVP":
        traza, Bc = rompe_ggh_modulo_svp(tarea.base, n)
        extra = ""
    else:
        traza, m, Bc, yc = rompe_ggh_modulo_cvp(tarea.base, tarea.objetivo, n)
        e = [a - b for a, b in zip(yc, traza.respuesta)]
        extra = f", objetivo centrado = {v(yc)}, mensaje = {v(m)}, error = {v(e)}"
    print(f"{tarea.nombre}: {tarea.problema} módulo n = {n} (hipótesis 4), respuesta = {v(traza.respuesta)}, "
          f"reducida módulo n = {v([x % n for x in traza.respuesta])}")
    print(f"    base centrada = {', '.join(map(v, Bc))}{extra}")
    print(f"    tras centrar, el mayor número es el {tamano_tras_centrar(tarea.base, n):.0%} de n/2; "
          f"Hadamard {cociente_hadamard(tarea.base):.3f} → {cociente_hadamard(Bc):.3f}")
    print(f"    {describe_traza(traza)}")
    print(f"    tiempos: {describe_tiempos(traza.tiempos)}")


if __name__ == "__main__":
    from lector import lee_tarea
    from solver import CARPETA_TAREAS, describe_tiempos, describe_traza, rutas_de
    from tareas import TareaGGH
    for ruta in rutas_de(sys.argv[1:] or [CARPETA_TAREAS / "Criptoan_lisis_GGH"]):
        tarea = lee_tarea(ruta)
        if isinstance(tarea, TareaGGH) and tarea.modulo is not None:
            resuelve_tarea_modulo(tarea, describe_traza, describe_tiempos)
