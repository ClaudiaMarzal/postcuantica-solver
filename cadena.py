"""La cadena: un solo ataque para GGH, LWE y Ring-LWE.

Los tres se reducen a buscar un punto de una celosía (CVP) o el vector más corto (SVP);
solo cambian la celosía y el objetivo (ver ggh.py, lwe.py y ring_lwe.py). El ataque es
escalonado: primero lo polinómico y, solo si no basta, lo exponencial.

    CVP (objetivo y)                                      SVP (sin objetivo)
    1. LLL: base buena de la misma celosía                1. LLL
    2. Babai: candidato p, a distancia d de y             2. candidato: la flecha más corta
    3. Certificado: si 4·d² < min ‖b*_i‖², p es el        3. Certificado: si ‖candidato‖² ≤
       más cercano (demostrado, sin buscar)                  min ‖b*_i‖², es el más corto
    4. Si no: enumerar la bola de radio d alrededor       4. Si no: enumerar la bola de radio
       de y y quedarse con el más cercano                    ‖candidato‖ alrededor del origen

Por qué la respuesta es siempre la óptima: LLL no cambia la celosía; o el certificado
demuestra el candidato, o la enumeración mira todos los puntos de la bola (el candidato
está dentro, así que ninguno de fuera puede ser mejor).

Convención: las bases se dan por filas (cada fila es un vector).
"""

from dataclasses import dataclass, field
from time import perf_counter

from celosias import babai, cota_inferior_lambda1, norma2, puntos_en_bola, reduce_lll


@dataclass
class Traza:
    """Resultados intermedios del ataque (la traza que pide el enunciado)."""
    problema: str               # "CVP" o "SVP"
    candidato: list             # punto de Babai (CVP) o flecha más corta de la base LLL (SVP)
    d2_candidato: int           # ‖objetivo − candidato‖² (CVP) o ‖candidato‖² (SVP)
    suelo: float                # min ‖b*_i‖² de la base LLL: cota inferior de λ1²
    certificado: bool           # el suelo demuestra que el candidato es el óptimo
    bola: list | None           # (d², punto) de la enumeración; None si no ha hecho falta
    nodos: int | None           # nodos que ha recorrido la enumeración
    respuesta: list             # el óptimo: el candidato o el mejor de la bola
    tiempos: dict = field(default_factory=dict)   # segundos de cada fase


def cronometra(tiempos, fase, funcion, *argumentos):
    """Ejecuta funcion(*argumentos) y apunta en tiempos[fase] los segundos que ha tardado."""
    inicio = perf_counter()
    resultado = funcion(*argumentos)
    tiempos[fase] = perf_counter() - inicio
    return resultado


def certifica_cvp(d2, suelo):
    """¿Demuestra el suelo que el punto de Babai es el más cercano?

    Otro punto q mejor que p estaría a menos de d del objetivo y, así,
    ‖p − q‖ < 2d. Si 2d < min ‖b*_i‖ ≤ λ1, eso es imposible: no hay puntos tan juntos.
    Se compara al cuadrado (4·d² < suelo) para no sacar raíces.
    """
    return 4 * d2 < suelo


def certifica_svp(d2, suelo):
    """¿Demuestra el suelo que el candidato es el vector más corto? ‖candidato‖² ≤ min ‖b*_i‖² ≤ λ1²."""
    return d2 <= suelo


def resuelve_cvp(base, objetivo):
    """Punto de la celosía más cercano al objetivo: LLL, Babai, certificado y, si hace falta, enumeración."""
    tiempos = {}
    reducida = cronometra(tiempos, "LLL", reduce_lll, base)
    candidato = cronometra(tiempos, "Babai", babai, reducida, objetivo)
    d2 = norma2([y - p for y, p in zip(objetivo, candidato)])
    suelo = cronometra(tiempos, "suelo", cota_inferior_lambda1, reducida)
    if certifica_cvp(d2, suelo):
        return Traza("CVP", candidato, d2, suelo, True, None, None, candidato, tiempos)
    bola, nodos = cronometra(tiempos, "enumeración", puntos_en_bola, reducida, d2, objetivo)
    return Traza("CVP", candidato, d2, suelo, False, bola, nodos, bola[0][1], tiempos)


def resuelve_svp(base):
    """Vector no nulo más corto de la celosía: LLL, certificado y, si hace falta, enumeración."""
    tiempos = {}
    reducida = cronometra(tiempos, "LLL", reduce_lll, base)
    candidato = min((fila for fila in reducida if any(fila)), key=norma2)
    d2 = norma2(candidato)
    suelo = cronometra(tiempos, "suelo", cota_inferior_lambda1, reducida)
    if certifica_svp(d2, suelo):
        return Traza("SVP", candidato, d2, suelo, True, None, None, candidato, tiempos)
    bola, nodos = cronometra(tiempos, "enumeración", puntos_en_bola, reducida, d2)
    return Traza("SVP", candidato, d2, suelo, False, bola, nodos, bola[0][1], tiempos)
