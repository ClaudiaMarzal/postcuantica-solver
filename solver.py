"""Trabajo de criptoanálisis postcuántico (GGH, LWE y Ring-LWE).

Lee los ficheros .txt de las tareas (tareas.py y lector.py) y rompe las tres familias
con un solo ataque, la cadena LLL → Babai → certificado → enumeración (cadena.py), CON certificado.
Solo cambian la celosía y el objetivo:

    GGH        L(B), objetivo c                        ggh.py
    LWE        Λ_q(A), objetivo b                      lwe.py
    Ring-LWE   Λ_q([I; rot(a)]), objetivo (0, b)       ring_lwe.py

Uso (necesita fpylll: dependencias en pyproject.toml):
    uv run solver.py                     # todas las tareas de Csd/Material_asociado
    uv run solver.py fichero.txt ...     # solo esos ficheros (o carpetas)
"""

import sys
from pathlib import Path

from celosias import cociente_hadamard
from ggh import rompe_ggh_cvp, rompe_ggh_svp
from ggh_modulo import resuelve_tarea_modulo
from lector import lee_tarea
from lwe import descifra_bit, rompe_lwe
from ring_lwe import descifra_palabra, rompe_ring_lwe
from tareas import TareaGGH, TareaLWE, TareaRingLWE

CARPETA_TAREAS = Path(__file__).resolve().parents[2] / "Csd" / "Material_asociado"


def rutas_de(argumentos):
    rutas = []
    for a in argumentos or [CARPETA_TAREAS]:
        a = Path(a)
        rutas += sorted(a.rglob("*.txt")) if a.is_dir() else [a]
    return rutas


def lee_tareas(argumentos):
    return [lee_tarea(ruta) for ruta in rutas_de(argumentos)]


# ---------------------------------------------------------------------------
# Cómo se imprime la traza (igual para las tres familias)
# ---------------------------------------------------------------------------

def comprobacion(respuesta, solucion):
    if solucion is None:
        return ""
    return "  ✓" if respuesta == solucion else f"  ✗ esperado {solucion}"


def describe_traza(traza):
    """Una línea con los pasos de la cadena: candidato, suelo, certificado y enumeración."""
    if traza.problema == "SVP":
        texto = f"LLL ‖b1‖² = {traza.d2_candidato}, suelo min ‖b*‖² = {traza.suelo:g}: "
        condicion = f"{traza.d2_candidato} ≤ {traza.suelo:g}"
    else:
        texto = f"Babai d² = {traza.d2_candidato}, suelo min ‖b*‖² = {traza.suelo:g}: "
        condicion = f"4·{traza.d2_candidato} < {traza.suelo:g}"
    if traza.certificado:
        return texto + f"{condicion}, óptimo demostrado sin enumerar"
    distancias = ", ".join(str(d2) for d2, _ in traza.bola)
    mejora = "el candidato ya era el óptimo" if traza.bola[0][0] == traza.d2_candidato else "mejora al candidato"
    return (texto + f"no certifica → bola d² ≤ {traza.d2_candidato}: {len(traza.bola)} punto(s) "
            f"(d² = {distancias}), {traza.nodos} nodos; {mejora}")


def describe_tiempos(tiempos):
    return ", ".join(f"{fase} {segundos * 1000:.2f} ms" for fase, segundos in tiempos.items())


def imprime_traza(traza):
    print(f"    {describe_traza(traza)}")
    print(f"    tiempos: {describe_tiempos(traza.tiempos)}")


# ---------------------------------------------------------------------------
# Una función por familia
# ---------------------------------------------------------------------------

def resuelve_ggh(tarea):
    """Imprime la solución del SVP o del CVP (y el mensaje), la traza y lo ortogonal que era la base pública."""
    if tarea.modulo is not None:
        resuelve_tarea_modulo(tarea, describe_traza, describe_tiempos)   # hipótesis 4, ver ggh_modulo.py
        return
    if tarea.problema == "SVP":
        traza, mensaje = rompe_ggh_svp(tarea.base), None
    else:
        traza, mensaje = rompe_ggh_cvp(tarea.base, tarea.objetivo)
    respuesta = traza.respuesta
    print(f"{tarea.nombre}: {tarea.problema}, n = {tarea.n}, Hadamard = {cociente_hadamard(tarea.base):.3f}, "
          f"respuesta = {respuesta}{comprobacion(respuesta, tarea.solucion)}")
    imprime_traza(traza)
    if mensaje is not None:
        print(f"    mensaje (pasos con la base pública): {mensaje}")


def resuelve_lwe(tarea):
    """Imprime s, el error encontrado, la respuesta (s o el bit) y la traza de la cadena."""
    rotura = rompe_lwe(tarea.A, tarea.b, tarea.q)
    cabecera = f"{tarea.nombre}: LWE, k = {tarea.k}, m = {tarea.m}, q = {tarea.q}"
    if rotura.s is None:
        print(f"{cabecera}: A·s ≡ punto no tiene solución única módulo q")
        imprime_traza(rotura.traza)
        return
    respuesta = descifra_bit(rotura.s, tarea.criptograma, tarea.q) if tarea.criptograma else rotura.s
    print(f"{cabecera}: s = {rotura.s}, max|e| = {max(map(abs, rotura.e))}, respuesta = {respuesta}"
          f"{comprobacion(respuesta, tarea.solucion)}")
    imprime_traza(rotura.traza)
    if rotura.todo_Zm:
        print("    ¡ojo! Λ_q(A) = Z^m (tantas ecuaciones como incógnitas): este s no tiene sentido")


def resuelve_ring_lwe(tarea):
    """Imprime s, el error, la palabra descifrada y la traza de la cadena en dimensión 2n."""
    rotura = rompe_ring_lwe(tarea.a, tarea.b, tarea.q)
    respuesta = descifra_palabra(rotura.s, tarea.criptograma, tarea.q) if tarea.criptograma else rotura.s
    if isinstance(respuesta, list) and tarea.criptograma:
        respuesta_texto = "".join(map(str, respuesta))
    else:
        respuesta_texto = str(respuesta)
    print(f"{tarea.nombre}: Ring-LWE, n = {tarea.n} (celosía de dimensión {2 * tarea.n}), q = {tarea.q}: "
          f"s = {rotura.s}, max|s| = {max(map(abs, rotura.s))}, max|e| = {max(map(abs, rotura.e))}, "
          f"respuesta = {respuesta_texto}{comprobacion(respuesta, tarea.solucion)}")
    imprime_traza(rotura.traza)
    if not rotura.coherente:
        print("    ¡ojo! el punto no es (s, a·s): revisar")


RESOLVEDORES = {TareaGGH: resuelve_ggh, TareaLWE: resuelve_lwe, TareaRingLWE: resuelve_ring_lwe}


if __name__ == "__main__":
    for tarea in lee_tareas(sys.argv[1:]):
        RESOLVEDORES[type(tarea)](tarea)
