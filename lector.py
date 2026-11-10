"""Lectura de los ficheros .txt de las tareas.

Las líneas que empiezan por «#» son comentarios (y además vienen con la codificación
rota), así que se ignoran. Lo que importa es el orden de las líneas de datos:

    GGH:       problema (SVP/CVP) | objetivo | dimensión n | n filas de la base
               | módulo o N/A | solución o N/D
    LWE:       problema (rotura/criptoanálisis) | criptograma o N/A | k | q
               | m filas de A | b | solución o N/D
    Ring-LWE:  problema | c1 | c2 | n | q | a | b | solución o N/D
"""

import re
import unicodedata
from pathlib import Path

from tareas import TareaGGH, TareaLWE, TareaRingLWE


class ErrorFormato(Exception):
    """El fichero no tiene la estructura esperada."""


# ---------------------------------------------------------------------------
# Utilidades de lectura
# ---------------------------------------------------------------------------

def lineas_de_datos(texto):
    """Líneas que no son comentarios ni están vacías, sin espacios sobrantes."""
    return [l.strip() for l in texto.splitlines() if l.strip() and not l.lstrip().startswith("#")]


def enteros(linea):
    """Todos los enteros de una línea; ignora llaves, comas y espacios (hay filas sin «{»)."""
    return [int(x) for x in re.findall(r"-?\d+", linea)]


def entero(linea):
    numeros = enteros(linea)
    if len(numeros) != 1:
        raise ErrorFormato(f"se esperaba un único número y hay {len(numeros)}: {linea!r}")
    return numeros[0]


def no_disponible(linea):
    return linea.upper() in ("N/A", "N/D")


def sin_tildes(texto):
    return "".join(c for c in unicodedata.normalize("NFKD", texto) if not unicodedata.combining(c))


# ---------------------------------------------------------------------------
# Detección del tipo
# ---------------------------------------------------------------------------

def tipo_por_palabra(datos):
    """Tipo según la primera línea de datos (el «problema a resolver»)."""
    problema = sin_tildes(datos[0]).lower()
    if problema in ("svp", "cvp"):
        return "GGH"
    if "ring" in problema:
        return "Ring-LWE"
    if problema.startswith(("rotura", "cripto")):  # también «criptograma» (Tier0 T04-T05)
        return "LWE"
    return None


def tipo_por_forma(datos):
    """Tipo según la estructura, sin mirar la palabra del problema.

    La dimensión es la primera línea con un solo número. Después:
      - GGH: vienen n filas de n números (base cuadrada).
      - LWE: viene q (otro número suelto), m filas de A, b y la solución.
      - Ring-LWE: viene q, y solo quedan tres líneas: a, b y la solución.
    """
    i = next((j for j in range(1, len(datos)) if len(enteros(datos[j])) == 1), None)
    if i is None or i + 1 >= len(datos):
        return None
    n = entero(datos[i])
    if len(enteros(datos[i + 1])) == n and n > 1:
        return "GGH"
    restantes = len(datos) - (i + 2)
    return "Ring-LWE" if restantes == 3 else "LWE"


def detectar_tipo(datos):
    por_palabra, por_forma = tipo_por_palabra(datos), tipo_por_forma(datos)
    if por_palabra and por_forma and por_palabra != por_forma:
        raise ErrorFormato(f"la palabra dice {por_palabra} pero la estructura parece {por_forma}")
    tipo = por_palabra or por_forma
    if tipo is None:
        raise ErrorFormato(f"no se reconoce el tipo de problema: {datos[0]!r}")
    return tipo


# ---------------------------------------------------------------------------
# Lectores de cada tipo
# ---------------------------------------------------------------------------

def lee_ggh(nombre, datos):
    problema = datos[0].upper()
    n = entero(datos[2])
    base = [enteros(l) for l in datos[3:3 + n]]
    resto = datos[3 + n:]
    if len(base) != n or any(len(fila) != n for fila in base):
        raise ErrorFormato(f"la base debería ser {n}×{n}")
    if len(resto) != 2:
        raise ErrorFormato(f"tras la base se esperaban 2 líneas (módulo y solución) y hay {len(resto)}")

    objetivo = None
    if problema == "CVP":
        objetivo = enteros(datos[1])
        if len(objetivo) != n:
            raise ErrorFormato(f"el objetivo tiene {len(objetivo)} coordenadas y la dimensión es {n}")

    modulo = None if no_disponible(resto[0]) else entero(resto[0])
    solucion = None if no_disponible(resto[1]) else enteros(resto[1])
    return TareaGGH(nombre, problema, objetivo, base, modulo, solucion)


def _cabecera_lwe(datos):
    """Separa las líneas del criptograma (entre el problema y la dimensión) del resto."""
    i = 1
    while i < len(datos) and not no_disponible(datos[i]) and len(enteros(datos[i])) > 1:
        i += 1
    lineas_cripto = [] if no_disponible(datos[i]) else datos[1:i]
    if no_disponible(datos[i]):
        i += 1
    return lineas_cripto, datos[i:]


def lee_lwe(nombre, datos):
    problema = "rotura" if sin_tildes(datos[0]).lower().startswith("rotura") else "criptoanalisis"
    lineas_cripto, resto = _cabecera_lwe(datos)
    k, q = entero(resto[0]), entero(resto[1])
    A = [enteros(l) for l in resto[2:-2]]
    b = enteros(resto[-2])

    if any(len(fila) != k for fila in A):
        raise ErrorFormato(f"todas las filas de A deberían tener k = {k} coeficientes")
    if len(b) != len(A):
        raise ErrorFormato(f"b tiene {len(b)} términos y A tiene {len(A)} filas")

    criptograma = None
    if lineas_cripto:
        numeros = enteros(" ".join(lineas_cripto))
        if len(numeros) != k + 1:
            raise ErrorFormato(f"el criptograma debería ser (u, v) con {k} + 1 números")
        criptograma = (numeros[:-1], numeros[-1])

    sol = resto[-1]
    if no_disponible(sol):
        solucion = None
    elif problema == "rotura":
        solucion = enteros(sol)
    else:
        solucion = entero(sol)
    return TareaLWE(nombre, problema, q, A, b, criptograma, solucion)


def lee_ring_lwe(nombre, datos):
    lineas_cripto, resto = _cabecera_lwe(datos)
    if len(resto) != 5:
        raise ErrorFormato(f"tras el criptograma se esperaban 5 líneas (n, q, a, b, solución) y hay {len(resto)}")
    n, q = entero(resto[0]), entero(resto[1])
    a, b = enteros(resto[2]), enteros(resto[3])
    if len(a) != n or len(b) != n:
        raise ErrorFormato(f"a y b deberían tener n = {n} coeficientes")

    criptograma = None
    if lineas_cripto:
        if len(lineas_cripto) != 2:
            raise ErrorFormato("el criptograma de Ring-LWE debería tener dos líneas (c1 y c2)")
        c1, c2 = enteros(lineas_cripto[0]), enteros(lineas_cripto[1])
        if len(c1) != n or len(c2) != n:
            raise ErrorFormato(f"c1 y c2 deberían tener n = {n} coeficientes")
        criptograma = (c1, c2)

    sol = resto[4]
    if no_disponible(sol):
        solucion = None
    elif set(sol) <= {"0", "1"}:     # palabra de bits escrita seguida: «11010100»
        solucion = [int(c) for c in sol]
    else:
        solucion = enteros(sol)
    return TareaRingLWE(nombre, "criptoanalisis", n, q, a, b, criptograma, solucion)


LECTORES = {"GGH": lee_ggh, "LWE": lee_lwe, "Ring-LWE": lee_ring_lwe}


def lee_tarea(ruta):
    ruta = Path(ruta)
    texto = ruta.read_bytes().decode("latin-1")   # latin-1 nunca falla al decodificar
    datos = lineas_de_datos(texto)
    if not datos:
        raise ErrorFormato("el fichero no tiene datos")
    tipo = detectar_tipo(datos)
    return LECTORES[tipo](ruta.name, datos)
