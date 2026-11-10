"""Tipos de tarea: GGH, LWE y Ring-LWE."""

from dataclasses import dataclass


@dataclass
class TareaGGH:
    nombre: str
    problema: str               # "SVP" o "CVP"
    objetivo: list | None       # punto objetivo (None en SVP)
    base: list                  # n vectores de n coordenadas, uno por fila
    modulo: int | None          # «reducción módulo n» (None si N/A)
    solucion: list | None       # None si N/D

    @property
    def n(self):
        return len(self.base)


@dataclass
class TareaLWE:
    nombre: str
    problema: str               # "rotura" o "criptoanalisis"
    q: int
    A: list                     # m filas de k coeficientes
    b: list                     # m términos independientes
    criptograma: tuple | None   # (u, v): u vector de k elementos, v entero
    solucion: list | int | None  # secreto s (rotura) o bit (criptoanálisis)

    @property
    def k(self):
        return len(self.A[0])

    @property
    def m(self):
        return len(self.A)


@dataclass
class TareaRingLWE:
    nombre: str
    problema: str               # "criptoanalisis"
    n: int                      # grado: se trabaja en Z_q[x]/(x^n + 1)
    q: int
    a: list                     # coeficientes de a, empezando por x^0
    b: list
    criptograma: tuple | None   # (c1, c2), cada uno con n coeficientes
    solucion: list | None       # palabra de bits
