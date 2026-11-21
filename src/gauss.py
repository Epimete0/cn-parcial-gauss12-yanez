# src/gauss.py
# Implementacion de eliminacion gaussiana:
# 1) Sin pivoteo (para mostrar el problema del pivote a(1,1) minusculo)
# 2) Con pivoteo parcial escalado (para recuperar mantisa y estabilidad)
# Se incluye sustitucion regresiva, calculo de factor de crecimiento y traza.

import csv
import numpy as np

def leer_matriz_csv(ruta_archivo):
    # Leemos la matriz A desde el archivo csv
    matriz = []
    with open(ruta_archivo, 'r', encoding='utf-8') as archivo:
        lector = csv.reader(archivo)
        for fila in lector:
            if not fila: # por si hay lineas vacias
                continue
            # Convertimos cada elemento a float de 64 bits (doble precision)
            fila_numeros = []
            for elemento in fila:
                fila_numeros.append(float(elemento.strip()))
            matriz.append(fila_numeros)
    return np.array(matriz, dtype=np.float64)


def leer_vector_csv(ruta_archivo):
    # Leemos el vector b desde el archivo csv
    vector = []
    with open(ruta_archivo, 'r', encoding='utf-8') as archivo:
        for linea in archivo:
            linea_limpia = linea.strip()
            if linea_limpia:
                vector.append(float(linea_limpia))
    return np.array(vector, dtype=np.float64)


def formatear_matriz_texto(A, paso):
    # Genera un bloque de texto formateado para la traza de eliminacion
    n = len(A)
    lineas = []
    lineas.append(f"============================================================")
    lineas.append(f"   ESTADO DE LA MATRIZ TRAS EL PASO {paso} DE ELIMINACIÓN")
    lineas.append(f"============================================================")
    for i in range(n):
        fila_str = " ".join(f"{A[i, j]:14.6e}" for j in range(n))
        lineas.append(f"Fila {i+1:2d} | {fila_str}")
    lineas.append("")
    return "\n".join(lineas)


def eliminacion_sin_pivote(A_entrada, b_entrada, contador=None):
    # Eliminacion gaussiana simple sin pivoteo.
    # Copiamos para no alterar los datos originales
    A = A_entrada.copy()
    b = b_entrada.copy()
    n = len(b)
    
    max_multiplicador = 0.0
    max_elemento_global = np.max(np.abs(A))
    max_elemento_inicial = max_elemento_global
    
    trazas = []
    
    # Recorremos las columnas de 0 a n-2 (los 11 pasos)
    for k in range(n - 1):
        # Pivote actual en la diagonal
        pivote = A[k, k]
        
        # Eliminamos hacia abajo desde la fila k+1 hasta la ultima
        for i in range(k + 1, n):
            # OJO: Calculamos el multiplicador una sola vez por fila (exigencia T1)
            # m_ik = a_ik / a_kk
            m = A[i, k] / pivote
            
            if contador:
                contador.registrar_division_multiplicador(k)
            
            # Guardamos el multiplicador mas grande que aparezca
            if abs(m) > max_multiplicador:
                max_multiplicador = abs(m)
                
            # Barrido de columnas desde la columna k en adelante (sin tocar ceros previos)
            for j in range(k, n):
                A[i, j] = A[i, j] - m * A[k, j]
                
            if contador:
                # Se actualizaron (n - k) columnas en la fila i
                contador.registrar_operacion_fila(k, n - k)
            
            # Actualizamos tambien el lado derecho b
            b[i] = b[i] - m * b[k]
            if contador:
                contador.registrar_operacion_b(k)
        
        # Revisamos cuanto crecieron los elementos de la matriz
        max_paso = np.max(np.abs(A))
        if max_paso > max_elemento_global:
            max_elemento_global = max_paso
            
        trazas.append(formatear_matriz_texto(A, k + 1))
        
    factor_crecimiento = max_elemento_global / max_elemento_inicial
    
    return A, b, max_multiplicador, factor_crecimiento, trazas


def eliminacion_pivote_escalado(A_entrada, b_entrada, contador=None):
    # Eliminacion gaussiana con pivoteo parcial escalado
    A = A_entrada.copy()
    b = b_entrada.copy()
    n = len(b)
    
    max_multiplicador = 0.0
    max_elemento_global = np.max(np.abs(A))
    max_elemento_inicial = max_elemento_global
    
    trazas = []
    
    # Vector de escala inicial: el maximo absoluto de cada fila
    # s_i = max_j |a_ij|
    s = np.zeros(n)
    for i in range(n):
        s[i] = np.max(np.abs(A[i, :]))
    
    for k in range(n - 1):
        # Busqueda del mejor pivote escalado entre las filas k a n-1
        # Comparamos |a_ik| / s_i
        fila_mejor = k
        mejor_ratio = abs(A[k, k]) / s[k]
        
        for i in range(k + 1, n):
            if contador:
                contador.registrar_comparacion_pivote(1)
            ratio_actual = abs(A[i, k]) / s[i]
            if ratio_actual > mejor_ratio:
                mejor_ratio = ratio_actual
                fila_mejor = i
                
        # Si la mejor fila no es la fila k, intercambiamos filas completas
        if fila_mejor != k:
            if contador:
                contador.registrar_intercambio_filas()
            # Intercambiamos filas en A, b y en el vector de escalas s
            A[[k, fila_mejor]] = A[[fila_mejor, k]]
            b[[k, fila_mejor]] = b[[fila_mejor, k]]
            s[[k, fila_mejor]] = s[[fila_mejor, k]]
            
        pivote = A[k, k]
        
        # Ahora eliminamos hacia abajo igual que antes
        for i in range(k + 1, n):
            # Calculamos el multiplicador una sola vez por fila
            m = A[i, k] / pivote
            
            if contador:
                contador.registrar_division_multiplicador(k)
                
            if abs(m) > max_multiplicador:
                max_multiplicador = abs(m)
                
            for j in range(k, n):
                A[i, j] = A[i, j] - m * A[k, j]
                
            if contador:
                contador.registrar_operacion_fila(k, n - k)
                
            b[i] = b[i] - m * b[k]
            if contador:
                contador.registrar_operacion_b(k)
                
        max_paso = np.max(np.abs(A))
        if max_paso > max_elemento_global:
            max_elemento_global = max_paso
            
        trazas.append(formatear_matriz_texto(A, k + 1))
        
    factor_crecimiento = max_elemento_global / max_elemento_inicial
    
    return A, b, max_multiplicador, factor_crecimiento, trazas


def sustitucion_regresiva(U, b, contador=None):
    # Resuelve el sistema triangular superior Ux = b hacia atras
    n = len(b)
    x = np.zeros(n, dtype=np.float64)
    
    # Desde la ultima fila (n-1) hacia arriba hasta la primera (0)
    for i in range(n - 1, -1, -1):
        suma_acumulada = 0.0
        operaciones_mult = 0
        operaciones_sum = 0
        
        for j in range(i + 1, n):
            suma_acumulada += U[i, j] * x[j]
            operaciones_mult += 1
            operaciones_sum += 1
            
        # x_i = (b_i - suma) / u_ii
        x[i] = (b[i] - suma_acumulada) / U[i, i]
        
        if contador:
            # En la fila i hay (n - 1 - i) multiplicaciones, (n - 1 - i) sumas/restas y 1 division
            contador.registrar_sustitucion(mults=operaciones_mult, 
                                          sumas=operaciones_sum, 
                                          divs=1)
            
    return x
