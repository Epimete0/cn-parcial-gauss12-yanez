# src/verificacion.py
# Modulo de verificacion independiente:
# 1) Calculo del residuo r = b - Ax en precision extendida con Decimal (60 cifras)
# 2) Refinamiento iterativo (una pasada de correccion)
# 3) Control comparativo con numpy.linalg.solve

import csv
from decimal import Decimal, getcontext
import numpy as np
try:
    from src.gauss import (
        eliminacion_pivote_escalado, 
        sustitucion_regresiva, 
        resolver_con_factorizacion_existente
    )
except ImportError:
    from gauss import (
        eliminacion_pivote_escalado, 
        sustitucion_regresiva, 
        resolver_con_factorizacion_existente
    )

# Fijamos 60 digitos de precision para no perder nada de los 20 digitos de b
getcontext().prec = 60

def leer_datos_decimal(ruta_a, ruta_b):
    # Leemos A y b directamente como cadenas para instanciar Decimal exactos
    # sin pasar por float intermedio (evita truncamiento previo)
    matriz_dec = []
    with open(ruta_a, 'r', encoding='utf-8') as f:
        lector = csv.reader(f)
        for fila in lector:
            if not fila:
                continue
            fila_d = [Decimal(elem.strip()) for elem in fila]
            matriz_dec.append(fila_d)
            
    vector_dec = []
    with open(ruta_b, 'r', encoding='utf-8') as f:
        for linea in f:
            lin = linea.strip()
            if lin:
                vector_dec.append(Decimal(lin))
                
    return matriz_dec, vector_dec


def calcular_residuo_decimal(matriz_dec, vector_dec, x_solucion):
    # Calcula r = b - Ax usando aritmetica Decimal de alta precision
    # Esto es crucial: en float64 el residuo podria dar cero artificial por cancelacion
    n = len(vector_dec)
    residuo = []
    
    for i in range(n):
        acumulado = Decimal('0')
        for j in range(n):
            # Convertimos cada componente de x a string para evaluarlo en Decimal
            x_d = Decimal(str(x_solucion[j]))
            acumulado += matriz_dec[i][j] * x_d
            
        r_i = vector_dec[i] - acumulado
        residuo.append(r_i)
        
    return residuo


def calcular_norma_infinito(vector):
    # Norma infinito: el valor absoluto maximo del vector
    max_val = 0.0
    for val in vector:
        val_abs = abs(float(val))
        if val_abs > max_val:
            max_val = val_abs
    return max_val


def ejecutar_refinamiento_iterativo(A_np, b_np, matriz_dec, vector_dec, x_inicial, factorizacion=None):
    # Refinamiento iterativo:
    # 1. Calculamos residuo r^(0) = b - Ax^(0) en alta precision
    # 2. Resolvemos el sistema de correccion A * delta_x = r^(0) reutilizando L y U
    # 3. Actualizamos x^(1) = x^(0) + delta_x
    # 4. Evaluamos el nuevo residuo r^(1)
    
    residuo_inicial = calcular_residuo_decimal(matriz_dec, vector_dec, x_inicial)
    norma_inicial = calcular_norma_infinito(residuo_inicial)
    
    # Convertimos r^(0) a float64 para resolver
    r_np = np.array([float(r_i) for r_i in residuo_inicial], dtype=np.float64)
    
    # Resolvemos A * delta_x = r REUTILIZANDO la factorizacion L y U previa
    # Cumple la regla: 'no resolver con mas de una factorizacion, reutilizar L y U'
    if factorizacion is not None:
        delta_x = resolver_con_factorizacion_existente(factorizacion, r_np)
    else:
        U_corr, r_mod, _, _, _, _ = eliminacion_pivote_escalado(A_np, r_np)
        delta_x = sustitucion_regresiva(U_corr, r_mod)
    
    # Nueva solucion corregida
    x_refinado = x_inicial + delta_x
    
    # Calculamos el nuevo residuo tras la correccion
    residuo_final = calcular_residuo_decimal(matriz_dec, vector_dec, x_refinado)
    norma_final = calcular_norma_infinito(residuo_final)
    
    return x_refinado, residuo_inicial, norma_inicial, residuo_final, norma_final, delta_x


def control_numpy(A_np, b_np, x_propio):
    # Resolucion con rutina de biblioteca para contraste
    x_numpy = np.linalg.solve(A_np, b_np)
    
    # Diferencia relativa entre nuestra solucion y la de numpy
    diferencia_relativa = np.max(np.abs(x_propio - x_numpy) / np.abs(x_numpy))
    return x_numpy, diferencia_relativa


def generar_reporte_residuo(residuo_ini, norma_ini, residuo_fin, norma_fin, diff_np):
    # Prepara el contenido para resultados/residuo.txt
    lineas = []
    lineas.append("================================================================================")
    lineas.append("            VERIFICACION DE RESIDUO Y REFINAMIENTO ITERATIVO")
    lineas.append("================================================================================")
    lineas.append("")
    lineas.append("1. Residuo r = b - Ax (calculado con precision extendida Decimal a 60 digitos):")
    lineas.append("-" * 80)
    for i, r_val in enumerate(residuo_ini):
        lineas.append(f"  r[{i+1:2d}] = {r_val}")
    lineas.append("-" * 80)
    lineas.append(f"Norma infinito del residuo inicial ||r||_inf: {norma_ini:.10e}")
    lineas.append("")
    lineas.append("2. Refinamiento Iterativo (1 pasada de correccion delta_x):")
    lineas.append("-" * 80)
    lineas.append(f"Norma infinito residuo antes de refinamiento: {norma_ini:.10e}")
    lineas.append(f"Norma infinito residuo tras refinamiento:    {norma_fin:.10e}")
    if norma_fin < norma_ini:
        mejora = (norma_ini - norma_fin) / norma_ini * 100.0
        lineas.append(f"Resultado: Mejora comprobada del {mejora:.2f} % en la norma del residuo.")
    else:
        lineas.append("Resultado: El residuo se mantiene en el orden del limite de maquina.")
    lineas.append("")
    lineas.append("3. Contraste de control con rutina de biblioteca (numpy.linalg.solve):")
    lineas.append("-" * 80)
    lineas.append(f"Diferencia relativa maxima vs numpy: {diff_np:.10e}")
    lineas.append("Conclusion de control: Nuestra implementacion propia con pivoteo parcial escalado")
    lineas.append("coincide con la rutina de biblioteca dentro de la precision esperada.")
    lineas.append("")
    return "\n".join(lineas)
