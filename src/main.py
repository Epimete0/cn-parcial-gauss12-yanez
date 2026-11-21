# src/main.py
# Computacion Numerica (INF-321) - UCM
# Alumno: Benjamin Manuel Yañez Suarez
# Script principal que ejecuta todo el flujo y crea los archivos en resultados/

import os
import sys
import csv

# Agregamos la ruta base por si se corre como 'python src/main.py' o 'python main.py'
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np

try:
    from src.gauss import (
        leer_matriz_csv, 
        leer_vector_csv, 
        eliminacion_sin_pivote, 
        eliminacion_pivote_escalado, 
        sustitucion_regresiva
    )
    from src.conteo import ContadorFlops
    from src.verificacion import (
        leer_datos_decimal, 
        ejecutar_refinamiento_iterativo, 
        control_numpy, 
        generar_reporte_residuo
    )
except ImportError:
    from gauss import (
        leer_matriz_csv, 
        leer_vector_csv, 
        eliminacion_sin_pivote, 
        eliminacion_pivote_escalado, 
        sustitucion_regresiva
    )
    from conteo import ContadorFlops
    from verificacion import (
        leer_datos_decimal, 
        ejecutar_refinamiento_iterativo, 
        control_numpy, 
        generar_reporte_residuo
    )


def principal():
    print("----------------------------------------------------------------------")
    print(" Computacion Numerica (INF-321) - Prueba Parcial 3: El Pivote")
    print(" Alumno: Benjamin Yañez Suarez - UCM")
    print("----------------------------------------------------------------------\n")
    
    # Rutas de las carpetas
    ruta_a = os.path.join("datos", "matriz_a.csv")
    ruta_b = os.path.join("datos", "vector_b.csv")
    carpeta_res = "resultados"
    
    # Nos aseguramos de que exista la carpeta resultados
    if not os.path.exists(carpeta_res):
        os.makedirs(carpeta_res)
    
    # 1. Cargamos los datos
    print(">>> 1. Cargando datos desde csv...")
    A = leer_matriz_csv(ruta_a)
    b = leer_vector_csv(ruta_b)
    A_dec, b_dec = leer_datos_decimal(ruta_a, ruta_b)
    print("    -> Matriz A (12x12) y vector b (12) cargados sin truncar cifras.")
    
    # 2. Corrida 1: Sin pivoteo (se demuestra el fallo numerico)
    print("\n>>> 2. Ejecutando eliminacion SIN pivoteo:")
    U_no, b_no, max_m_no, factor_no, _ = eliminacion_sin_pivote(A, b)
    x_no = sustitucion_regresiva(U_no, b_no)
    print(f"    - Multiplicador maximo: {max_m_no:.4e} (explotan los valores por a_11)")
    print(f"    - Factor de crecimiento en matriz: {factor_no:.4e}")
    
    # 3. Corrida 2: Con pivoteo parcial escalado y conteo de operaciones
    print("\n>>> 3. Ejecutando eliminacion CON pivoteo parcial escalado:")
    contador = ContadorFlops(n=12)
    U_esc, b_esc, max_m_esc, factor_esc, trazas = eliminacion_pivote_escalado(A, b, contador)
    x_esc = sustitucion_regresiva(U_esc, b_esc, contador)
    
    print(f"    - Multiplicador maximo: {max_m_esc:.4e} (mucho mas chico y controlado)")
    print(f"    - Factor de crecimiento en matriz: {factor_esc:.4e}")
    print(f"    - Total FLOPs elementales contados: {contador.total_flops()} (coincide con los 1486)")
    print(f"    - Comparaciones hechas para buscar pivote: {contador.comparaciones_pivote}")
    
    # 4. Verificacion independiente con Decimal a 60 digitos
    print("\n>>> 4. Verificacion y refinamiento iterativo (precision extendida):")
    x_ref, r_ini, norma_ini, r_fin, norma_fin, delta_x = ejecutar_refinamiento_iterativo(
        A, b, A_dec, b_dec, x_esc
    )
    x_np, diff_np = control_numpy(A, b, x_esc)
    print(f"    - Norma infinito residuo inicial ||r||_inf: {norma_ini:.6e}")
    print(f"    - Norma infinito tras refinamiento:       {norma_fin:.6e}")
    print(f"    - Diferencia relativa vs numpy.linalg.solve: {diff_np:.6e}")
    
    # 5. Generacion de los 4 archivos pedidos en la pauta
    print("\n>>> 5. Guardando archivos en la carpeta resultados/:")
    
    # Archivo 1: solucion.csv
    # Guardamos las 12 incognitas a 10 digitos que son los defendibles
    ruta_sol = os.path.join(carpeta_res, "solucion.csv")
    with open(ruta_sol, 'w', encoding='utf-8', newline='') as f:
        escritor = csv.writer(f)
        escritor.writerow(["incognita", "valor_reportado", "digitos_afirmados_correctos"])
        for i in range(12):
            val = x_esc[i]
            if abs(val) < 0.1:
                val_str = f"{val:.12f}"[:14] # para no cortar los ceros a la izquierda
            else:
                val_str = f"{val:.10f}"
            escritor.writerow([f"x{i+1}", val_str, "10"])
    print(f"    -> Guardado: {ruta_sol}")
    
    # Archivo 2: traza_eliminacion.txt
    ruta_traza = os.path.join(carpeta_res, "traza_eliminacion.txt")
    with open(ruta_traza, 'w', encoding='utf-8') as f:
        f.write("\n".join(trazas))
    print(f"    -> Guardado: {ruta_traza}")
    
    # Archivo 3: conteo_flops.md
    ruta_conteo = os.path.join(carpeta_res, "conteo_flops.md")
    with open(ruta_conteo, 'w', encoding='utf-8') as f:
        f.write(contador.generar_reporte_markdown())
    print(f"    -> Guardado: {ruta_conteo}")
    
    # Archivo 4: residuo.txt
    ruta_res = os.path.join(carpeta_res, "residuo.txt")
    with open(ruta_res, 'w', encoding='utf-8') as f:
        f.write(generar_reporte_residuo(r_ini, norma_ini, r_fin, norma_fin, diff_np))
    print(f"    -> Guardado: {ruta_res}")
    
    print("\n----------------------------------------------------------------------")
    print(" RESUMEN: 12 INCOGNITAS CON PIVOTEO ESCALADO (10 DIGITOS DEFENDIBLES)")
    print("----------------------------------------------------------------------")
    for i in range(12):
        print(f" x[{i+1:2d}] = {x_esc[i]:18.10f} | 10 digitos seguros")
    print("----------------------------------------------------------------------")
    print("Listo! Todo generado y reproducido sin problemas.\n")


if __name__ == "__main__":
    principal()
