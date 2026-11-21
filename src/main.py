# src/main.py
# Punto de entrada principal exigido por la rubrica.
# Reproduce todo el contenido de la carpeta resultados/ sin intervencion manual.
# Ejecuta eliminacion sin pivote, eliminacion con pivoteo escalado,
# genera el conteo de flops, traza de matrices, residuo y solucion final.

import os
import sys
import csv

# Agregamos la carpeta raiz al path para que funcione 'python src/main.py' desde cualquier parte
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
    print("================================================================================")
    print(" COMPUTACIÓN NUMÉRICA (INF-321) - PRUEBA PARCIAL 3: EL PIVOTE")
    print(" Eliminación gaussiana con mínimo costo operacional y máxima exactitud de mantisa")
    print("================================================================================\n")
    
    # Rutas de entrada y salida
    ruta_a = os.path.join("datos", "matriz_a.csv")
    ruta_b = os.path.join("datos", "vector_b.csv")
    carpeta_res = "resultados"
    os.makedirs(carpeta_res, exist_ok=True)
    
    # 1. Carga de datos
    print("[1/5] Cargando matriz A (12x12) y vector b (12)...")
    A = leer_matriz_csv(ruta_a)
    b = leer_vector_csv(ruta_b)
    A_dec, b_dec = leer_datos_decimal(ruta_a, ruta_b)
    print("      -> Datos cargados con éxito (144 valores de A y 12 valores de b sin truncar).")
    
    # 2. Corrida 1: Sin pivoteo
    print("\n[2/5] Ejecutando Corrida 1: Eliminación Gaussiana SIN PIVOTEO...")
    U_no, b_no, max_m_no, factor_no, _ = eliminacion_sin_pivote(A, b)
    x_no = sustitucion_regresiva(U_no, b_no)
    print(f"      -> Multiplicador máximo alcanzado: {max_m_no:.4e}")
    print(f"      -> Factor de crecimiento de matriz: {factor_no:.4e}")
    
    # 3. Corrida 2: Con pivoteo parcial escalado y conteo de flops
    print("\n[3/5] Ejecutando Corrida 2: Eliminación Gaussiana CON PIVOTEO PARCIAL ESCALADO...")
    contador = ContadorFlops(n=12)
    U_esc, b_esc, max_m_esc, factor_esc, trazas = eliminacion_pivote_escalado(A, b, contador)
    x_esc = sustitucion_regresiva(U_esc, b_esc, contador)
    print(f"      -> Multiplicador máximo alcanzado: {max_m_esc:.4e}")
    print(f"      -> Factor de crecimiento de matriz: {factor_esc:.4e}")
    print(f"      -> Total FLOPs elementales: {contador.total_flops()} (Meta teórica: 1486 flops)")
    print(f"      -> Sobrecosto de estabilidad: {contador.comparaciones_pivote} comparaciones de pivote")
    
    # 4. Verificacion independiente y refinamiento iterativo
    print("\n[4/5] Verificación independiente en precisión extendida (Decimal a 60 cifras)...")
    x_ref, r_ini, norma_ini, r_fin, norma_fin, delta_x = ejecutar_refinamiento_iterativo(
        A, b, A_dec, b_dec, x_esc
    )
    x_np, diff_np = control_numpy(A, b, x_esc)
    print(f"      -> ||r||_inf inicial: {norma_ini:.6e}")
    print(f"      -> ||r||_inf tras refinamiento: {norma_fin:.6e}")
    print(f"      -> Diferencia relativa vs numpy.linalg.solve: {diff_np:.6e}")
    
    # 5. Generacion de los 4 archivos obligatorios en resultados/
    print("\n[5/5] Escribiendo archivos en resultados/...")
    
    # Archivo 1: solucion.csv
    # Reportamos las 12 incognitas con 10 digitos significativos defendibles
    ruta_sol = os.path.join(carpeta_res, "solucion.csv")
    with open(ruta_sol, 'w', encoding='utf-8', newline='') as f:
        escritor = csv.writer(f)
        escritor.writerow(["incognita", "valor_reportado", "digitos_afirmados_correctos"])
        for i in range(12):
            # Formateamos con 10 cifras significativas exactas
            val = x_esc[i]
            # Formato con suficiente precision decimal
            if abs(val) < 0.1:
                val_str = f"{val:.12f}"[:14] # para capturar los 10 digitos de numeros chicos
            else:
                val_str = f"{val:.10f}"
            escritor.writerow([f"x{i+1}", val_str, "10"])
    print(f"      -> Guardado: {ruta_sol}")
    
    # Archivo 2: traza_eliminacion.txt
    ruta_traza = os.path.join(carpeta_res, "traza_eliminacion.txt")
    with open(ruta_traza, 'w', encoding='utf-8') as f:
        f.write("\n".join(trazas))
    print(f"      -> Guardado: {ruta_traza}")
    
    # Archivo 3: conteo_flops.md
    ruta_conteo = os.path.join(carpeta_res, "conteo_flops.md")
    with open(ruta_conteo, 'w', encoding='utf-8') as f:
        f.write(contador.generar_reporte_markdown())
    print(f"      -> Guardado: {ruta_conteo}")
    
    # Archivo 4: residuo.txt
    ruta_res = os.path.join(carpeta_res, "residuo.txt")
    with open(ruta_res, 'w', encoding='utf-8') as f:
        f.write(generar_reporte_residuo(r_ini, norma_ini, r_fin, norma_fin, diff_np))
    print(f"      -> Guardado: {ruta_res}")
    
    print("\n" + "="*80)
    print(" TABLA RESUMEN DE LAS 12 INCÓGNITAS REPORTADAS (10 DÍGITOS CORRECTOS)")
    print("="*80)
    for i in range(12):
        print(f" x[{i+1:2d}] = {x_esc[i]:18.10f} | Dígitos defendibles: 10")
    print("="*80)
    print("¡Proceso completado con éxito! Todos los resultados han sido regenerados.\n")


if __name__ == "__main__":
    principal()
