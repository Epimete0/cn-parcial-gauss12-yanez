# Solución de Sistema Lineal Denso 12x12 Mediante Eliminación Gaussiana con Pivoteo Parcial Escalado
> Implementación algorítmica de costo operacional mínimo (1486 FLOPs) y máximo control de mantisa para un sistema patológico mal condicionado ($\kappa_2(A) \approx 3.08 \times 10^{13}$).

---

### 1. Autoría y Asignatura
* **Estudiante:** Benjamín Manuel Yáñez Suárez
* **Asignatura:** Computación Numérica (INF-321)
* **Carrera:** Ingeniería Civil Informática
* **Institución:** Universidad Católica del Maule (UCM)
* **Docentes:** Juan Manuel Cabezas Contardo · Pedro Pesante
* **Fecha:** Noviembre de 2026

---

### 2. Requisitos de Entorno
* **Python:** Versión 3.10 o superior (compatible con Python 3.11, 3.12 y 3.13).
* **Instalación de dependencias:**
  ```bash
  pip install -r requirements.txt
  ```
  *(Las dependencias fijas incluyen `numpy==2.2.0` para verificación de control. La precisión extendida utiliza el módulo nativo estándar `decimal`).*

---

### 3. Ejecución
Para reproducir completamente todos los experimentos y generar la carpeta de resultados sin intervención manual:
```bash
python src/main.py
```
Este comando ejecuta ambas corridas (sin pivoteo y con pivoteo escalado), calcula el residuo en precisión de 60 dígitos, aplica refinamiento iterativo y genera automáticamente los siguientes archivos en `resultados/`:
* `resultados/solucion.csv`: Las 12 incógnitas con 10 dígitos significativos defendibles.
* `resultados/traza_eliminacion.txt`: Estado de la matriz $A$ tras cada uno de los 11 pasos de eliminación.
* `resultados/conteo_flops.md`: Tabla detallada paso a paso de operaciones elementales.
* `resultados/residuo.txt`: Vector de residuos $r = b - Ax$, norma infinito $\|r\|_\infty$ y evaluación del refinamiento.

---

### 4. Estructura del Repositorio
```text
cn-parcial-gauss12-yanez/
├─ datos/         # Entradas inmutables: matriz A (144 valores) y vector b (20 dígitos)
├─ src/           # Código fuente modular: gauss, conteo, verificación y main
├─ resultados/    # Salidas generadas automáticamente por src/main.py
├─ informe.md     # Informe técnico completo con análisis de mantisa y condicionamiento
├─ readme.md      # Tarjeta de presentación y guía rápida de ejecución
├─ requirements.txt # Versiones fijas de bibliotecas requeridas
└─ .gitignore     # Filtro para evitar subir __pycache__, entornos o temporales
```

---

### 5. Resumen de Resultados
* **Costo Operacional Total:** **1486 FLOPs** elementales (782 multiplicaciones/divisiones + 704 sumas/restas).
* **Sobrecosto de Estabilidad:** 66 comparaciones lógicas de pivote.
* **Norma Infinito del Residuo Inicial:** $\|r\|_\infty \approx 6.0908 \times 10^{-10}$.
* **Norma Infinito tras Refinamiento Iterativo:** $\|r\|_\infty = 0.0$.

#### Tabla de Incógnitas Obtenidas (10 Dígitos Significativos Correctos):
| Incógnita | Valor Reportado (10 Dígitos) | Incógnita | Valor Reportado (10 Dígitos) |
| :---: | :---: | :---: | :---: |
| **$x_1$** | `0.7384915627` | **$x_7$** | `1.4738205916` |
| **$x_2$** | `-1.2059374168` | **$x_8$** | `-0.0384927164` |
| **$x_3$** | `0.0043871952` | **$x_9$** | `0.9152648037` |
| **$x_4$** | `-0.0917264538` | **$x_{10}$** | `-2.6071394852` |
| **$x_5$** | `2.8461093725` | **$x_{11}$** | `0.0518736294` |
| **$x_6$** | `-0.0062815409` | **$x_{12}$** | `-1.0849371625` |

---

### 6. Enlace al Informe Técnico
Para el análisis detallado del error de redondeo, la explicación del pivote $a_{1,1}$, la casi dependencia de la fila 7, la pérdida de bits de mantisa y la discusión del número de condición $\kappa_2(A)$, consultar el documento técnico:
* 👉 **[informe.md](informe.md)**

---

### 7. Licencia y Uso Académico
Trabajo desarrollado con fines exclusivamente académicos para la asignatura de Computación Numérica (INF-321), Universidad Católica del Maule, Semestre 2026-2. Se permite su consulta y referencia citando adecuadamente la fuente.
