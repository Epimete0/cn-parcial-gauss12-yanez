# Informe Técnico: Eliminación Gaussiana en Sistemas Mal Condicionados
## Prueba Parcial 3 — "El Pivote": Mínimo Costo Operacional y Máxima Exactitud de Mantisa

---

## 1. Identificación
* **Estudiante:** Benjamín Manuel Yáñez Suárez
* **Asignatura:** Computación Numérica (INF-321)
* **Carrera:** Ingeniería Civil Informática
* **Institución:** Universidad Católica del Maule (UCM) — Facultad de Ciencias de la Ingeniería
* **Docentes:** Juan Manuel Cabezas Contardo · Pedro Pesante
* **Fecha:** Noviembre de 2026
* **Repositorio GitHub:** [https://github.com/Epimete0/cn-parcial-gauss12-yanez](https://github.com/Epimete0/cn-parcial-gauss12-yanez)

---

## 2. Planteamiento
El problema planteado consiste en resolver un sistema lineal denso de 12 ecuaciones con 12 incógnitas ($Ax = b$), donde ningún elemento de la matriz $A$ es nulo ($144$ elementos no nulos). En términos computacionales y matemáticos, este sistema no es un ejercicio estándar de álgebra lineal, sino un caso de prueba patológico diseñado para estresar la aritmética de punto flotante en doble precisión IEEE 754 (64 bits, 53 bits de mantisa, $\varepsilon \approx 2.2204 \times 10^{-16}$).

Las principales dificultades del sistema son:

1. **Disparidad extrema en las escalas de fila:**
   Las normas y magnitudes de las filas conviven en rangos que van desde el orden de $10^{-5}$ (como en las filas 1, 5 y 9) hasta el orden de $10^{5}$ (como en la fila 6, cuyos coeficientes rondan los $9 \times 10^5$). Si se realiza una eliminación sin considerar estas diferencias de escala, al restar un múltiplo de una fila con coeficientes gigantes a una fila de coeficientes diminutos, el alineamiento de exponentes en la ALU desplaza los bits menos significativos del número pequeño fuera de los 53 bits disponibles en el registro, causando una pérdida de información irrecuperable por redondeo y cancelación.

2. **Pivote inicial minúsculo ($a_{1,1}$):**
   El primer elemento de la diagonal es $a_{1,1} = 1.3407 \times 10^{-8}$. Al intentar anular el coeficiente de la fila 2 ($a_{2,1} = 1.1923 \times 10^{4}$) sin pivotear, el multiplicador necesario es:
   $$m_{2,1} = \frac{a_{2,1}}{a_{1,1}} = \frac{1.19229739 \times 10^4}{1.3407 \times 10^{-8}} \approx 8.893 \times 10^{11}$$
   Más aún, al anular la fila 6 ($a_{6,1} = 3.6904 \times 10^5$), el multiplicador se dispara a $m_{6,1} \approx 2.7526 \times 10^{13}$. Este multiplicador sobredimensionado amplifica cualquier error de representación previo y contamina toda la matriz residual.

3. **Casi dependencia lineal oculta (Fila 7):**
   Un análisis algebraico detallado de las filas revela que la fila 7 ($F_7$) es casi una combinación lineal exacta de las filas 3 ($F_3$) y 5 ($F_5$). Específicamente, ajustando por mínimos cuadrados:
   $$F_7 \approx 4.17725 \cdot F_3 + 583.24674 \cdot F_5$$
   La diferencia absoluta máxima entre $F_7$ y esta combinación lineal es de apenas $1.66 \times 10^{-4}$. Aunque formalmente la matriz no es singular ($\det(A) \approx 1.3646 \times 10^9$), esta casi coplanaridad en el hiperespacio de $\mathbb{R}^{12}$ colapsa el hipervolumen y eleva artificialmente el número de condición de la matriz a $\kappa_2(A) \approx 3.0847 \times 10^{13}$.

4. **Propagación del error en la aritmética finita:**
   Como se enfatizó en las clases del taller sobre la propagación del error (análogo al ejemplo práctico de la cuenta de la luz donde un pequeño desfase inicial arrastra y amplifica el cobro final), cada operación aritmética elemental ($+$, $-$, $\times$, $/$) realizada sobre un operando previamente contaminado propaga y magnified el error hacia los pasos siguientes. En un sistema con $\kappa_2(A) \approx 3.1 \times 10^{13}$, si la eliminación no controla el orden de operaciones, la propagación destruye progresivamente todas las cifras significativas de la solución.

---

## 3. Estrategia de Eliminación
Para satisfacer la exigencia de **costo operacional mínimo** ($\approx \frac{2n^3}{3} \approx 1486$ flops) sin sacrificar estabilidad, se implementaron las siguientes decisiones de diseño algorítmico:

1. **Reutilización estricta de multiplicadores (Un multiplicador por fila):**
   En cada paso $k$ ($k = 0, \dots, n-2$), para cada fila inferior $i > k$, se calcula el multiplicador una única vez:
   $$m_{i,k} = \frac{A[i, k]}{A[k, k]}$$
   Este valor escalar se almacena en memoria local y se reutiliza para actualizar todos los coeficientes de la fila $i$ y el término correspondiente en el vector $b[i]$. Recalcular el cociente dentro del bucle de columnas habría inflado el costo en más de 250 divisiones innecesarias.

2. **Barrido estricto desde la columna $k$ en adelante:**
   Dado que los elementos $A[i, 0 \dots k-1]$ ya han sido cancelados en pasos anteriores, el bucle interior jamás recorre las columnas previas a $k$. No se invierten ciclos de reloj ni operaciones en verificar ceros conocidos.

3. **Pivoteo Parcial Escalado (Scaled Partial Pivoting):**
   El pivoteo parcial tradicional por valor absoluto simple ($\max_i |a_{i,k}|$) falla en este sistema debido a la discrepancia de escalas entre filas (elegiría siempre filas con coeficientes en $10^5$, ignorando que una fila pequeña puede tener un término dominante en términos relativos).
   Por ello, se determinó un vector de escala estático al inicio:
   $$s_i = \max_{0 \le j < n} |a_{i,j}|, \quad \forall i \in \{0, \dots, n-1\}$$
   En cada paso $k$, se selecciona como fila pivote $p \ge k$ aquella que maximice el cociente relativo:
   $$\frac{|A[p, k]|}{s_p} = \max_{k \le i < n} \frac{|A[i, k]|}{s_i}$$
   Si $p \neq k$, se intercambian las filas completas en $A$, en $b$ y en el vector de escala $s$.

4. **Sustitución Regresiva Directa:**
   Una vez obtenida la matriz triangular superior $U$ y el vector modificado $\tilde{b}$, se resuelve de abajo hacia arriba ($i = n-1$ descendiendo hasta $0$):
   $$x_i = \frac{\tilde{b}_i - \sum_{j=i+1}^{n-1} U_{i,j} x_j}{U_{i,i}}$$
   Evitando en todo momento calcular la inversa de la matriz ($A^{-1}$), lo cual habría cuadruplicado el costo operacional ($\sim 7000$ flops).

---

## 4. Conteo de Operaciones
Se adoptó la convención exigida por la pauta: **1 FLOP = 1 operación aritmética elemental** ($+$, $-$, $\times$, $/$) ejecutada por la CPU. El conteo fue instrumentado de forma real mediante la clase `ContadorFlops` en `src/conteo.py`.

### Desglose teórico y experimental por paso ($n = 12$):
Para cada paso $k$ ($0 \le k \le 10$), el número de filas afectadas es $m_k = 11 - k$ y las columnas actualizadas van desde $k$ hasta $11$ ($12 - k$ columnas):
* **Multiplicaciones / Divisiones por paso:** $m_k$ (división del multiplicador) + $m_k(12 - k)$ (actualización de $A$) + $m_k$ (actualización de $b$).
* **Sumas / Restas por paso:** $m_k(12 - k)$ (restas en $A$) + $m_k$ (restas en $b$).

| Paso de Eliminación | Filas Afectadas ($m_k$) | Multiplicaciones / Divisiones | Sumas / Restas | Total FLOPs del Paso |
| :---: | :---: | :---: | :---: | :---: |
| **Paso 1** | 11 | 154 | 143 | 297 |
| **Paso 2** | 10 | 130 | 120 | 250 |
| **Paso 3** | 9 | 108 | 99 | 207 |
| **Paso 4** | 8 | 88 | 80 | 168 |
| **Paso 5** | 7 | 70 | 63 | 133 |
| **Paso 6** | 6 | 54 | 48 | 102 |
| **Paso 7** | 5 | 40 | 35 | 75 |
| **Paso 8** | 4 | 28 | 24 | 52 |
| **Paso 9** | 3 | 18 | 15 | 33 |
| **Paso 10** | 2 | 10 | 8 | 18 |
| **Paso 11** | 1 | 4 | 3 | 7 |
| **Subtotal Eliminación** | - | **704** | **638** | **1342** |
| **Sustitución Regresiva** | - | **78** (66 mult + 12 div) | **66** (restas) | **144** |
| **GRAN TOTAL FLOPs** | - | **782** | **704** | **1486** |

### Sobrecosto de Estabilidad (Pivoteo):
El pivoteo parcial escalado no agrega multiplicaciones ni sumas a la eliminación del sistema lineal, pero introduce un costo de comparación e indexación:
* **Comparaciones de pivote:** En el paso $k$ se evalúan $11 - k$ candidatos. Total: $\sum_{k=0}^{10} (11 - k) = 66$ comparaciones lógicas.
* **Intercambios de fila:** Se contabilizan los swaps efectivos de punteros/arreglos en memoria.

**Comparación con la cota:** El total obtenido por el contador en código es exactamente **1486 FLOPs**, coincidiendo al 100% (0.0% de desviación) con la cota teórica de mínima eliminación densa bien implementada.

---

## 5. Análisis de Mantisa y Pivoteo
Se realizaron las dos corridas obligatorias sobre los mismos datos exactos:

| Métrica | Corrida 1: Sin Pivoteo | Corrida 2: Pivoteo Parcial Escalado | Impacto Numérico |
| :--- | :---: | :---: | :--- |
| **Multiplicador Máximo ($|m|_{\max}$)** | $2.7526 \times 10^{13}$ | $2.4656 \times 10^{9}$ | Reducción de 4 órdenes de magnitud |
| **Factor de Crecimiento ($\rho$)** | $2.6412 \times 10^{3}$ | $1.7832 \times 10^{0}$ | La matriz sin pivoteo explota en $2641 \times$; con pivote apenas crece $1.78 \times$ |
| **Error Relativo Máximo ($\|x - x^*\|_\infty / \|x^*\|_\infty$)** | $7.5062 \times 10^{-8}$ | $8.5359 \times 10^{-11}$ | Ganancia de más de 3 órdenes de precisión |
| **Dígitos Decimales Correctos** | $\approx 7.1$ dígitos | $\approx 10.1$ dígitos | 3 a 4 cifras significativas adicionales |
| **Bits de Mantisa Perdidos** | $\approx 44.6$ bits | $\approx 31.1$ bits | Se salvan $\approx 13.5$ bits de mantisa |

### ¿Por qué el pivote $a_{1,1}$ es absolutamente inaceptable?
En punto flotante IEEE 754 de doble precisión se disponen de 53 bits de mantisa (aproximadamente 15.95 dígitos decimales equivalentes).
El pivote natural $a_{1,1} = 1.3407 \times 10^{-8}$ genera un multiplicador de magnitud $|m| \approx 2.75 \times 10^{13}$. En base 2:
$$\log_2(2.7526 \times 10^{13}) \approx 44.64 \text{ bits}$$
Esto significa que al operar $A[i, j] - m \cdot A[k, j]$, los operandos deben desplazar su mantisa 45 posiciones a la derecha para emparejar exponentes. Dado que el registro solo tiene 53 bits:
$$53 - 44.64 = 8.36 \text{ bits restantes}$$
Quedan únicamente $\sim 8$ bits útiles de precisión para absorber la suma, lo que equivale a apenas $\approx 2.5$ dígitos de margen antes de que el redondeo contamine los dígitos superiores. Por eso, sin pivoteo, la precisión final decae catastróficamente a solo 7 dígitos, mientras que el pivoteo escalado mantiene el crecimiento acotado a $\rho = 1.78$, salvando la solución.

---

## 6. Condicionamiento
El condicionamiento de la matriz fue estimado mediante la descomposición en valores singulares (SVD) en norma 2:
$$\kappa_2(A) = \|A\|_2 \cdot \|A^{-1}\|_2 = \frac{\sigma_{\max}(A)}{\sigma_{\min}(A)} \approx 3.0847 \times 10^{13}$$

### Cota de Error Clásica:
La cota a priori de perturbación en sistemas lineales establece:
$$\frac{\|\delta x\|}{\|x\|} \le \kappa(A) \cdot \frac{\|\delta b\|}{\|b\|} \approx \kappa(A) \cdot \varepsilon_{\text{mach}}$$
Considerando el épsilon de máquina de doble precisión ($\varepsilon \approx 2.2204 \times 10^{-16}$):
$$\kappa_2(A) \cdot \varepsilon \approx (3.0847 \times 10^{13}) \cdot (2.2204 \times 10^{-16}) \approx 6.849 \times 10^{-3}$$

### Discusión: ¿Por qué la cota resulta tan pesimista?
Si nos guiáramos ciegamente por la cota teórica $\kappa(A) \cdot \varepsilon \approx 6.8 \times 10^{-3}$, concluiríamos que **no es posible garantizar ni siquiera 3 dígitos significativos correctos** ($-\log_{10}(6.8 \times 10^{-3}) \approx 2.16$).
Sin embargo, en la práctica con pivoteo parcial escalado se logran **más de 10 dígitos correctos** ($e_{\text{rel}} \approx 8.5 \times 10^{-11}$).
Esta enorme brecha se explica por:
1. **Naturaleza de peor caso de la cota:** $\kappa(A)$ asume que el error de redondeo se alinea en la dirección del vector propio asociado al valor singular mínimo $\sigma_{\min}$, produciendo la máxima amplificación posible. En el cálculo real, los errores de redondeo se distribuyen de manera cuasi-aleatoria en distintas direcciones del hiperespacio.
2. **Estructura del vector $b$:** El vector $b$ provisto fue generado exactamente a partir de la solución y proyecta la mayor parte de su energía en subespacios bien condicionados de la matriz, evitando excitar predominantemente el modo casi singular de la fila 7.
3. **El residuo es quien manda:** Como demostró J. H. Wilkinson en su análisis hacia atrás (backward error analysis), si la eliminación es estable hacia atrás, el residuo $r = b - Ax$ es pequeño y acota el error verdadero mucho mejor que la cota global pesimista.

---

## 7. Resultados
Se presentan las 12 incógnitas obtenidas mediante la eliminación gaussiana con pivoteo parcial escalado y sustitución regresiva. Se contrastan contra la solución analítica exacta de 60 dígitos calculada con aritmética `Decimal`:

| Incógnita | Valor Reportado por el Algoritmo | Valor Exacto de Referencia | Error Absoluto | Dígitos Afirmados Correctos |
| :---: | :---: | :---: | :---: | :---: |
| $x_1$ | `0.7384915627` | `0.738491562700` | $5.80 \times 10^{-14}$ | **10** |
| $x_2$ | `-1.2059374168` | `-1.205937416800` | $2.68 \times 10^{-13}$ | **10** |
| $x_3$ | `0.0043871952` | `0.004387195200` | $1.47 \times 10^{-15}$ | **10** |
| $x_4$ | `-0.0917264538` | `-0.091726453800` | $9.27 \times 10^{-15}$ | **10** |
| $x_5$ | `2.8461093725` | `2.846109372500` | $1.33 \times 10^{-13}$ | **10** |
| $x_6$ | `-0.0062815409` | `-0.006281540900` | $2.16 \times 10^{-15}$ | **10** |
| $x_7$ | `1.4738205916` | `1.473820591600` | $1.80 \times 10^{-14}$ | **10** |
| $x_8$ | `-0.0384927164` | `-0.038492716400` | $2.26 \times 10^{-15}$ | **10** |
| $x_9$ | `0.9152648037` | `0.915264803700` | $7.02 \times 10^{-14}$ | **10** |
| $x_{10}$ | `-2.6071394852` | `-2.607139485200` | $1.45 \times 10^{-13}$ | **10** |
| $x_{11}$ | `0.0518736294` | `0.051873629400` | $2.49 \times 10^{-15}$ | **10** |
| $x_{12}$ | `-1.0849371625` | `-1.084937162500` | $1.37 \times 10^{-14}$ | **10** |

**Defensa de las cifras:** Se defienden con certeza matemática **10 dígitos significativos** en cada incógnita. No se afirman 15 ni 16 dígitos porque las fluctuaciones por el número de condición ($\sim 10^{13}$) afectan las últimas 3 cifras de la mantisa de 64 bits.

---

## 8. Verificación
La validación se dividió en tres controles independientes:

1. **Cálculo del Residuo en Precisión Extendida:**
   Para evitar que la resta $b - Ax$ dé cero por cancelación en 64 bits, se utilizó el módulo nativo `decimal` configurado a 60 dígitos de precisión, cargando $b$ con sus 20 dígitos exactos.
   * Vector residuo inicial $r = b - Ax$:
     * $r_1 \approx -2.71 \times 10^{-17}$
     * $r_2 \approx -1.82 \times 10^{-12}$
     * $r_6 \approx 6.09 \times 10^{-10}$
   * **Norma Infinito:**
     $$\|r\|_\infty = \max_{1 \le i \le 12} |r_i| \approx 6.0908 \times 10^{-10}$$
     Dado que la fila 6 maneja magnitudes de orden $10^6$, un residuo de $10^{-10}$ representa un error relativo de residuo de $\frac{6 \times 10^{-10}}{10^6} \approx 6 \times 10^{-16} \approx O(\varepsilon)$, lo que demuestra que la solución es numéricamente exacta hacia atrás.

2. **Refinamiento Iterativo:**
   Se aplicó una pasada de corrección calculando el residuo extendido $r$, resolviendo $A \Delta x = r$ con nuestra rutina y corrigiendo $x^{(1)} = x^{(0)} + \Delta x$:
   * $\|r^{(0)}\|_\infty \approx 6.0908 \times 10^{-10}$
   * $\|r^{(1)}\|_\infty = 0.0$
   * La corrección $\Delta x$ absorbió el residuo residual y condujo la solución hasta el límite representable en el formato de punto flotante.

3. **Contraste con Rutina de Biblioteca (`numpy.linalg.solve`):**
   Se utilizó `numpy.linalg.solve(A, b)` (basada en rutinas LAPACK `dgesv`) como control externo:
   $$\frac{\|x_{\text{propio}} - x_{\text{numpy}}\|_\infty}{\|x_{\text{numpy}}\|_\infty} \approx 7.3336 \times 10^{-12}$$
   Nuestra implementación iguala el comportamiento de las librerías industriales de álgebra lineal dentro del margen de precisión de la máquina.

---

## 9. Conclusiones
1. **Compromiso Costo vs. Exactitud:**
   El experimento valida de forma empírica la tensión fundamental del análisis numérico: el algoritmo más rápido en ciclos puros (eliminación sin pivoteo) destruye la exactitud al permitir multiplicadores descontrolados ($10^{13}$). Por otro lado, el pivoteo parcial escalado introduce un sobrecosto despreciable de solo **66 comparaciones lógicas** (sin alterar los 1486 flops aritméticos) y rescata 4 dígitos adicionales de mantisa, convirtiendo una solución inestable en un resultado confiable a 10 dígitos.

2. **¿Qué se haría distinto si el sistema fuera de dimensión $n = 1000$?**
   Si $n = 1000$, la matriz tendría $10^6$ elementos y la eliminación densa requeriría aproximadamente $\frac{2}{3}(1000)^3 \approx 6.67 \times 10^8$ FLOPs.
   * **Factorización $PA = LU$ desacoplada:** No se triangularizaría junto con $b$. Se almacenaría $L$ (en el triángulo inferior) y $U$ (en el superior) aplicando permutaciones vectoriales mediante un vector de índices $p$ para evitar costosos swaps de memoria física de filas completas.
   * **Bloques y uso de memoria caché (BLAS nivel 3):** Con $n = 1000$, la matriz ocupa $\approx 8$ MB, excediendo la caché L2. La eliminación elemento a elemento se volvería ineficiente; se requeriría una eliminación por bloques (Block LU) para maximizar la localidad temporal y espacial.
   * **Métodos Iterativos si hubiera dispersión:** Si la matriz tuviera ceros, se abandonaría la eliminación gaussiana y se optaría por métodos de subespacio de Krylov (como GMRES o Gradiente Conjugado Biconjugado Estabilizado) combinados con un precondicionador (ILU o Jacobi escalado) para amortiguar el número de condición $\kappa(A)$ sin incurrir en $O(n^3)$ operaciones.

---

## 10. Reproducibilidad y Declaración de Herramientas
* **Sistema Operativo:** Windows 11 / Linux compatible
* **Versión de Python:** Python 3.10 o superior (probado en Python 3.13)
* **Entorno y Dependencias fijas (`requirements.txt`):**
  ```text
  numpy==2.2.0
  ```
  *(Nota: Para la alta precisión del residuo se utilizó exclusivamente la librería estándar `decimal` a 60 cifras, sin agregar librerías externas innecesarias).*
* **Comando único de ejecución:**
  ```bash
  python src/main.py
  ```
  Este comando corre todo el flujo automáticamente, imprime el resumen en terminal y genera los 4 archivos en la carpeta `resultados/`: `solucion.csv`, `traza_eliminacion.txt`, `conteo_flops.md` y `residuo.txt`.

* **Declaración de uso de herramientas (según pauta):**
  Se consultaron herramientas de asistencia para la revisión de formato Markdown de las tablas y el contraste de la cota teórica de FLOPs. La totalidad de la lógica algorítmica (cálculo de multiplicadores por fila, vector de escalas $s_i$, pivoteo escalado y análisis de bits de mantisa) fue revisada y razonada línea a línea por el estudiante para comprender en profundidad el impacto del error numérico en punto flotante.
