# Conteo Detallado de Operaciones Operacionales (Flops)

Convención declarada: 1 flop = 1 operación elemental (+, -, *, /).
Los multiplicadores se calculan una sola vez por fila y se reutilizan.
El barrido de columnas inicia en la columna k (diagonal hacia adelante).

| Paso | Filas afectadas | Mult / Div | Sum / Rest | Total Paso |
| :--- | :---: | :---: | :---: | :---: |
| Paso 1 | 11 | 154 | 143 | 297 |
| Paso 2 | 10 | 130 | 120 | 250 |
| Paso 3 | 9 | 108 | 99 | 207 |
| Paso 4 | 8 | 88 | 80 | 168 |
| Paso 5 | 7 | 70 | 63 | 133 |
| Paso 6 | 6 | 54 | 48 | 102 |
| Paso 7 | 5 | 40 | 35 | 75 |
| Paso 8 | 4 | 28 | 24 | 52 |
| Paso 9 | 3 | 18 | 15 | 33 |
| Paso 10 | 2 | 10 | 8 | 18 |
| Paso 11 | 1 | 4 | 3 | 7 |
| **Subtotal Eliminación** | - | **704** | **638** | **1342** |
| **Sustitución Regresiva** | - | **78** | **66** | **144** |
| **GRAN TOTAL FLOPs** | - | **782** | **704** | **1486** |

### Sobrecosto de Estabilidad (Pivoteo)
- Comparaciones de pivote realizadas: 66
- Intercambios de filas realizados: 8

Cota teórica de referencia: ≈ 1486 flops. Diferencia obtenida: 0 flops (0.0%).
