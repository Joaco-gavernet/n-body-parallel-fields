import os
import re
import numpy as np
import argparse

# Argumento para cantidad de cuerpos
parser = argparse.ArgumentParser(description="Comparar resultados N-Body para un valor de N específico")
parser.add_argument("N", type=int, help="Cantidad de cuerpos a comparar")
args = parser.parse_args()
N_OBJETIVO = args.N

# Detectar archivos .log que empiezan con n_body y terminan en .log
todos_los_logs = [f for f in os.listdir('.') if f.startswith('n_body') and f.endswith('.log')]

# Diccionario para mantener orden de prioridad
orden_claves = []
archivos = {}
tiempos = {}

for archivo in sorted(todos_los_logs):
    if match := re.match(r'n_body_(\d+)\.log$', archivo):
        N = int(match.group(1))
        if N != N_OBJETIVO:
            continue
        clave = f'secuencial-{N:06d}'
    elif match := re.match(r'n_body_pthread_(\d+)_(\d+)\.log$', archivo):
        T, N = map(int, match.groups())
        if N != N_OBJETIVO:
            continue
        clave = f'pthreads-{T:03d}-{N:06d}'
    elif match := re.match(r'n_body_hybrid_(\d+)_(\d+)_(\d+)\.log$', archivo):
        P, T, N = map(int, match.groups())
        if N != N_OBJETIVO:
            continue
        clave = f'hybrid-{P:03d}-{T:03d}-{N:06d}'
    else:
        continue

    archivos[clave] = archivo
    orden_claves.append(clave)

    with open(archivo, 'r') as f:
        for linea in f:
            if 'Tiempo en segundos' in linea:
                match = re.search(r"([\d.]+)", linea)
                if match:
                    tiempo = np.float128(match.group(1))
                    tiempos[clave] = tiempo
                break

orden_claves = list(dict.fromkeys(orden_claves))


def cargar_posiciones(path):
    cuerpos = {}
    patron = re.compile(
        r"Cuerpo (\d+): (?:\(([-+eE0-9.]+), ([-+eE0-9.]+), ([-+eE0-9.]+)\)|px=([-+eE0-9.]+), py=([-+eE0-9.]+), pz=([-+eE0-9.]+))")
    with open(path, 'r') as f:
        for linea in f:
            match = patron.match(linea.strip())
            if match:
                id_cuerpo = int(match.group(1))
                # Check which format was matched
                if match.group(2) is not None:  # Sequential format
                    px = np.float128(match.group(2))
                    py = np.float128(match.group(3))
                    pz = np.float128(match.group(4))
                else:  # Hybrid format
                    px = np.float128(match.group(5))
                    py = np.float128(match.group(6))
                    pz = np.float128(match.group(7))
                cuerpos[id_cuerpo] = (px, py, pz)
    return cuerpos


datos = {nombre: cargar_posiciones(
    archivos[nombre]) for nombre in orden_claves}

comparados = set()
pares = []
for i, a in enumerate(orden_claves):
    for b in orden_claves[i+1:]:
        if (a, b) not in comparados and (b, a) not in comparados:
            pares.append((a, b))
            comparados.add((a, b))

nombre_log = f"comparaciones_{N_OBJETIVO}.log"
with open(nombre_log, "w") as log:
    for a, b in pares:
        log.write(f"Comparando: {a} vs {b}\n")
        posiciones_a = datos[a]
        posiciones_b = datos[b]

        ids_comunes = set(posiciones_a.keys()) & set(posiciones_b.keys())

        if not ids_comunes:
            log.write("No hay cuerpos en común entre estos dos archivos.\n\n")
            continue

        max_dif_x = (None, 0.0)
        max_dif_y = (None, 0.0)
        max_dif_z = (None, 0.0)
        suma_dif_x = 0.0
        suma_dif_y = 0.0
        suma_dif_z = 0.0
        suma_error_pct_x = 0.0
        suma_error_pct_y = 0.0
        suma_error_pct_z = 0.0
        eps = 1e-12

        log.write(
            f"{'ID':<3} {'ΔX':>10} {'ΔY':>10} {'ΔZ':>12} {'eX':>10} {'eY':>10} {'eZ':>10}\n")
        for i in sorted(ids_comunes):
            x1, y1, z1 = posiciones_a[i]
            x2, y2, z2 = posiciones_b[i]

            dx = abs(x1 - x2)
            dy = abs(y1 - y2)
            dz = abs(z1 - z2)

            suma_dif_x += dx
            suma_dif_y += dy
            suma_dif_z += dz

            ex = (dx / (abs(x1) + eps)) * 100
            ey = (dy / (abs(y1) + eps)) * 100
            ez = (dz / (abs(z1) + eps)) * 100

            suma_error_pct_x += ex
            suma_error_pct_y += ey
            suma_error_pct_z += ez

            if dx > max_dif_x[1]:
                max_dif_x = (i, dx)
            if dy > max_dif_y[1]:
                max_dif_y = (i, dy)
            if dz > max_dif_z[1]:
                max_dif_z = (i, dz)

            if ex > 0.0000009 or ey > 0.0000009 or ez > 0.0000009:
                log.write(
                    f"{i:<5} {dx:12.6e} {dy:12.6e} {dz:12.6e} {ex:.6f}% {ey:.6f}% {ez:.6f}%\n")

        num_cuerpos = len(ids_comunes)
        prom_x = suma_dif_x / num_cuerpos
        prom_y = suma_dif_y / num_cuerpos
        prom_z = suma_dif_z / num_cuerpos

        prom_pct_x = suma_error_pct_x / num_cuerpos
        prom_pct_y = suma_error_pct_y / num_cuerpos
        prom_pct_z = suma_error_pct_z / num_cuerpos

        if a.startswith('secuencial') or b.startswith('secuencial'):
            base = a if a.startswith('secuencial') else b
            otro = b if a == base else a
            tiempo_seq = tiempos.get(base, 0.0)
            tiempo_otro = tiempos.get(otro, 0.0)
            speedup = tiempo_seq / tiempo_otro if tiempo_otro > 0 else 0

            procesos = 1
            hilos = 1
            if otro.startswith('pthreads'):
                match = re.search(r'pthreads-(\d+)-', otro)
                if match:
                    hilos = int(match.group(1))
            elif otro.startswith('hybrid'):
                match = re.search(r'hybrid-(\d+)-(\d+)-', otro)
                if match:
                    procesos = int(match.group(1))
                    hilos = int(match.group(2))

            eficiencia = speedup / (hilos * procesos) * 100

            log.write("\nResumen:\n\n")
            log.write(f"Cantidad de cuerpos procesados: {num_cuerpos}\n")
            log.write(f"Tiempo secuencial: {tiempo_seq:.6f} s\n")
            log.write(f"Tiempo {otro}: {tiempo_otro:.6f} s\n\n\n")
            log.write(f"Speedup: {speedup:.4f}\n")
            log.write(f"Eficiencia: {eficiencia:.4f}%\n\n\n")
        else:
            log.write("\nResumen:\n")
            log.write(f"Cantidad de cuerpos procesados: {num_cuerpos}\n\n\n")

        log.write(
            f"Máxima diferencia en X: cuerpo {max_dif_x[0]} con ΔX = {max_dif_x[1]:.6e}\n")
        log.write(
            f"Máxima diferencia en Y: cuerpo {max_dif_y[0]} con ΔY = {max_dif_y[1]:.6e}\n")
        log.write(
            f"Máxima diferencia en Z: cuerpo {max_dif_z[0]} con ΔZ = {max_dif_z[1]:.6e}\n\n")
        log.write(f"Promedio ΔX = {prom_x:.6e}\n")
        log.write(f"Promedio ΔY = {prom_y:.6e}\n")
        log.write(f"Promedio ΔZ = {prom_z:.6e}\n\n")
        log.write(f"Error porcentual promedio ΔX% = {prom_pct_x:.6f}%\n")
        log.write(f"Error porcentual promedio ΔY% = {prom_pct_y:.6f}%\n")
        log.write(f"Error porcentual promedio ΔZ% = {prom_pct_z:.6f}%\n\n")
        log.write("-" * 60 + "\n\n")
