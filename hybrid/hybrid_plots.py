import matplotlib.pyplot as plt
import numpy as np

# Datos
N_values = np.array([512, 1024, 2048, 4096])

# Tiempos de ejecución (en segundos)
tiempos_secuencial = [11.3079, 45.527, 182.3186, 726.365]
tiempos_4_threads = [3.4906, 13.4793, 52.3246, 207.6854]
tiempos_8_threads = [2.2709, 7.4536, 27.5177, 120.9985]
tiempos_16_threads = [1.9655, 5.0805, 18.0252, 65.5426]

# Speedup
speedup_4 = [3.2395, 3.3775, 3.4844, 3.4974]
speedup_8 = [4.9795, 6.1081, 6.6255, 6.0031]
speedup_16 = [5.7532, 8.9611, 10.1147, 11.0823]

# Eficiencia
eficiencia_4 = [0.8099, 0.8444, 0.8711, 0.8744]
eficiencia_8 = [0.6224, 0.7635, 0.8282, 0.7504]
eficiencia_16 = [0.3596, 0.5601, 0.6322, 0.6926]

# Plot 1: Tiempos de ejecución
plt.figure(figsize=(6, 4))
plt.plot(N_values, tiempos_secuencial, 'o-', label='Secuencial')
plt.plot(N_values, tiempos_4_threads, 'o-', label='4 threads')
plt.plot(N_values, tiempos_8_threads, 'o-', label='8 threads')
plt.plot(N_values, tiempos_16_threads, 'o-', label='16 threads')
plt.xscale('log', base=2)
plt.xticks(N_values, N_values)
plt.ylabel('Tiempo (s)')
plt.xlabel('Tamaño del problema (N)')
plt.title('Tiempos de ejecución')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('plots/tiempos_ejecucion.png', dpi=300, bbox_inches='tight')
plt.close()

# Plot 2: Speedup
plt.figure(figsize=(6, 4))
plt.plot(N_values, speedup_4, 'o-', label='4 threads')
plt.plot(N_values, speedup_8, 'o-', label='8 threads')
plt.plot(N_values, speedup_16, 'o-', label='16 threads')
plt.xscale('log', base=2)
plt.xticks(N_values, N_values)
plt.ylabel('Speedup')
plt.xlabel('Tamaño del problema (N)')
plt.title('Speedup')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('plots/speedup.png', dpi=300, bbox_inches='tight')
plt.close()

# Plot 3: Eficiencia
plt.figure(figsize=(6, 4))
plt.plot(N_values, eficiencia_4, 'o-', label='4 threads')
plt.plot(N_values, eficiencia_8, 'o-', label='8 threads')
plt.plot(N_values, eficiencia_16, 'o-', label='16 threads')
plt.xscale('log', base=2)
plt.xticks(N_values, N_values)
plt.xlabel('Tamaño del problema (N)')
plt.ylabel('Eficiencia')
plt.title('Eficiencia')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('plots/eficiencia.png', dpi=300, bbox_inches='tight')
plt.close()

print("Plots have been saved as:")
print("- tiempos_ejecucion.png")
print("- speedup.png")
print("- eficiencia.png")
