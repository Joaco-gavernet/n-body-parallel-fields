import matplotlib.pyplot as plt
import numpy as np

# Datos
N_values = np.array([512, 1024, 2048, 4096])

# Tiempos de ejecución (en segundos)
tiempos_secuencial = [11.3079, 45.527, 182.3186, 726.365]
tiempos_4_threads = [3.3074, 12.7459, 49.8091, 200.0725]
tiempos_8_threads = [1.9905, 7.0129, 26.2362, 100.8778]

# Speedup
speedup_4 = [3.41897, 3.571894, 3.660347, 3.630509]
speedup_8 = [5.680934, 6.491894, 6.949124, 7.200444]

# Eficiencia
eficiencia_4 = [0.854742, 0.892973, 0.915087, 0.907627]
eficiencia_8 = [0.710117, 0.811487, 0.86864, 0.900056]

# Plot 1: Tiempos de ejecución
plt.figure(figsize=(6, 4))
plt.plot(N_values, tiempos_secuencial, 'o-', label='Secuencial')
plt.plot(N_values, tiempos_4_threads, 'o-', label='4 threads')
plt.plot(N_values, tiempos_8_threads, 'o-', label='8 threads')
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
