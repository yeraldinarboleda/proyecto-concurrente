import matplotlib.pyplot as plt
import os
from Bio import SeqIO
import time
import numpy as np

def create_unique_filename(base_dir, base_name, extension):
    """Genera un nombre de archivo único en el directorio especificado"""
    counter = 1
    while True:
        filename = f"{base_name}_{counter}.{extension}"
        full_path = os.path.join(base_dir, filename)
        if not os.path.exists(full_path):
            return full_path
        counter += 1

def calculate_speedup(elapsed_times):
    """Calcula el speedup de los tiempos de ejecución."""
    seq_time = elapsed_times[0]  # Supone que el primer tiempo es el secuencial completo
    min_time = 1e-6  # Un valor mínimo para evitar divisiones por cero
    speedup = [seq_time / max(t, min_time) for t in elapsed_times]
    return speedup

def calculate_efficiency(speedup, num_processes):
    """Calcula la eficiencia utilizando el speed-up y la cantidad de procesos"""
    efficiency = [s / p for s, p in zip(speedup, num_processes)]
    return efficiency

def calculate_scalability(speedup, num_processes):
    """Calcula la escalabilidad utilizando el speed-up y la cantidad de procesos"""
    scalability = [s / p for s, p in zip(speedup, num_processes)]
    return scalability

def save_performance_graphs(elapsed_times, output_dir, num_processes):
    """Genera y guarda las gráficas de aceleración, eficiencia y escalabilidad"""
    
    speedup = calculate_speedup(elapsed_times)
    efficiency = calculate_efficiency(speedup, num_processes)
    scalability = calculate_scalability(speedup, num_processes)

    plt.figure()
    plt.plot(num_processes, speedup)
    plt.xlabel('Número de Procesos')
    plt.ylabel('Aceleración')
    plt.title('Aceleración')
    plt.savefig(create_unique_filename(output_dir, 'aceleracion', 'png'))
    plt.close()

    plt.figure()
    plt.plot(num_processes, efficiency)
    plt.xlabel('Número de Procesos')
    plt.ylabel('Eficiencia')
    plt.title('Eficiencia')
    plt.savefig(create_unique_filename(output_dir, 'eficiencia', 'png'))
    plt.close()

    plt.figure()
    plt.plot(num_processes, scalability)
    plt.xlabel('Número de Procesos')
    plt.ylabel('Escalabilidad')
    plt.title('Escalabilidad')
    plt.savefig(create_unique_filename(output_dir, 'escalabilidad', 'png'))
    plt.close()
    
def read_fasta(file_path):
    """Lee una o más secuencias de un archivo FASTA y las concatena"""
    sequences = []
    with open(file_path, "r") as file:
        for seq_record in SeqIO.parse(file, "fasta"):
            sequences.append(str(seq_record.seq))
    return ''.join(sequences)

def estimate_remaining_time(start_time, current_step, total_steps):
    """Estima el tiempo restante para completar la ejecución"""
    elapsed_time = time.time() - start_time
    if elapsed_time < 1e-6:
        return float('inf')
    steps_per_second = current_step / elapsed_time
    remaining_steps = total_steps - current_step
    if steps_per_second <= 0:
        return float('inf')
    estimated_remaining_time = remaining_steps / steps_per_second
    return estimated_remaining_time

