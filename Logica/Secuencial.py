import gc
import psutil
import os
import argparse
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import time
import Utilidades as util
import matplotlib

# Cambiar el backend a Agg
matplotlib.use('Agg')

def read_fasta(file_path):
    """Lee una o más secuencias de un archivo FASTA y las concatena"""
    sequences = []  # Inicializa una lista vacía para almacenar las secuencias
    with open(file_path, "r") as file:  # Abre el archivo FASTA en modo de lectura
        for seq_record in SeqIO.parse(file, "fasta"):  # Itera sobre cada registro de secuencia en el archivo
            sequences.append(str(seq_record.seq))  # Convierte la secuencia a una cadena y la añade a la lista
    return ''.join(sequences)  # Concatena todas las secuencias y las devuelve como una sola cadena


def estimate_remaining_time(start_time, current_step, total_steps):
    """Estima el tiempo restante para completar la ejecución"""
    elapsed_time = time.time() - start_time  # Calcula el tiempo transcurrido desde el inicio
    if elapsed_time < 1e-6:  # Previene la división por cero
        return float('inf')
    steps_per_second = current_step / elapsed_time  # Calcula la tasa de pasos por segundo
    remaining_steps = total_steps - current_step  # Calcula los pasos restantes
    if steps_per_second <= 0:  # Si no hay progreso, evita la división por cero
        return float('inf')
    estimated_remaining_time = remaining_steps / steps_per_second  # Estima el tiempo restante
    return estimated_remaining_time  # Devuelve el tiempo restante estimado


def monitor_memory(used_bytes, max_system_percentage):
    """Monitoriza el uso de memoria y devuelve True si se excede el umbral"""
    process = psutil.Process(os.getpid())  # Obtiene el proceso actual
    memory_info = process.memory_info()  # Obtiene la información de memoria del proceso
    return memory_info.rss > used_bytes or psutil.virtual_memory().percent > max_system_percentage
    # Devuelve True si el uso de memoria del proceso excede used_bytes o si el porcentaje de memoria del sistema excede max_system_percentage


def create_dotplot(seq1, seq2, output_dir, window_size=500, step_size=100, memory_limit=5*1024**3, max_system_memory=85):
    """Crea y dibuja un dotplot por partes usando una ventana deslizante"""
    len1, len2 = len(seq1), len(seq2)  # Obtiene las longitudes de las secuencias
    total_steps = ((len1 - window_size) // step_size + 1) * ((len2 - window_size) // step_size + 1)
    # Calcula el número total de pasos para el dotplot
    current_step = 0  # Inicializa el contador de pasos
    start_time = time.time()  # Registra el tiempo de inicio
    last_save_time = start_time  # Registra el último tiempo de guardado
    save_interval = 60 * 60  # Intervalo de guardado de 60 minutos
    image_counter = 0  # Contador de imágenes guardadas

    os.makedirs(output_dir, exist_ok=True)  # Crea el directorio de salida si no existe
    elapsed_times = []  # Lista para almacenar los tiempos transcurridos de cada paso

    for i in range(0, len1 - window_size + 1, step_size):  # Itera sobre la primera secuencia con el tamaño de paso
        for j in range(0, len2 - window_size + 1, step_size):  # Itera sobre la segunda secuencia con el tamaño de paso

            if monitor_memory(memory_limit, max_system_memory):  # Monitoriza el uso de memoria
                print("Uso de memoria cerca del límite. Liberando memoria.")
                gc.collect()  # Libera la memoria recolectando basura
                continue  # Continúa al siguiente paso si el uso de memoria es alto

            step_start_time = time.time()  # Registra el tiempo de inicio del paso

            # Convertir las secuencias a arrays de bytes (más compacto que arrays de caracteres)
            window1 = np.frombuffer(seq1[i:i+window_size].encode('utf-8'), dtype=np.uint8)
            window2 = np.frombuffer(seq2[j:j+window_size].encode('utf-8'), dtype=np.uint8)

            dotplot = (window1[:, None] == window2).astype(np.int8)  # Genera el dotplot como una matriz de comparación

            fig, ax = plt.subplots(figsize=(10, 10))  # Crea una figura y un eje para el gráfico
            ax.imshow(dotplot, cmap='gray_r', extent=(j, j+window_size, i, i+window_size), origin='lower')  # Muestra el dotplot

            step_elapsed_time = time.time() - step_start_time  # Calcula el tiempo transcurrido para este paso
            elapsed_times.append(step_elapsed_time)  # Añade el tiempo transcurrido a la lista

            current_step += 1  # Incrementa el contador de pasos
            if current_step % 2000 == 0:  # Cada 2000 pasos, realiza acciones adicionales
                elapsed_time = time.time() - start_time  # Calcula el tiempo total transcurrido
                remaining_time = estimate_remaining_time(start_time, current_step, total_steps)
                # Estima el tiempo restante
                print(f"Progreso: {current_step}/{total_steps} pasos completados. Tiempo transcurrido: {elapsed_time:.2f} segundos. Tiempo restante estimado: {remaining_time:.2f} segundos")

                current_time = time.time()  # Obtiene el tiempo actual
                if current_time - last_save_time >= save_interval:  # Si ha pasado el intervalo de guardado
                    image_counter += 1  # Incrementa el contador de imágenes
                    image_path = util.create_unique_filename(output_dir, 'dotplot', 'png')  # Crea un nombre de archivo único
                    plt.savefig(image_path)  # Guarda la imagen del dotplot
                    print(f"Imagen intermedia guardada en {image_path}")
                    last_save_time = current_time  # Actualiza el último tiempo de guardado

                gc.collect()

            plt.close(fig)  # Cierra la figura para liberar memoria

    final_image_path = os.path.join(output_dir, 'dotplot_final.png')  # Define el path del dotplot final
    plt.savefig(final_image_path)  # Guarda el dotplot final
    print(f"Dotplot final guardado en {final_image_path}")
    gc.collect()  # Libera la memoria recolectando basura


def main(seq1_path, seq2_path, memory_limit, max_system_memory):
    overall_start_time = time.time()
    
    print("Leyendo secuencias...")
    read_start_time = time.time()
    seq1 = read_fasta(seq1_path)
    seq2 = read_fasta(seq2_path)
    read_elapsed_time = time.time() - read_start_time
    print(f"Tiempo de carga de datos: {read_elapsed_time:.2f} segundos")
    
    output_dir = os.path.join('.', 'Imagenes', 'Secuencial')
    
    print("Creando dotplot...")
    dotplot_start_time = time.time()
    create_dotplot(seq1, seq2, output_dir, memory_limit=memory_limit, max_system_memory=max_system_memory)
    dotplot_elapsed_time = time.time() - dotplot_start_time
    print(f"Tiempo de generación de imagen: {dotplot_elapsed_time:.2f} segundos")
    
    overall_elapsed_time = time.time() - overall_start_time
    print(f"Tiempo total de ejecución: {overall_elapsed_time:.2f} segundos")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Genera un dotplot a partir de dos secuencias en formato FASTA.')
    parser.add_argument('seq1', help='Ruta del archivo FASTA de la primera secuencia')
    parser.add_argument('seq2', help='Ruta del archivo FASTA de la segunda secuencia')
    parser.add_argument('--memory_limit', type=int, default=5, help='Límite de memoria (GB)')
    parser.add_argument('--max_system_memory', type=int, default=85, help='Porcentaje máximo de uso de memoria del sistema (RAM)')
    
    args = parser.parse_args()
    
    main(args.seq1, args.seq2, args.memory_limit * 1024**3, args.max_system_memory)

    
#EJECUTAR python Logica\Secuencial.py Archivos\E_coli.fna Archivos\Salmonella.fna
#EJECUTAR EN UBUNTU python3 Logica/Secuencial.py Archivos/E_coli.fna Archivos/Salmonella.fna