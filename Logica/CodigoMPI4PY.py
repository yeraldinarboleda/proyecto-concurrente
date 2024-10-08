import gc
import psutil
import os
import argparse
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import time
from mpi4py import MPI
import Utilidades as util
import matplotlib
##pruebas ...
print("Inicio del script")  # Mensaje de depuración

# Cambiar el backend a Agg
matplotlib.use('Agg')

def read_fasta(file_path):
    print(f"Leyendo archivo FASTA: {file_path}")  # Mensaje de depuración
    sequences = []
    with open(file_path, "r") as file:
        for seq_record in SeqIO.parse(file, "fasta"):
            sequences.append(str(seq_record.seq))
    return ''.join(sequences)

def estimate_remaining_time(start_time, current_step, total_steps):
    elapsed_time = time.time() - start_time
    if elapsed_time < 1e-6:
        return float('inf')
    steps_per_second = current_step / elapsed_time
    remaining_steps = total_steps - current_step
    if steps_per_second <= 0:
        return float('inf')
    estimated_remaining_time = remaining_steps / steps_per_second
    return estimated_remaining_time

def monitor_memory(used_bytes, max_system_percentage):
    process = psutil.Process(os.getpid())
    memory_info = process.memory_info()
    return memory_info.rss > used_bytes or psutil.virtual_memory().percent > max_system_percentage

def create_dotplot(seq1, seq2, output_dir, window_size=500, step_size=100, memory_limit=5*1024**3, max_system_memory=85):
    """
    Crea un dotplot a partir de dos secuencias de ADN.

    Parameters:
    - seq1 (str): La primera secuencia de ADN.
    - seq2 (str): La segunda secuencia de ADN.
    - output_dir (str): Directorio de salida donde se guardarán las imágenes del dotplot.
    - window_size (int): Tamaño de la ventana deslizante para el dotplot (por defecto: 500).
    - step_size (int): Tamaño del paso para deslizar la ventana (por defecto: 100).
    - memory_limit (int): Límite de memoria en bytes para controlar el uso de memoria (por defecto: 5 GB).
    - max_system_memory (int): Porcentaje máximo de memoria del sistema permitido (por defecto: 85).

    Returns:
    None

    El dotplot es una representación gráfica de las similitudes entre dos secuencias de ADN. 
    Este algoritmo utiliza la biblioteca MPI para realizar el cálculo de manera paralela en múltiples procesos.
    El dotplot se guarda como una imagen en el directorio de salida especificado.
    """
    len1, len2 = len(seq1), len(seq2)
    total_steps = ((len1 - window_size) // step_size + 1) * ((len2 - window_size) // step_size + 1)
    current_step = 0
    start_time = time.time()
    last_save_time = start_time
    save_interval = 60 * 60  # 60 minutos en segundos (1 hora)
    image_counter = 0

    os.makedirs(output_dir, exist_ok=True)
    elapsed_times = []

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    print(f"Proceso {rank} iniciando create_dotplot")  # Mensaje de depuración

    for i in range(rank, len1 - window_size + 1, step_size * size):
        for j in range(0, len2 - window_size + 1, step_size):

            if monitor_memory(memory_limit, max_system_memory):
                print(f"Proceso {rank}: Uso de memoria cerca del límite. Liberando memoria.")
                gc.collect()
                continue

            step_start_time = time.time()

            window1 = np.frombuffer(seq1[i:i+window_size].encode('utf-8'), dtype=np.uint8)
            window2 = np.frombuffer(seq2[j:j+window_size].encode('utf-8'), dtype=np.uint8)

            dotplot = (window1[:, None] == window2).astype(np.int8)

            fig, ax = plt.subplots(figsize=(10, 10))
            ax.imshow(dotplot, cmap='gray_r', extent=(j, j+window_size, i, i+window_size), origin='lower')

            step_elapsed_time = time.time() - step_start_time
            elapsed_times.append(step_elapsed_time)

            current_step += 1
            if current_step % 2000 == 0:
                elapsed_time = time.time() - start_time
                remaining_time = estimate_remaining_time(start_time, current_step, total_steps)
                print(f"Proceso {rank}: Progreso: {current_step}/{total_steps} pasos completados. Tiempo transcurrido: {elapsed_time:.2f} segundos. Tiempo restante estimado: {remaining_time:.2f} segundos")

                current_time = time.time()
                if current_time - last_save_time >= save_interval:
                    image_counter += 1
                    image_path = util.create_unique_filename(output_dir, f'dotplot_{rank}', 'png')
                    plt.savefig(image_path)
                    print(f"Proceso {rank}: Imagen intermedia guardada en {image_path}")
                    last_save_time = current_time
                    util.desempeño(elapsed_times, output_dir)

                gc.collect()

            plt.close(fig)  # Cierra la figura para liberar memoria

    if rank == 0:
        final_image_path = os.path.join(output_dir, 'dotplot_final.png')
        plt.savefig(final_image_path)
        print(f"Dotplot final guardado en {final_image_path}")
        gc.collect()

def main(seq1_path, seq2_path, memory_limit, max_system_memory):
    overall_start_time = time.time()
    
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    
    if rank == 0:
        print("Leyendo secuencias...")
    read_start_time = time.time()
    seq1 = read_fasta(seq1_path)
    seq2 = read_fasta(seq2_path)
    read_elapsed_time = time.time() - read_start_time
    if rank == 0:
        print(f"Tiempo de carga de datos: {read_elapsed_time:.2f} segundos")
    
    output_dir = os.path.join('.', 'Imagenes', 'MPI4PY')
    
    if rank == 0:
        print("Creando dotplot...")
    dotplot_start_time = time.time()
    create_dotplot(seq1, seq2, output_dir, memory_limit=memory_limit, max_system_memory=max_system_memory)
    dotplot_elapsed_time = time.time() - dotplot_start_time
    if rank == 0:
        print(f"Tiempo de generación de imagen: {dotplot_elapsed_time:.2f} segundos")
    
    overall_elapsed_time = time.time() - overall_start_time
    if rank == 0:
        print(f"Tiempo total de ejecución: {overall_elapsed_time:.2f} segundos")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Genera un dotplot a partir de dos secuencias en formato FASTA.')
    parser.add_argument('seq1', help='Ruta del archivo FASTA de la primera secuencia')
    parser.add_argument('seq2', help='Ruta del archivo FASTA de la segunda secuencia')
    parser.add_argument('--memory_limit', type=int, default=5, help='Límite de memoria (GB)')
    parser.add_argument('--max_system_memory', type=int, default=85, help='Porcentaje máximo de uso de memoria del sistema (RAM)')
    
    args = parser.parse_args()
    
    main(args.seq1, args.seq2, args.memory_limit * 1024**3, args.max_system_memory)



#EJECUTAR#mpiexec -n 2 python Logica\CodigoMPI4PY.py Archivos/E_coli.fna Archivos/Salmonella.fna
#EJECUTAR#mpiexec -n 2 python Logica\CodigoMPI4PY.py Archivos/prueba1.fna Archivos/prueba2.fna

