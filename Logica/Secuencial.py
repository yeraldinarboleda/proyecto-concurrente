import argparse
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import time
import os
import gc
import psutil

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
    
    if elapsed_time < 1e-6:  # Avoid division by zero or very small elapsed time
        return float('inf')  # Return infinity since no meaningful estimation can be made
    
    steps_per_second = current_step / elapsed_time
    remaining_steps = total_steps - current_step
    
    if steps_per_second <= 0:
        return float('inf')  # Return infinity if steps_per_second is zero or negative
    
    estimated_remaining_time = remaining_steps / steps_per_second
    return estimated_remaining_time

def create_unique_filename(base_dir, base_name, extension):
    """Genera un nombre de archivo único en el directorio especificado"""
    counter = 1
    while True:
        filename = f"{base_name}_{counter}.{extension}"
        full_path = os.path.join(base_dir, filename)
        if not os.path.exists(full_path):
            return full_path
        counter += 1

def monitor_memory(used_percentage):
    """Monitoriza el uso de memoria y devuelve True si se excede el umbral"""
    memory_info = psutil.virtual_memory()
    return memory_info.percent > used_percentage

def create_dotplot(seq1, seq2, output_dir, window_size=300, step_size=100, memory_limit=80):
    """Crea y dibuja un dotplot por partes usando una ventana deslizante"""
    len1, len2 = len(seq1), len(seq2)
    
    total_steps = ((len1 - window_size) // step_size + 1) * ((len2 - window_size) // step_size + 1)
    current_step = 0
    
    start_time = time.time()
    last_save_time = start_time
    save_interval = 5 * 60  # 10 minutos en segundos
    image_counter = 0
    
    os.makedirs(output_dir, exist_ok=True)
    
    elapsed_times = []
    
    for i in range(0, len1 - window_size + 1, step_size):
        for j in range(0, len2 - window_size + 1, step_size):
            step_start_time = time.time()
            
            window1 = np.frombuffer(seq1[i:i+window_size].encode(), dtype=np.uint8)
            window2 = np.frombuffer(seq2[j:j+window_size].encode(), dtype=np.uint8)
            
            dotplot = (window1[:, None] == window2).astype(np.int8)
            
            fig, ax = plt.subplots()
            ax.imshow(dotplot, cmap='gray_r', extent=(j, j+window_size, i, i+window_size), origin='lower')
            plt.close(fig)  # Close the figure to release memory
            
            step_elapsed_time = time.time() - step_start_time
            elapsed_times.append(step_elapsed_time)
            
            current_step += 1
            if current_step % 500 == 0:  # Print progress every 500 steps
                elapsed_time = time.time() - start_time
                remaining_time = estimate_remaining_time(start_time, current_step, total_steps)
                print(f"Progreso: {current_step}/{total_steps} pasos completados. Tiempo transcurrido: {elapsed_time:.2f} segundos. Tiempo restante estimado: {remaining_time:.2f} segundos")
                
                current_time = time.time()
                if current_time - last_save_time >= save_interval:
                    image_counter += 1
                    image_path = create_unique_filename(output_dir, 'dotplot', 'png')
                    fig, ax = plt.subplots()
                    ax.imshow(dotplot, cmap='gray_r', extent=(j, j+window_size, i, i+window_size), origin='lower')
                    fig.savefig(image_path)
                    plt.close(fig)
                    print(f"Imagen intermedia guardada en {image_path}")
                    last_save_time = current_time
                    
                    # Guardar las gráficas de desempeño, aceleración, eficiencia y escalabilidad
                    save_performance_graphs(elapsed_times, image_counter, output_dir)
            
            # Liberar memoria manualmente y forzar la recolección de basura
            del dotplot, window1, window2
            gc.collect()
            
            # Monitorizar el uso de memoria
            if monitor_memory(memory_limit):
                print("Advertencia: Uso de memoria excedido. Reduciendo tamaño de ventana y paso.")
                window_size //= 2
                step_size //= 2
    
    print("Generación de dotplot completada.")
    
def save_performance_graphs(elapsed_times, image_counter, output_dir):
    """Genera y guarda las gráficas de desempeño, aceleración, eficiencia y escalabilidad"""
    # Gráfica de desempeño
    plt.figure()
    plt.plot(elapsed_times)
    plt.xlabel('Paso')
    plt.ylabel('Tiempo (s)')
    plt.title('Desempeño')
    plt.savefig(create_unique_filename(output_dir, 'desempeno', 'png'))
    plt.close()
    
    # Placeholder: Métricas de aceleración, eficiencia y escalabilidad
    # Aquí se pueden implementar cálculos de aceleración, eficiencia y escalabilidad basados en elapsed_times
    acceleration = [1] * len(elapsed_times)  # Placeholder, debe calcularse correctamente
    efficiency = [1] * len(elapsed_times)  # Placeholder, debe calcularse correctamente
    scalability = [1] * len(elapsed_times)  # Placeholder, debe calcularse correctamente
    
    # Gráfica de aceleración
    plt.figure()
    plt.plot(acceleration)
    plt.xlabel('Paso')
    plt.ylabel('Aceleración')
    plt.title('Aceleración')
    plt.savefig(create_unique_filename(output_dir, 'aceleracion', 'png'))
    plt.close()
    
    # Gráfica de eficiencia
    plt.figure()
    plt.plot(efficiency)
    plt.xlabel('Paso')
    plt.ylabel('Eficiencia')
    plt.title('Eficiencia')
    plt.savefig(create_unique_filename(output_dir, 'eficiencia', 'png'))
    plt.close()
    
    # Gráfica de escalabilidad
    plt.figure()
    plt.plot(scalability)
    plt.xlabel('Paso')
    plt.ylabel('Escalabilidad')
    plt.title('Escalabilidad')
    plt.savefig(create_unique_filename(output_dir, 'escalabilidad', 'png'))
    plt.close()

def main(seq1_path, seq2_path, memory_limit):
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
    create_dotplot(seq1, seq2, output_dir, memory_limit=memory_limit)
    dotplot_elapsed_time = time.time() - dotplot_start_time
    print(f"Tiempo de generación de imagen: {dotplot_elapsed_time:.2f} segundos")
    
    overall_elapsed_time = time.time() - overall_start_time
    print(f"Tiempo total de ejecución: {overall_elapsed_time:.2f} segundos")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Genera un dotplot a partir de dos secuencias en formato FASTA.')
    parser.add_argument('seq1', help='Ruta del archivo FASTA de la primera secuencia')
    parser.add_argument('seq2', help='Ruta del archivo FASTA de la segunda secuencia')
    parser.add_argument('--memory_limit', type=int, default=80, help='Porcentaje máximo de uso de memoria (RAM)')
    
    args = parser.parse_args()
    
    main(args.seq1, args.seq2, args.memory_limit)


    
#EJECUTAR python Logica\Secuencial.py Archivos\E_coli.fna Archivos\Salmonella.fna
#EJECUTAR EN UBUNTU python3 Logica/Secuencial.py Archivos/E_coli.fna Archivos/Salmonella.fna
