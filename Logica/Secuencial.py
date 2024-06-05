import argparse
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import time
import os

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

def create_dotplot(seq1, seq2, window_size=500, step_size=50):
    """Crea y dibuja un dotplot por partes usando una ventana deslizante"""
    len1, len2 = len(seq1), len(seq2)
    fig, ax = plt.subplots()
    
    total_steps = ((len1 - window_size) // step_size + 1) * ((len2 - window_size) // step_size + 1)
    current_step = 0
    
    start_time = time.time()
    last_save_time = start_time
    save_interval = 5 * 60  # 5 minutos en segundos
    image_counter = 0
    
    output_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'imagenes')
    os.makedirs(output_dir, exist_ok=True)
    
    for i in range(0, len1 - window_size + 1, step_size):
        for j in range(0, len2 - window_size + 1, step_size):
            window1 = np.array(list(seq1[i:i+window_size]))
            window2 = np.array(list(seq2[j:j+window_size]))
            
            dotplot = (window1[:, None] == window2).astype(np.int8)
            
            ax.imshow(dotplot, cmap='gray_r', extent=(j, j+window_size, i, i+window_size), origin='lower')
            
            current_step += 1
            if current_step % 100 == 0:  # Print progress every 100 steps
                elapsed_time = time.time() - start_time
                remaining_time = estimate_remaining_time(start_time, current_step, total_steps)
                print(f"Progreso: {current_step}/{total_steps} pasos completados. Tiempo transcurrido: {elapsed_time:.2f} segundos. Tiempo restante estimado: {remaining_time:.2f} segundos")
                
                current_time = time.time()
                if current_time - last_save_time >= save_interval:
                    image_counter += 1
                    image_path = os.path.join(output_dir, f'dotplot_{image_counter}.png')
                    plt.savefig(image_path)
                    print(f"Imagen intermedia guardada en {image_path}")
                    last_save_time = current_time
    
    ax.set_xlabel('Secuencia 2')
    ax.set_ylabel('Secuencia 1')
    ax.set_title('Dotplot')
    plt.show()

def main(seq1_path, seq2_path):
    start_time = time.time()
    print("Leyendo secuencias...")
    seq1 = read_fasta(seq1_path)
    seq2 = read_fasta(seq2_path)
    
    print("Creando dotplot...")
    create_dotplot(seq1, seq2)
    
    elapsed_time = time.time() - start_time
    print(f"Tiempo total de ejecución: {elapsed_time:.2f} segundos")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Genera un dotplot a partir de dos secuencias en formato FASTA.')
    parser.add_argument('seq1', help='Archivo de la primera secuencia en formato FASTA')
    parser.add_argument('seq2', help='Archivo de la segunda secuencia en formato FASTA')
    
    args = parser.parse_args()
    
    main(args.seq1, args.seq2)


    
#EJECUTAR python Logica\Secuencial.py Archivos\E_coli.fna Archivos\Salmonella.fna
#EJECUTAR EN UBUNTU python3 Logica/Secuencial.py Archivos/E_coli.fna Archivos/Salmonella.fna
