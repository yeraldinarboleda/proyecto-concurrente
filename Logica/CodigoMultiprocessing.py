import multiprocessing as mp
import numpy as np
import matplotlib.pyplot as plt
import time
import os
import Utilidades as Util

    
def draw_dotplot(matrix, fig_name='dotplot.png'):
    """
    Dibuja un gráfico de puntos basado en la matriz dada.

    Parámetros:
    matrix (numpy.ndarray): La matriz que representa el gráfico de puntos.
    fig_name (str): El nombre del archivo de salida (por defecto es 'dotplot.png').
    """
    plt.figure(figsize=(10, 10))
    plt.imshow(matrix, cmap='Greys', aspect='auto')
    plt.ylabel("Salmonela")
    plt.xlabel("E coli")
    plt.savefig(fig_name)
    plt.close()

def worker(args):
    """
    Esta función toma una tupla de argumentos que incluye un índice 'i', una secuencia 'secuencia1' y otra secuencia 'secuencia2'.
    Luego, crea una fila de ceros del mismo tamaño que 'secuencia2'.
    A continuación, recorre la secuencia 'secuencia2' y si el elemento en la posición 'i' de 'secuencia1' es igual al elemento en la posición 'j' de 'secuencia2',
    asigna el valor 1 a la posición 'j' de la fila.
    Finalmente, devuelve la fila resultante.
    """
    i, secuencia1, secuencia2 = args
    row = np.zeros(len(secuencia2), dtype=np.int8)
    for j in range(len(secuencia2)):
        if secuencia1[i] == secuencia2[j]:
            row[j] = 1
    return row

def parallel_dotplot(secuencia1, secuencia2, n_cpu=mp.cpu_count(), batch_size=1000):
    """
    Calcula el dotplot de manera paralela utilizando multiprocessing.

    Parámetros:
    secuencia1 (str): La primera secuencia.
    secuencia2 (str): La segunda secuencia.
    n_cpu (int): El número de procesadores a utilizar (por defecto es el número de núcleos de la CPU).
    batch_size (int): El tamaño del lote de procesamiento (por defecto es 1000).

    Retorna:
    numpy.ndarray: La matriz del dotplot resultante.
    """
    total_steps = len(secuencia1)
    start_time = time.time()
    elapsed_times = []
    results = []
    
    with mp.Pool(processes=n_cpu) as pool:
        for i in range(0, total_steps, batch_size):
            batch_secuencia1 = secuencia1[i:i+batch_size]
            batch_results = pool.map(worker, [(j, batch_secuencia1, secuencia2) for j in range(len(batch_secuencia1))])
            results.extend(batch_results)
            
            elapsed_time = time.time() - start_time
            elapsed_times.append(elapsed_time)
            percentage = min((i + batch_size) / total_steps * 100, 100)
            print(f"Progreso: {percentage:.2f}%")
            print(f"Tiempo transcurrido: {elapsed_time:.2f} segundos")
            estimated_remaining_time = Util.estimate_remaining_time(start_time, i + batch_size, total_steps)
            print(f"Tiempo estimado restante: {estimated_remaining_time:.2f} segundos")
            
            # Guardar punto de control de la imagen del dotplot parcial
            partial_matrix = np.vstack(results).astype(np.int8)
            draw_dotplot(partial_matrix, f'Imagenes/Multiprocessing/dotplot_parcial_{i + batch_size}_{n_cpu}.png')

            # Guardar gráficas de rendimiento hasta el momento
            Util.save_performance_graphs(elapsed_times, 'Imagenes/Multiprocessing', list(range(1, len(elapsed_times) + 1)))
            # Limpiar resultados para liberar memoria
            results = []

    return np.vstack(results).astype(np.int8)

def main():
    base_dir = 'Archivos'
    output_dir = 'Imagenes/Multiprocessing'
    os.makedirs(output_dir, exist_ok=True)
    
    secuencia1_path = os.path.join(base_dir, 'Salmonella.fna')
    secuencia2_path = os.path.join(base_dir, 'E_coli.fna')
    
    print("Leyendo secuencias de archivos FASTA...")
    secuencia1 = Util.read_fasta(secuencia1_path)
    secuencia2 = Util.read_fasta(secuencia2_path)
    len1, len2 = len(secuencia1), len(secuencia2)
    
    print(f"Secuencia 1: {len1} caracteres")
    print(f"Secuencia 2: {len2} caracteres")
    
    n_proc = [1, mp.cpu_count()] # Ajustar según la cantidad de núcleos de la CPU
    strong_times = []
    weak_times = []

    # Medir rendimiento en ejecución fuerte
    for i in n_proc:
        print(f"Iniciando dotplot con {i} procesadores...")
        start_time = time.time()
        dotplot = parallel_dotplot(secuencia1, secuencia2, i)
        end_time = time.time()
        elapsed_time = end_time - start_time
        strong_times.append(elapsed_time)
        print(f"Tiempo con {i} procesadores: {elapsed_time:.2f} segundos")
    
    # Medir rendimiento en ejecución débil
    for i in n_proc:
        secuencia_long = secuencia2 * i
        print(f"Iniciando dotplot con {i} procesadores y secuencia escalada...")
        start_time = time.time()
        dotplot = parallel_dotplot(secuencia_long, secuencia_long, i)
        end_time = time.time()
        elapsed_time = end_time - start_time
        weak_times.append(elapsed_time)
        print(f"Tiempo con {i} procesadores (secuencia escalada): {elapsed_time:.2f} segundos")

    # Guardar dotplot final
    draw_dotplot(dotplot, os.path.join(output_dir, 'dotplot_final.png'))

    # Guardar gráficas de rendimiento
    Util.save_performance_graphs(strong_times, output_dir, n_proc)
    Util.save_performance_graphs(weak_times, output_dir, n_proc)

if __name__ == '__main__':
    main()

# python Logica\CodigoMultiprocessing.py