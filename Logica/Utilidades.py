import matplotlib.pyplot as plt
import os

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
    speedup = [seq_time / max(time, min_time) for time in elapsed_times]
    return speedup

def calculate_efficiency(speedup, num_processes):
    """Calcula la eficiencia utilizando el speed-up y la cantidad de procesos"""
    efficiency = [s / num_processes for s in speedup]
    return efficiency

def calculate_scalability(speedup, num_processes):
    """Calcula la escalabilidad utilizando el speed-up y la cantidad de procesos"""
    scalability = [s / num_processes for s in speedup]
    return scalability



def save_performance_graphs(elapsed_times, output_dir):
    """Genera y guarda las gráficas de desempeño, aceleración, eficiencia y escalabilidad"""
    plt.figure()
    plt.plot(elapsed_times)
    plt.xlabel('Paso')
    plt.ylabel('Tiempo (s)')
    plt.title('Desempeño')
    plt.savefig(create_unique_filename(output_dir, 'desempeno', 'png'))
    plt.close()

    speedup = calculate_speedup(elapsed_times)
    # efficiency = calculate_efficiency(speedup, num_processes)
    # scalability = calculate_scalability(speedup, num_processes)

    plt.figure()
    plt.plot(speedup)
    plt.xlabel('Paso')
    plt.ylabel('Aceleración')
    plt.title('Aceleración')
    plt.savefig(create_unique_filename(output_dir, 'aceleracion', 'png'))
    plt.close()

    # plt.figure()
    # plt.plot(efficiency)
    # plt.xlabel('Paso')
    # plt.ylabel('Eficiencia')
    # plt.title('Eficiencia')
    # plt.savefig(create_unique_filename(output_dir, 'eficiencia', 'png'))
    # plt.close()

    # plt.figure()
    # plt.plot(scalability)
    # plt.xlabel('Paso')
    # plt.ylabel('Escalabilidad')
    # plt.title('Escalabilidad')
    # plt.savefig(create_unique_filename(output_dir, 'escalabilidad', 'png'))
    # plt.close()