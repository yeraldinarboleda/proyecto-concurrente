import numpy as np
import cv2
from multiprocessing import Pool

    
def aplicar_filtro(input_matrix, output_path):
    """
    Aplica un filtro personalizado a una matriz de entrada y guarda la imagen filtrada en la ruta de salida.

    Parámetros:
    - input_matrix: Matriz de entrada que representa una imagen.
    - output_path: Ruta de salida donde se guardará la imagen filtrada.

    El código realiza las siguientes operaciones:
    1. Define un kernel personalizado.
    2. Aplica el filtro 2D utilizando el kernel personalizado a la matriz de entrada.
    3. Normaliza la matriz filtrada en el rango de 0 a 127.
    4. Aplica un umbral a la matriz normalizada para obtener una imagen binaria.
    5. Guarda la imagen binaria en la ruta de salida.
    6. Muestra la imagen filtrada en una ventana.
    7. Espera a que se presione una tecla para cerrar la ventana.

    No devuelve ningún valor.
    """
    
    custom_kernel = np.array([[1, 0, 0],
                              [0, 1, 0],
                              [0, 0, 1]])

    filtered_matrix = cv2.filter2D(input_matrix, -1, custom_kernel)

    normalized_matrix = cv2.normalize(filtered_matrix, None, 0, 127, cv2.NORM_MINMAX)

    threshold_value = 126
    _, thresholded_matrix = cv2.threshold(normalized_matrix, threshold_value, 255, cv2.THRESH_BINARY)

    cv2.imwrite(output_path, thresholded_matrix)
    cv2.imshow('Filtered Image', thresholded_matrix)
    cv2.waitKey(0)
    cv2.destroyAllWindows()
    
# Cargar la imagen en una matriz
image_path = 'Imagenes\Secuencial\dotplot_17.png'  

#image_path = 'Imagenes\Multiprocessing\dotplot_parcial_1200_2.png' 
matrix = cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)  # Cargar la imagen en escala de grises

# Aplicar el filtro a la matriz
filtered_image_path = 'Imagenes\\Secuencial\\dotplot_17_filtered.png'
aplicar_filtro(matrix, filtered_image_path)

