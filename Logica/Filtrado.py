import numpy as np
import cv2
from multiprocessing import Pool

def detect_diagonals(sub_image):
    """
    Esta función detecta las diagonales en una subimagen dada.

    Parámetros:
    - sub_image: La subimagen en la que se desea detectar las diagonales.

    Retorna:
    - La subimagen con las diagonales detectadas dibujadas en ella.
    """
    
    # Convertir la subimagen a escala de grises
    gray = cv2.cvtColor(sub_image, cv2.COLOR_BGR2GRAY)
    # Detectar bordes
    edges = cv2.Canny(gray, 50, 150, apertureSize=3)
    # Detectar líneas usando la Transformada de Hough
    lines = cv2.HoughLinesP(edges, 1, np.pi / 180, threshold=50, minLineLength=20, maxLineGap=5)
    # Dibujar las líneas detectadas
    if lines is not None:
        for line in lines:
            x1, y1, x2, y2 = line[0]
            cv2.line(sub_image, (x1, y1), (x2, y2), (0, 0, 0), 2)
    return sub_image

def process_image_parallel(image, num_processes):
    """
    Procesa una imagen dividiéndola en subimágenes y las procesa en paralelo utilizando multiprocessing.

    Args:
        image (numpy.ndarray): La imagen a procesar.
        num_processes (int): El número de procesos a utilizar para el procesamiento en paralelo.

    Returns:
        numpy.ndarray: La imagen resultante después de procesar las subimágenes en paralelo.
    """
    
    # Dividir la imagen en subimágenes
    height, width = image.shape[:2]
    sub_height = height // num_processes

    # Crear subimágenes considerando toda la altura
    sub_images = [image[i * sub_height: (i + 1) * sub_height if i != num_processes - 1 else height] for i in range(num_processes)]

    # Usar multiprocessing para procesar las subimágenes en paralelo
    with Pool(processes=num_processes) as pool:
        processed_sub_images = pool.map(detect_diagonals, sub_images)

    # Crear una imagen en blanco para los resultados
    result_image = np.zeros_like(image)

    # Combinar las subimágenes procesadas en la imagen de resultado
    for i, sub_image in enumerate(processed_sub_images):
        start_row = i * sub_height
        end_row = start_row + sub_image.shape[0]
        result_image[start_row:end_row, :] = sub_image

    return result_image

def main():
    # Cargar la imagen del dotplot
    #image_path = 'Imagenes\Secuencial\dotplot_1.png'
    image_path = 'Imagenes\\empo\\dotplot_parcial_1200_2.png'
    image = cv2.imread(image_path)
    
    if image is None:
        print(f"Error: la imagen en {image_path} no se pudo cargar.")
        return

    # Procesar la imagen
    num_processes = 4  # Ajusta este valor según tus necesidades
    processed_image = process_image_parallel(image, num_processes)

    # Redimensionar la imagen para que se ajuste a la pantalla
    screen_res = 1280, 720  # Ajusta esta resolución a la de tu pantalla
    scale_width = screen_res[0] / processed_image.shape[1]
    scale_height = screen_res[1] / processed_image.shape[0]
    scale = min(scale_width, scale_height)
    window_width = int(processed_image.shape[1] * scale)
    window_height = int(processed_image.shape[0] * scale)

    cv2.namedWindow('Processed Image', cv2.WINDOW_NORMAL)
    cv2.resizeWindow('Processed Image', window_width, window_height)

    # Visualizar la imagen procesada
    cv2.imshow('Processed Image', processed_image)
    cv2.waitKey(0)
    cv2.destroyAllWindows()

if __name__ == '__main__':
    main()
