import numpy as np
from scipy.ndimage import convolve
import matplotlib.pyplot as plt
from PIL import Image
import concurrent.futures

# Cargar la imagen
img = Image.open('Imagenes\Secuencial\dotplot_29.png').convert('L')
img_array = np.array(img)

# Definir kernels para detectar líneas diagonales
kernel1 = np.array([[2, -1, -1],
                    [-1, 2, -1],
                    [-1, -1, 2]])

kernel2 = np.array([[-1, -1, 2],
                    [-1, 2, -1],
                    [2, -1, -1]])

def apply_filter(img_section, kernel):
    return convolve(img_section, kernel)

# Dividir la imagen en secciones para procesamiento paralelo con superposición
def split_image_with_overlap(image, num_splits, overlap):
    sections = []
    step = image.shape[0] // num_splits
    for i in range(num_splits):
        start = i * step
        end = (i + 1) * step if i < num_splits - 1 else image.shape[0]
        start = max(0, start - overlap)
        end = min(image.shape[0], end + overlap)
        sections.append(image[start:end])
    return sections

num_splits = 4  # Número de divisiones para procesamiento paralelo
overlap = 10    # Superposición de 10 píxeles
image_sections = split_image_with_overlap(img_array, num_splits, overlap)

# Procesar las secciones de la imagen en paralelo
filtered_sections1 = []
filtered_sections2 = []
with concurrent.futures.ThreadPoolExecutor() as executor:
    futures1 = [executor.submit(apply_filter, section, kernel1) for section in image_sections]
    futures2 = [executor.submit(apply_filter, section, kernel2) for section in image_sections]
    filtered_sections1 = [f.result() for f in futures1]
    filtered_sections2 = [f.result() for f in futures2]

# Eliminar la superposición y combinar las secciones filtradas
filtered_img1 = np.vstack([section[overlap:-overlap] if i not in (0, num_splits-1) else section for i, section in enumerate(filtered_sections1)])
filtered_img2 = np.vstack([section[overlap:-overlap] if i not in (0, num_splits-1) else section for i, section in enumerate(filtered_sections2)])

# Mostrar las imágenes filtradas
plt.figure(figsize=(10, 10))
plt.subplot(1, 2, 1)
plt.imshow(filtered_img1, cmap='gray')
plt.title('Filtered Image - Diagonal 1')

plt.subplot(1, 2, 2)
plt.imshow(filtered_img2, cmap='gray')
plt.title('Filtered Image - Diagonal 2')

plt.show()
