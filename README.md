﻿# proyecto-concurrente:

## Estructura del Proyecto

- **Archivos**: Contiene los archivos .fna de Salmonella y E. coli.
- **Imagenes**: Contiene los dotplots y las gráficas de desempeño.
- **Logica**: Contiene el código fuente de las diferentes implementaciones.

## PARA CLONAR EL REPOSITORIO
### git clone https://github.com/yeraldinarboleda/proyecto-concurrente.git

## EJECUTAR EL CODIGO SECUENCIAL EN LA TERMINAL:
### python Logica\Secuencial.py Archivos\E_coli.fna Archivos\Salmonella.fna

## EJECUTAR EL CODIGO MULTIPROCESSING EN LA TERMINAL:
### python Logica\CodigoMultiprocessing.py

## EJECUTAR EL CODIGO MPI EN LA TERMINAL:
### mpiexec -n 2 python Logica\CodigoMPI4PY.py Archivos/E_coli.fna Archivos/Salmonella.fna

## PARA EJECUTAR LA FUNCION DE FILTRADO SE CORRE ESTANDO EN EL ARCHIVO Filtrado.py
### Ademas se tiene que pasar la ruta de la imagen en el main para la variable image_path a la cual se le va a hacer el filtrado 

## LIBRERÍAS USADAS
- gc (Garbage Collection interface)
- psutil (process and system utilities)
- os (Miscellaneous operating system interfaces)
- argparse (Parser for command-line options, arguments, and subcommands)
- Bio (Biopython - Computational biology library)
- numpy (Numerical Python)
- matplotlib.pyplot (Matplotlib - Plotting library)
- time (Time access and conversions)
- mpi4py (MPI for Python)
- multiprocessing (Process-based parallelism)
- matplotlib (Matplotlib - Plotting library)
- cv2 (OpenCV - Open Source Computer Vision Library)

# INSTALACIÓN DE LIBRERIAS
- pip install psutil
- pip install biopython
- pip install numpy
- pip install matplotlib
- pip install mpi4py
- pip install opencv-python
