import csv
import os

def process_points_from_csv_to_txt(csv_filepath, output_txt_filepath="coordenadas.txt"):
    """
    Lee un archivo CSV con puntos (etiqueta, x, y) y guarda las coordenadas X e Y
    en dos cadenas de texto separadas en un archivo TXT.

    Args:
        csv_filepath (str): La ruta al archivo CSV de entrada.
        output_txt_filepath (str): La ruta donde se guardará el archivo TXT de salida.
                                   Por defecto es 'coordenadas.txt'.
    """
    x_coords = []
    y_coords = []

    try:
        # 1. Leer el archivo CSV
        with open(csv_filepath, mode='r', newline='', encoding='utf-8') as csvfile:
            reader = csv.DictReader(csvfile)
            
            # Verificar si las columnas 'x' y 'y' existen
            if 'x' not in reader.fieldnames or 'y' not in reader.fieldnames:
                print(f"Error: El archivo CSV '{csv_filepath}' debe contener las columnas 'x' y 'y'.")
                return

            for row in reader:
                try:
                    # Intentar convertir las coordenadas a float
                    x_coords.append(float(row['x']))
                    y_coords.append(float(row['y']))
                except ValueError:
                    print(f"Advertencia: No se pudieron convertir las coordenadas en la fila: {row}. Se omitirá esta fila.")
                    continue

        # 2. Formatear las coordenadas en cadenas de texto
        x_string = f"x={x_coords}"
        y_string = f"y={y_coords}"

        # 3. Guardar las cadenas en un archivo TXT
        with open(output_txt_filepath, mode='w', encoding='utf-8') as txtfile:
            txtfile.write(x_string + "\n") # Escribir la cadena X, seguida de un salto de línea
            txtfile.write(y_string + "\n") # Escribir la cadena Y
        
        print(f"Coordenadas X e Y guardadas exitosamente en '{output_txt_filepath}'.")

    except FileNotFoundError:
        print(f"Error: El archivo CSV '{csv_filepath}' no fue encontrado. Verifica la ruta.")
    except Exception as e:
        print(f"Ocurrió un error inesperado: {e}")

# --- Cómo usar el script ---
if __name__ == "__main__":
    # **IMPORTANTE:** Reemplaza 'puntos_geogebra.csv' con la ruta real de tu archivo CSV.
    # Asegúrate de que este archivo esté en la misma carpeta que el script o proporciona la ruta completa.
    input_csv_file = 'd:/Escritorio/alfredo/puntos_geogebra.csv' # <--- ¡Modifica esta línea con la ruta de tu CSV!
    output_text_file = 'd:/Escritorio/alfredo/coordenadas_xy.txt'
    # Ejecutar la función principal para procesar el CSV
    process_points_from_csv_to_txt(input_csv_file, output_text_file)