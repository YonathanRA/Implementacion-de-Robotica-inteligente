import xml.etree.ElementTree as ET
import csv
import os

def extract_points_from_xml_file_to_csv(xml_filepath, output_csv_filepath="points_from_xml.csv"):
    """
    Lee un archivo XML de GeoGebra (geogebra.xml), extrae los puntos y los guarda en un archivo CSV.

    Args:
        xml_filepath (str): La ruta al archivo XML (por ejemplo, 'geogebra.xml').
        output_csv_filepath (str): La ruta donde se guardará el archivo CSV.
                                   Por defecto es 'points_from_xml.csv'.
    """
    points_data = []

    try:
        # 1. Parsear el archivo XML
        tree = ET.parse(xml_filepath)
        root = tree.getroot()

        # 2. Buscar y extraer los datos de los puntos
        # La ruta XPath busca elementos de tipo 'point' dentro de 'construction'.
        # Esto asegura que se extraigan en el orden en que aparecen en el XML.
        for element in root.findall(".//construction/element[@type='point']"):
            label = element.get('label')
            coords_element = element.find('coords')
            
            if coords_element is not None:
                x = coords_element.get('x')
                y = coords_element.get('y')
                
                if label and x is not None and y is not None:
                    points_data.append({'label': label, 'x': float(x), 'y': float(y)})

        # 3. Generar el archivo CSV
        if points_data:
            with open(output_csv_filepath, 'w', newline='', encoding='utf-8') as csvfile:
                fieldnames = ['label', 'x', 'y']
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

                writer.writeheader() # Escribir la fila de encabezado
                for point in points_data:
                    writer.writerow(point) # Escribir cada fila de punto
            print(f"Puntos extraídos exitosamente de '{xml_filepath}' y guardados en '{output_csv_filepath}'.")
        else:
            print(f"No se encontraron puntos en el archivo XML '{xml_filepath}'.")

    except FileNotFoundError:
        print(f"Error: El archivo XML '{xml_filepath}' no fue encontrado. Verifica la ruta.")
    except ET.ParseError as e:
        print(f"Error al parsear el XML de '{xml_filepath}': {e}. El archivo podría estar corrupto o mal formado.")
    except Exception as e:
        print(f"Ocurrió un error inesperado: {e}")

# --- Cómo usar el script ---
if __name__ == "__main__":
    # **IMPORTANTE:** Reemplaza 'geogebra.xml' con la ruta real de tu archivo XML.
    # Asegúrate de que este archivo esté en la misma carpeta que el script o proporciona la ruta completa.
    xml_input_file = 'd:/Escritorio/alfredo/geogebra.xml' # <--- ¡Modifica esta línea con la ruta de tu archivo XML!
    output_csv_file = 'd:/Escritorio/alfredo/puntos_geogebra.csv'
    # Ejecutar la función principal para extraer los puntos
    extract_points_from_xml_file_to_csv(xml_input_file, output_csv_file)