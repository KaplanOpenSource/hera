import zipfile

import os
import zipfile

from hera.utils.jsonutils import loadJSON


def add_directory_to_zip(zipf, folder_path, zip_path=""):
    """
    Recursively adds a folder to the zip file.
    zip_path is the path inside the zip.
    """
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            full_path = os.path.join(root, file)
            # Compute the relative path to keep folder structure
            relative_path = os.path.relpath(full_path, folder_path)
            zipf.write(full_path, os.path.join(zip_path, relative_path))


def zip_items(zip_filename, items):
    """
    Create a zip file and add files, directories, dictionaries.


    Parameters
    ----------
    zip_filename : str
        The string name. If the file does not end with .zip add it.

    items : list of [string,dict]
        If string: must be either a file on disk or a name of dict.
        if dict - adds all the keys as files whole content is the value.

    Returns
    -------

    """
    if not zip_filename.endswith(".zip"):
        zip_filename = f"{zip_filename}.zip"

    with zipfile.ZipFile(zip_filename, 'w', compression=zipfile.ZIP_DEFLATED) as zipf:
        for item in items:
            if isinstance(item,str):
                if os.path.isfile(item):
                    # Existing file
                    zipf.write(item, arcname=os.path.basename(item))
                elif os.path.isdir(item):
                    # Directory recursively
                    for root, dirs, files in os.walk(item):
                        for file in files:
                            full_path = os.path.join(root, file)
                            relative_path = os.path.relpath(full_path, start=os.path.dirname(item))
                            zipf.write(full_path, arcname=relative_path)
                else:
                    raise(f"Type Error: The item {item} is not a file or a directory")
            elif isinstance(item, dict):
                # Dict → each key is filename, value is content
                for filename, content in item.items():
                    if not isinstance(filename, str):
                        raise TypeError(f"Dict key must be a string, got {type(filename)}")
                    if not isinstance(content, str):
                        raise TypeError(f"Dict value must be a string, got {type(content)}")
                    zipf.writestr(filename, content)
            else:
                raise TypeError(f"Unsupported item type: {type(item)}")


def list_json_files_in_zip(zip_path):
    """Return a list of dicts with name and parsed content for each JSON file in a zip."""
    jsonFiles = []

    # Open the zip file in read mode
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        # List all file names inside the zip
        for file_name in zip_ref.namelist():
            if ".json" in file_name:

                with zip_ref.open(file_name) as file:
                    try:
                        content = file.read().decode('utf-8')
                    except UnicodeDecodeError:
                        print("Binary file, skipped reading content.")

                jsonFiles.append(dict(name=file_name,content=loadJSON((content))))

    return jsonFiles