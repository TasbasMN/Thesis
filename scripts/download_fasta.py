import os
import requests
from tqdm import tqdm

base_url = "https://ftp.ensembl.org/pub/grch37/release-110/fasta/homo_sapiens/dna/"
chromosome_numbers = list(range(1, 23)) + ["MT", "X", "Y"]
download_directory = "fasta/grch37/"  # Replace with the desired directory path

# Create the download directory if it doesn't exist
os.makedirs(download_directory, exist_ok=True)

for chromosome_number in chromosome_numbers:
    file_url = f"{base_url}Homo_sapiens.GRCh37.dna.chromosome.{chromosome_number}.fa.gz"
    file_name = file_url.split("/")[-1]
    file_path = os.path.join(download_directory, file_name)

    response = requests.get(file_url, stream=True)

    total_size = int(response.headers.get("content-length", 0))
    block_size = 1024  # 1 KB

    progress_bar = tqdm(total=total_size, unit="B", unit_scale=True, desc=file_name, ncols=80)

    with open(file_path, "wb") as file:
        for data in response.iter_content(block_size):
            progress_bar.update(len(data))
            file.write(data)

    progress_bar.close()
