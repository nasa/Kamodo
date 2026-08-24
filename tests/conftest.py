import pytest
import requests
import tarfile
from pathlib import Path

BASE_URL = "https://ccmc.gsfc.nasa.gov/publicData/KAMODO_TestData"

@pytest.fixture(scope="session")
def get_test_data():
    """
    A factory fixture that returns a function to download and extract 
    specific model test data on demand.
    """
    def _download_and_extract(model_name):
        tests_dir = Path(__file__).parent
        data_dir = tests_dir / 'TestData'
        model_dir = data_dir / model_name
        tar_path = data_dir / f"{model_name}.tgz"

        # If the model's folder already exists, skip the download!
        if model_dir.exists():
            return model_dir

        # Ensure the TestData directory exists
        data_dir.mkdir(exist_ok=True)
        
        url = f"{BASE_URL}/{model_name}.tgz"
        print(f"\n[Fixture] Downloading {model_name} data from {url}...")
        
        response = requests.get(url, stream=True)
        # Raise an error if the URL is wrong or the file is missing (e.g. 404)
        response.raise_for_status() 

        # Write the tarball to disk
        with open(tar_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)

        print(f"\n[Fixture] Extracting {model_name}.tgz...")
        # Open and extract the gzipped tarball
        with tarfile.open(tar_path, "r:gz") as tar:
            # Note: in Python 3.12+, adding filter='data' is recommended for security, 
            # but standard extractall() works fine for trusted internal data.
            tar.extractall(path=data_dir)

        # Clean up the .tgz file to save disk space
        tar_path.unlink()
        
        return model_dir

    return _download_and_extract

