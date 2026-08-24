import argparse
import requests
import tarfile
import sys
import re
from pathlib import Path

BASE_URL = "https://ccmc.gsfc.nasa.gov/publicData/KAMODO_TestData"

def format_size(size_bytes):
    """Converts bytes to a human-readable MB or GB format."""
    if size_bytes == 0:
        return "Unknown Size"
    if size_bytes >= 1024**3:
        return f"{size_bytes / (1024**3):.1f} GB"
    else:
        return f"{size_bytes / (1024**2):.1f} MB"

def get_available_models(include_sizes=False):
    """Scrapes the CCMC directory page to find available .tgz files and optionally their sizes."""
    try:
        response = requests.get(BASE_URL + "/", timeout=10)
        response.raise_for_status()
        
        # Regex to find all href links that end in .tgz
        matches = re.findall(r'href="([^"]+)\.tgz"', response.text)
        models = sorted(list(set(matches)))
        
        if not include_sizes:
            return models
            
        # Fetch file sizes via HTTP HEAD requests
        model_info = []
        for m in models:
            try:
                head_resp = requests.head(f"{BASE_URL}/{m}.tgz", timeout=5)
                size_bytes = int(head_resp.headers.get('content-length', 0))
                model_info.append((m, format_size(size_bytes)))
            except requests.exceptions.RequestException:
                model_info.append((m, "Unknown Size"))
                
        return model_info
    except requests.exceptions.RequestException as e:
        print(f"[x] Could not fetch available models from {BASE_URL}. Error: {e}")
        return None

def print_custom_help():
    """Prints the script description, dynamically scraped models, and sizes."""
    print("\n" + "="*55)
    print(" Kamodo Test Data Downloader")
    print("="*55)
    print("This script downloads and extracts test data for Kamodo space weather models.")
    print("The data is pulled directly from the NASA CCMC public data server and saved")
    print("into the 'tests/TestData' directory, making it immediately ready for pytest.\n")
    
    print("Fetching available models and sizes from the server (this takes a second)...")
    model_info = get_available_models(include_sizes=True)
    
    if model_info:
        print("\nAvailable models on the server:")
        for m, size in model_info:
            print(f"  \u2022 {m:<15} ({size})")
    else:
        print("\nCould not retrieve the list of models dynamically. Please manually check:")
        print(BASE_URL)
        
    print("\nUsage:")
    print("  python download_TestData.py                 (Downloads default: ADELPHI)")
    print("  python download_TestData.py <model> <model> (Downloads specifically named models)")
    print("  python download_TestData.py ALL             (Downloads EVERY available model)")
    print("  python download_TestData.py help            (Shows this help message)\n")
    print("WARNING: Using 'ALL' or specifying heavy models may download 100GB+ of data.\n")

def download_and_extract(model_name, data_dir):
    """Downloads and extracts a specific model's test data."""
    model_dir = data_dir / model_name
    tar_path = data_dir / f"{model_name}.tgz"

    # Skip if already downloaded
    if model_dir.exists():
        print(f"[\u2713] {model_name} data already exists in {model_dir}. Skipping.")
        return

    # Ensure the TestData directory exists
    data_dir.mkdir(parents=True, exist_ok=True)
    
    url = f"{BASE_URL}/{model_name}.tgz"
    print(f"[\u2193] Downloading {model_name} from {url}...")
    
    try:
        response = requests.get(url, stream=True)
        response.raise_for_status() 
    except requests.exceptions.RequestException as e:
        print(f"[x] Failed to download {model_name}: {e}")
        return

    # Write the tarball to disk
    with open(tar_path, 'wb') as f:
        for chunk in response.iter_content(chunk_size=8192):
            f.write(chunk)

    print(f"[\u2935] Extracting {model_name}.tgz...")
    with tarfile.open(tar_path, "r:gz") as tar:
        if hasattr(tarfile, 'data_filter'):
            tar.extractall(path=data_dir, filter='data')
        else:
            tar.extractall(path=data_dir)

    # Clean up the .tgz file to save space
    tar_path.unlink()
    print(f"[\u2713] Successfully processed {model_name}.\n")

if __name__ == "__main__":
    # Intercept help requests to show the dynamic scrape screen
    if any(arg.lower() in ['help', '-h', '--help'] for arg in sys.argv):
        print_custom_help()
        sys.exit(0)

    parser = argparse.ArgumentParser(add_help=False) 
    parser.add_argument(
        "models", 
        nargs="*", 
        default=["ADELPHI"], 
    )
    
    args = parser.parse_args()
    
    # Process the ALL flag
    if any(m.upper() == "ALL" for m in args.models):
        print("Fetching list of all available models...")
        all_models = get_available_models(include_sizes=False)
        if all_models:
            args.models = all_models
            print(f"Queueing {len(args.models)} models for download...\n")
        else:
            print("[x] Failed to fetch the model list. Exiting.")
            sys.exit(1)
    
    # Anchor TestData directory relative to where this script is located
    root_dir = Path(__file__).resolve().parent
    target_data_dir = root_dir / "tests" / "TestData"
    
    print(f"Target data directory: {target_data_dir}\n")
    
    for model in args.models:
        download_and_extract(model, target_data_dir)

