"""Downloads and extracts UniBind data for select celltypes, combines into sites file."""

from json import load
from urllib import request
from collections import OrderedDict
from coreapi import Client

# Globals
PROFILE = snakemake.params.profile  # type: ignore
INSTALL_DIR = snakemake.params.install_dir  # type: ignore

###
# Functions
###


def download_meta(profile: str = PROFILE) -> OrderedDict:
    """Return UniBind API metadata for provided JASPAR profile"""
    # Initialize a client & load the schema document
    client = Client()
    schema = client.get("https://unibind.uio.no/api/v1/docs")
    action = ["datasets", "list"]

    # For example
    params = {
        "collection": "Robust",
        "page_size": 1000000,
        "jaspar_id": profile,
        "species": "Homo sapiens",
    }
    # Return profile metadata
    return client.action(schema, action, params=params)


def download_damos(metadata: OrderedDict) -> None:
    """Downloads all DAMO bedfiles from proivded metadata"""
    # URLs for download from metadata
    result_urls = [result["url"] for result in metadata["results"]]

    # Process data from links, need to pull out the bedfile
    for result_url in result_urls:
        with request.urlopen(result_url) as url:
            data = load(url)
            damo = data["tfbs"][0]["DAMO"][0]["bed_url"]
            name = damo.split("/")[-1]
            request.urlretrieve(damo, f"{INSTALL_DIR}/{name}")


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    # Download
    metadata = download_meta()
    download_damos(metadata)
