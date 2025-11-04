"""
This file is for IDEs and linters to understand the `snakemake` object.
It should not be executed by Snakemake.
"""
from typing import List, Dict

class Snakemake:
    input: List[str]
    output: List[str]
    wildcards: Dict[str, str]
    config: Dict
    params: Dict
    log: List[str]
    threads: int

# Create a dummy snakemake object
snakemake = Snakemake()