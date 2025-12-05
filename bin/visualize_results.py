#! /usr/bin/env python3

"""
After ensembling data, one final cell annotation remains, which can now be visualized.
In addition, having metadata for the cells from the single cell & atac seq data, more detailed visualizations can be created.
Furthermore, the ensembled cell annotation can be compared to the annotation of the individual tools.
# TODO: Load all the data from the different steps
# TODO: Combine data when necessary
# TODO: Visualize the results
"""

import pandas as pd
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import muon as mu
import tangram as tg
import argparse as ap



if __name__ == "__main__":
    # DUMMY
    parser = ap.ArgumentParser(description= "Here, the results will be visualized.")
    parser.add_argument()

    args = parser.parse_args()