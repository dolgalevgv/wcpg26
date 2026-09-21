#!/usr/bin/env python

import sys
from pathlib import Path

import pandas as pd


donors_path, donors_vcf_path = [Path(path) for path in sys.argv[1:3]]
donors = pd.read_csv(donors_path, dtype=str, usecols=["vcf_id", "donor_id"])

donors_vcf = donors[["vcf_id", "donor_id"]]

donors_vcf.to_csv(donors_vcf_path, index=False, header=False, sep=" ")
