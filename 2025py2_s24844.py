#!/usr/bin/env python3

from Bio import Entrez, SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import sys

class NCBIRetriever:
    def __init__(self, email, api_key, min_len, max_len):
        Entrez.email = email
        Entrez.api_key = api_key
        self.min_len = min_len
        self.max_len = max_len

    def search(self, taxid):
        term = f"txid{taxid}[Organism]"
        handle = Entrez.esearch(db="nucleotide", term=term, usehistory="y", retmax=0)
        result = Entrez.read(handle)
        self.webenv = result["WebEnv"]
        self.query_key = result["QueryKey"]
        return int(result["Count"])

    def fetch(self, total=1000):
        records = []
        for start in range(0, total, 100):
            handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text",
                                   retstart=start, retmax=100,
                                   webenv=self.webenv, query_key=self.query_key)
            for record in SeqIO.parse(handle, "genbank"):
                seq_len = len(record.seq)
                if self.min_len <= seq_len <= self.max_len:
                    records.append({
                        "accession": record.id,
                        "length": seq_len,
                        "description": record.description
                    })
            handle.close()
        return records

def save_csv(data, filename):
    df = pd.DataFrame(data)
    df.to_csv(filename, index=False)

def plot_lengths(data, filename):
    df = pd.DataFrame(data).sort_values(by="length", ascending=False)
    plt.figure(figsize=(12, 6))
    plt.plot(df["accession"], df["length"], marker="o")
    plt.xticks(rotation=90, fontsize=6)
    plt.xlabel("Accession")
    plt.ylabel("Sequence Length")
    plt.title("Sequence Lengths Sorted by Size")
    plt.tight_layout()
    plt.savefig(filename)

def main():
    email = input("Enter your email: ")
    api_key = input("Enter your NCBI API key: ")
    taxid = input("Enter taxonomic ID: ")
    min_len = int(input("Min sequence length: "))
    max_len = int(input("Max sequence length: "))

    retriever = NCBIRetriever(email, api_key, min_len, max_len)
    total = retriever.search(taxid)
    print(f"Total records found: {total}. Downloading and filtering...")

    data = retriever.fetch(total=min(1000, total))
    print(f"Filtered records: {len(data)}")

    csv_file = f"taxid_{taxid}_filtered.csv"
    png_file = f"taxid_{taxid}_plot.png"
    save_csv(data, csv_file)
    plot_lengths(data, png_file)

    print(f"Saved CSV to {csv_file}")
    print(f"Saved plot to {png_file}")

if __name__ == "__main__":
    main()
