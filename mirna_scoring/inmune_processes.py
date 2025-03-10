import requests
from Bio import Entrez

# Set your email here (required by NCBI)
Entrez.email = "your_email@example.com"


def fetch_go_terms(go_term):
    """Fetch genes associated with a given GO term."""
    url = f"http://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{go_term}/descendants"
    response = requests.get(url, headers={"Accept": "application/json"})

    if response.status_code == 200:
        data = response.json()
        # print (data['results'][0]['descendants'])
        terms = [term for term in data['results'][0]['descendants']]
        terms.append(go_term)  # Include the root term itself
        return terms
    else:
        raise Exception("Failed to fetch GO terms")


def fetch_genes_for_go_terms(go_terms):
    """Fetch genes associated with a list of GO terms."""
    genes = set()
    for term in go_terms:
        handle = Entrez.esearch(db="gene", term=f"{term}[GO]", retmax=1000)
        record = Entrez.read(handle)
        genes.update(record["IdList"])
    return genes


def main():
    go_term = "GO:0002376"  # Immune system process
    # go_terms = fetch_go_terms(go_term)
    genes = fetch_genes_for_go_terms(go_term)

    print(f"Found {len(genes)} genes associated with immune processes:")
    for gene in genes:
        print(gene)
    return genes