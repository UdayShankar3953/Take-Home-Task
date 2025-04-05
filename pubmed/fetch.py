from typing import List, Dict, Any
from Bio import Entrez
from config import Entrez_email

Entrez.email = Entrez_email

def fetch_pubmed_ids(query: str, max_results: int) -> List[str]:
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]

def fetch_pubmed_details(pubmed_ids: List[str]) -> List[Dict[str, Any]]:
    ids = ",".join(pubmed_ids)
    handle = Entrez.efetch(db="pubmed", id=ids, rettype="medline", retmode="xml")
    records = Entrez.read(handle)
    handle.close()
    return records["PubmedArticle"]
