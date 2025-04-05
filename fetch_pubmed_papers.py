import argparse
import csv
import re
from time import sleep
from typing import List, Dict, Any, Optional
from Bio import Entrez

# Set your email here (required by Entrez)
Entrez.email = "212g1a3953@gmail.com"  # Replace with your actual email

# Define keywords for non-academic (pharmaceutical/biotech) affiliation detection.
ACADEMIC_KEYWORDS = ["university", "college", "institute", "school", "department", "centre", "center", "hospital", "clinic", "faculty", "lab"]
PHARMA_KEYWORDS = ["pharma", "biotech", "therapeutics", "inc", "ltd", "llc", "gmbh", "co.", "corporation", "industries"]

def fetch_pubmed_ids(query: str, max_results: int) -> List[str]:
    """Fetches PubMed IDs for the given query using PubMed's ESearch API."""
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]

def fetch_pubmed_details(pubmed_ids: List[str]) -> List[Dict[str, Any]]:
    """Fetches detailed information about the given PubMed IDs using EFetch."""
    ids = ",".join(pubmed_ids)
    handle = Entrez.efetch(db="pubmed", id=ids, rettype="medline", retmode="xml")
    records = Entrez.read(handle)
    handle.close()
    return records["PubmedArticle"]

def extract_email(text: str) -> Optional[str]:
    """Extracts an email address from a string using regex."""
    match = re.search(r"[\w\.-]+@[\w\.-]+", text)
    return match.group(0) if match else None

def is_non_academic_affiliation(affiliation: str) -> bool:
    """Determines if an affiliation is non-academic by checking for academic keywords."""
    # If any academic keyword is present, consider it academic.
    return not any(keyword.lower() in affiliation.lower() for keyword in ACADEMIC_KEYWORDS)

def filter_pharma_affiliation(affiliation: str) -> bool:
    """Checks if the affiliation string contains pharma/biotech related keywords."""
    return any(keyword.lower() in affiliation.lower() for keyword in PHARMA_KEYWORDS)

def extract_paper_data(article: Dict[str, Any]) -> Dict[str, Any]:
    """
    Extracts required fields from a PubMedArticle XML element.
    Returns a dictionary with:
      - PubmedID
      - Title
      - Publication Date
      - Non-academicAuthor(s)
      - CompanyAffiliation(s)
      - Corresponding Author Email
    """
    medline = article["MedlineCitation"]
    article_data = medline["Article"]

    # PubmedID
    pubmed_id = str(medline.get("PMID", "N/A"))

    # Title
    title = article_data.get("ArticleTitle", "N/A")

    # Publication Date: Try to extract date from DateCreated; fallback if necessary.
    pub_date_elem = article_data.get("Journal", {}).get("JournalIssue", {}).get("PubDate", {})
    pub_date = f"{pub_date_elem.get('Year', '')}-{pub_date_elem.get('Month', '')}-{pub_date_elem.get('Day', '')}".strip("-")
    if not pub_date:
        pub_date = "N/A"

    non_academic_authors: List[str] = []
    company_affiliations: List[str] = []
    corresponding_email = ""

    authors = article_data.get("AuthorList", [])
    for author in authors:
        # Check if author has affiliation information.
        if "AffiliationInfo" in author and author["AffiliationInfo"]:
            affiliation = author["AffiliationInfo"][0].get("Affiliation", "")
            # Heuristic: if affiliation is non-academic, then include author.
            if is_non_academic_affiliation(affiliation):
                # Build full name if possible.
                full_name = f"{author.get('ForeName', '').strip()} {author.get('LastName', '').strip()}".strip()
                if full_name:
                    non_academic_authors.append(full_name)
                # If affiliation contains pharma/biotech keywords, add to company affiliations.
                if filter_pharma_affiliation(affiliation):
                    company_affiliations.append(affiliation)
                # Extract email if not already found.
                if not corresponding_email:
                    email = extract_email(affiliation)
                    if email:
                        corresponding_email = email

    return {
        "PubmedID": pubmed_id,
        "Title": title,
        "Publication Date": pub_date,
        "Non-academicAuthor(s)": "; ".join(non_academic_authors) if non_academic_authors else "N/A",
        "CompanyAffiliation(s)": "; ".join(company_affiliations) if company_affiliations else "N/A",
        "Corresponding Author Email": corresponding_email if corresponding_email else "N/A"
    }

def save_results_to_csv(results: List[Dict[str, Any]], output_file: str) -> None:
    """Saves the list of result dictionaries to a CSV file."""
    fieldnames = [
        "PubmedID", "Title", "Publication Date",
        "Non-academicAuthor(s)", "CompanyAffiliation(s)",
        "Corresponding Author Email"
    ]
    with open(output_file, mode="w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for result in results:
            writer.writerow(result)

def main():
    parser = argparse.ArgumentParser(
        description="Fetch research papers from PubMed with non-academic (pharma/biotech) affiliations."
    )
    parser.add_argument("--query", required=True, help="Search query for PubMed (supports full query syntax)")
    parser.add_argument("--max-results", type=int, default=10, help="Maximum number of papers to fetch")
    parser.add_argument("--file", default="", help="Output CSV file (if omitted, prints to console)")
    parser.add_argument("--debug", action="store_true", help="Print debug information during execution")
    args = parser.parse_args()

    if args.debug:
        print(f"[DEBUG] Query: {args.query}, Max Results: {args.max_results}")

    # Fetch PubMed IDs
    pubmed_ids = fetch_pubmed_ids(args.query, args.max_results)
    if args.debug:
        print(f"[DEBUG] Found {len(pubmed_ids)} PubMed IDs.")

    # Fetch paper details
    papers = fetch_pubmed_details(pubmed_ids)
    if args.debug:
        print(f"[DEBUG] Fetched details for {len(papers)} papers.")

    # Extract and filter paper data
    results = []
    for article in papers:
        sleep(0.3)  # To avoid overloading the API
        data = extract_paper_data(article)
        results.append(data)
        if args.debug and data["CompanyAffiliation(s)"] != "N/A":
            print(f"[DEBUG] {data['Title'][:60]}... | Companies: {data['CompanyAffiliation(s)']}")

    # If output file is specified, save CSV; otherwise, print the results.
    if args.file:
        save_results_to_csv(results, args.file)
        print(f"Saved {len(results)} results to {args.file}")
    else:
        for result in results:
            print(result)

if __name__ == "__main__":
    main()
