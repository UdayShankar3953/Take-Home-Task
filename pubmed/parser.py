import re
from typing import Dict, Any, List, Optional
from .utils import is_non_academic_affiliation, filter_pharma_affiliation, extract_email

def extract_paper_data(article: Dict[str, Any]) -> Dict[str, Any]:
    medline = article["MedlineCitation"]
    article_data = medline["Article"]

    pubmed_id = str(medline.get("PMID", "N/A"))
    title = article_data.get("ArticleTitle", "N/A")

    pub_date_elem = article_data.get("Journal", {}).get("JournalIssue", {}).get("PubDate", {})
    pub_date = f"{pub_date_elem.get('Year', '')}-{pub_date_elem.get('Month', '')}-{pub_date_elem.get('Day', '')}".strip("-")
    if not pub_date:
        pub_date = "N/A"

    non_academic_authors, company_affiliations, corresponding_email = [], [], ""
    authors = article_data.get("AuthorList", [])

    for author in authors:
        if "AffiliationInfo" in author and author["AffiliationInfo"]:
            affiliation = author["AffiliationInfo"][0].get("Affiliation", "")
            if is_non_academic_affiliation(affiliation):
                full_name = f"{author.get('ForeName', '').strip()} {author.get('LastName', '').strip()}".strip()
                if full_name:
                    non_academic_authors.append(full_name)
                if filter_pharma_affiliation(affiliation):
                    company_affiliations.append(affiliation)
                if not corresponding_email:
                    email = extract_email(affiliation)
                    if email:
                        corresponding_email = email

    return {
        "PubmedID": pubmed_id,
        "Title": title,
        "Publication Date": pub_date,
        "Non-academicAuthor(s)": "; ".join(non_academic_authors) or "N/A",
        "CompanyAffiliation(s)": "; ".join(company_affiliations) or "N/A",
        "Corresponding Author Email": corresponding_email or "N/A"
    }
