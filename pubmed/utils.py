import re

ACADEMIC_KEYWORDS = ["university", "college", "institute", "school", "department", "centre", "center", "hospital", "clinic", "faculty", "lab"]
PHARMA_KEYWORDS = ["pharma", "biotech", "therapeutics", "inc", "ltd", "llc", "gmbh", "co.", "corporation", "industries"]

def extract_email(text: str):
    match = re.search(r"[\w\.-]+@[\w\.-]+", text)
    return match.group(0) if match else None

def is_non_academic_affiliation(affiliation: str) -> bool:
    return not any(keyword in affiliation.lower() for keyword in ACADEMIC_KEYWORDS)

def filter_pharma_affiliation(affiliation: str) -> bool:
    return any(keyword in affiliation.lower() for keyword in PHARMA_KEYWORDS)
