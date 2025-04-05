import csv
from typing import List, Dict

def save_results_to_csv(results: List[Dict[str, str]], output_file: str):
    fieldnames = [
        "PubmedID", "Title", "Publication Date",
        "Non-academicAuthor(s)", "CompanyAffiliation(s)",
        "Corresponding Author Email"
    ]
    with open(output_file, mode="w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in results:
            writer.writerow(row)
