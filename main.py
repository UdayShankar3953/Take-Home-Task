import argparse
from time import sleep
from pubmed.fetch import fetch_pubmed_ids, fetch_pubmed_details
from pubmed.parser import extract_paper_data
from output.writer import save_results_to_csv

def main():
    parser = argparse.ArgumentParser(description="Fetch non-academic PubMed papers")
    parser.add_argument("--query", required=True, help="PubMed search query")
    parser.add_argument("--max-results", type=int, default=10)
    parser.add_argument("--file", default="")
    parser.add_argument("--debug", action="store_true")
    args = parser.parse_args()

    pubmed_ids = fetch_pubmed_ids(args.query, args.max_results)
    if args.debug:
        print(f"[DEBUG] Found {len(pubmed_ids)} PubMed IDs.")

    papers = fetch_pubmed_details(pubmed_ids)
    results = []

    for article in papers:
        sleep(0.3)
        data = extract_paper_data(article)
        results.append(data)
        if args.debug and data["CompanyAffiliation(s)"] != "N/A":
            print(f"[DEBUG] {data['Title'][:60]}... | Companies: {data['CompanyAffiliation(s)']}")

    if args.file:
        save_results_to_csv(results, args.file)
        print(f"Saved {len(results)} results to {args.file}")
    else:
        for result in results:
            print(result)

if __name__ == "__main__":
    main()
