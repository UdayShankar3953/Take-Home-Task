
# 🧬 PubMed Research Paper Fetcher

A modular Python CLI tool that fetches biomedical research papers from PubMed based on a search query and filters them to include only those with at least one author affiliated with a pharmaceutical or biotech company.



## 🚀 Features

- 🔍 Search PubMed articles via the Entrez API (Biopython)
- 🧪 Filter out academic authors using keyword heuristics
- 🏢 Identify papers with non-academic (pharma/biotech) contributors
- 📤 Export structured results to CSV
- ⚙️ Command-line interface with useful flags
- 🪄 Debug mode for step-by-step trace



## 📁 Project Structure

```
pubmed_fetcher/
├── main.py              # CLI entry point
├── config.py            # Entrez API config (email)
├── pubmed/
│   ├── fetch.py         # ESearch and EFetch logic
│   ├── parser.py        # XML parsing and author filtering
│   └── utils.py         # Helper functions (email, heuristics)
└── output/
    └── writer.py        # CSV writer
```



## 🛠️ Installation

> Requires Python 3.7+

1. Clone the repo:

```bash
git clone https://github.com/udayshankar3953/pubmed-fetcher.git
cd pubmed-fetcher
```

2. Create a virtual environment (optional but recommended):

```bash
python -m venv venv
source venv/bin/activate  # or venv\Scripts\activate on Windows
```

3. Install dependencies:

```bash
pip install -r requirements.txt
```

---

## ⚙️ Usage

```bash
python main.py --query "<search term>" --max-results <number> [--file output.csv] [--debug]
```

### Arguments:

| Flag            | Description                                      |
|-----------------|--------------------------------------------------|
| `--query`       | PubMed search query (e.g., `"cancer AND AI"`)   |
| `--max-results` | Number of results to fetch (e.g., `20`)         |
| `--file`        | (Optional) CSV output filename                  |
| `--debug`       | (Optional) Enable verbose debug logging         |

### Example:

```bash
python main.py --query "cancer AND chemotherapy" --max-results 10 --file results.csv --debug
```

---

## 🧠 Filtering Logic

- **Academic Affiliation (excluded):** university, institute, hospital, school, clinic, lab
- **Pharma/Biotech Affiliation (included):** pharma, biotech, inc, ltd, therapeutics, gmbh, llc

Any paper with at least one non-academic author is included in the final result.

---

## 📝 Sample Output

| PubmedID  | Title                     | Non-academicAuthor(s) | CompanyAffiliation(s) | Email              |
|-----------|---------------------------|------------------------|------------------------|--------------------|
| 12345678  | Sample Title              | Jane Doe               | Acme Biotech Inc.      | jane@acmebio.com   |

---

## 📚 Dependencies

- [Biopython](https://biopython.org/)
- Python standard libraries (`argparse`, `csv`, `xml`, `re`, `time`)

Install them with:

```bash
pip install biopython
```


## 📽️ Demo

Watch the demo video 👉 [Insert video link here]



## 👤 Author

**Udyavara Uday Shankar**  
📧 udayshankar.udyavara@gmail.com  
🔗 [LinkedIn](http://www.linkedin.com/in/udyavaraudayshankar/)  
💻 [GitHub](https://github.com/udayshankar3953)

---

## 🙌 Acknowledgements

- Thanks to Aganitha AI for the opportunity.
