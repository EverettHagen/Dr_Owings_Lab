import pandas as pd
from Bio import Entrez
import time
import re

# === Settings ===
Entrez.email = "your_email@example.com"  # Required by NCBI
INPUT_CSV = "shortened.csv"
OUTPUT_CSV = "enriched_OTUs_full.csv"
DELAY = 0.3  # seconds between NCBI requests
MAX_ABSTRACTS = 5  # abstracts per genus for disease extraction

# === Load CSV ===
df = pd.read_csv(INPUT_CSV)

# === Helper functions ===
def fill_default(value, default="Unknown"):
    return value if pd.notna(value) and str(value).strip() else default

def get_ncbi_lineage(genus):
    try:
        handle = Entrez.esearch(db="taxonomy", term=genus)
        record = Entrez.read(handle)
        handle.close()
        if record["IdList"]:
            tax_id = record["IdList"][0]
            handle2 = Entrez.efetch(db="taxonomy", id=tax_id)
            tax_info = Entrez.read(handle2)
            handle2.close()
            lineage = tax_info[0]["Lineage"]
            scientific_name = tax_info[0]["ScientificName"]
            return f"{scientific_name} ({lineage})"
        else:
            return "Unknown"
    except:
        return "Error"

def pubmed_search_count(genus, keyword):
    """Return True if PubMed has articles for genus + keyword."""
    try:
        handle = Entrez.esearch(db="pubmed", term=f"{genus}[Organism] AND {keyword}", retmax=1)
        record = Entrez.read(handle)
        handle.close()
        return int(record["Count"]) > 0
    except:
        return False

def extract_diseases_from_abstracts(genus):
    """Return comma-separated diseases found in PubMed abstracts."""
    diseases_found = set()
    try:
        handle = Entrez.esearch(db="pubmed", term=f"{genus}[Organism] AND disease", retmax=MAX_ABSTRACTS)
        record = Entrez.read(handle)
        handle.close()
        pmids = record["IdList"]
        if not pmids:
            return ""
        
        handle2 = Entrez.efetch(db="pubmed", id=pmids, rettype="abstract", retmode="text")
        abstracts_text = handle2.read()
        handle2.close()
        
        disease_keywords = [
            "asthma", "diabetes", "cancer", "periodontitis", "cholera",
            "tuberculosis", "food poisoning", "gastroenteritis", "infection",
            "chronic granulomatous disease", "COPD"
        ]
        for keyword in disease_keywords:
            if re.search(rf"\b{keyword}\b", abstracts_text, re.IGNORECASE):
                diseases_found.add(keyword.title())
        
        time.sleep(DELAY)
        return ", ".join(diseases_found)
    
    except:
        return ""

# === Fill defaults for blank columns ===
columns_to_fill = [
    "Pathogenic Species/Strains in Genus?",
    "Culturable?",
    "Environmental Associations",
    "Genus Contains Antimicrobial Resistant Strains?",
    "Associated with Food-borne Illness?",
    "Disease Associations",
    "Plant Pathogen?",
    "Endosymbiont?"
]

for col in columns_to_fill:
    df[col] = df[col].apply(lambda x: fill_default(x))

# === Enrich taxonomy ===
genera = df["Genus"].unique()
print("Fetching NCBI taxonomy lineage...")
taxonomy_dict = {genus: get_ncbi_lineage(genus) for genus in genera}
df["NCBI_Taxonomy"] = df["Genus"].map(taxonomy_dict)

# === PubMed enrichment ===
print("Performing PubMed searches for all relevant columns...")

search_columns = {
    "Pathogenic Species/Strains in Genus?": "pathogen OR pathogenic",
    "Disease Associations": "disease",
    "Environmental Associations": "environment OR soil OR water",
    "Associated with Food-borne Illness?": "food OR fermentation",
    "Genus Contains Antimicrobial Resistant Strains?": "antimicrobial resistant OR antibiotic resistant",
    "Culturable?": "cultured OR culturable",
    "Plant Pathogen?": "plant pathogen OR plant disease",
    "Endosymbiont?": "endosymbiont OR symbiont"
}

# Perform searches
search_results = {}
for col, keyword in search_columns.items():
    print(f"Searching for {col}...")
    search_results[col] = {genus: pubmed_search_count(genus, keyword) for genus in genera}

# Map results back to dataframe
for col in search_columns.keys():
    df[col] = df.apply(
        lambda row: "Yes" if search_results[col].get(row["Genus"], False) else row[col],
        axis=1
    )

# === Extract associated diseases ===
print("Extracting disease names from PubMed abstracts...")
df["Associated Diseases"] = ""
for idx, row in df.iterrows():
    if row["Disease Associations"] == "Yes":
        diseases = extract_diseases_from_abstracts(row["Genus"])
        df.at[idx, "Associated Diseases"] = diseases

# === Save enriched CSV ===
df.to_csv(OUTPUT_CSV, index=False)
print(f"Fully enriched CSV saved to {OUTPUT_CSV}")
