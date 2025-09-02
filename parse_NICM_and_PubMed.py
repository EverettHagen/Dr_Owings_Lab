import pandas as pd
from Bio import Entrez
import re
import time

Entrez.email = "your_email@example.com"
INPUT_CSV = "shortened.csv"
OUTPUT_CSV = "enriched_OTUs.csv"
DELAY = 0.3
MAX_ABSTRACTS = 5

df = pd.read_csv(INPUT_CSV)

def fill_default(value, default="Unknown"):
    return value if pd.notna(value) and str(value).strip() else default

def get_species_from_pubmed(genus, keywords):
    """Return a comma-separated list of species in genus associated with keywords."""
    species_set = set()
    try:
        handle = Entrez.esearch(db="pubmed", term=f"{genus}[Organism] AND ({keywords})", retmax=MAX_ABSTRACTS)
        record = Entrez.read(handle)
        handle.close()
        pmids = record["IdList"]
        if not pmids:
            return ""
        handle2 = Entrez.efetch(db="pubmed", id=pmids, rettype="abstract", retmode="text")
        abstracts_text = handle2.read()
        handle2.close()

        # Extract binomial names (Genus species)
        species_matches = re.findall(rf'\b{genus} [a-z]+', abstracts_text)
        for match in species_matches:
            species_set.add(match.strip())
        time.sleep(DELAY)
        return ", ".join(species_set)
    except:
        return ""

# Fill blanks
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

# Keywords for each property
property_keywords = {
    "Pathogenic": "pathogen OR pathogenic",
    "Antimicrobial Resistant": "antimicrobial resistant OR antibiotic resistant",
    "Culturable": "cultured OR culturable",
    "Plant Pathogen": "plant pathogen OR plant disease",
    "Endosymbiont": "endosymbiont OR symbiont"
}

# Create new columns for species names
for prop in property_keywords.keys():
    new_col = f"{prop} Species"
    df[new_col] = ""
    print(f"Fetching species for {prop}...")
    for idx, row in df.iterrows():
        genus = row["Genus"]
        species_list = get_species_from_pubmed(genus, property_keywords[prop])
        df.at[idx, new_col] = species_list
        # Update Yes/No columns if any species found
        if species_list:
            df.at[idx, f"{prop}?"] = "Yes"

# Save enriched CSV
df.to_csv(OUTPUT_CSV, index=False)
print(f"CSV with species-level enrichment saved as {OUTPUT_CSV}")
