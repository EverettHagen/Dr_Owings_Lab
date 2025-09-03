import pandas as pd
from Bio import Entrez
import re
import time

#example email works fine, keep as is
Entrez.email = "your_email@example.com"
#spreadsheet you want to look at
INPUT_CSV = "output.csv"
#output spreadsheet
OUTPUT_CSV = "enriched_OTUs.csv"
DELAY = 0.4

#number of abstracts you want the computer to look through 
MAX_ABSTRACTS = 20

df = pd.read_csv(INPUT_CSV)

#If the computer can't find relevant info, it will fill the cell with 'unknown'
def fill_default(value, default="Unknown"):
    return value if pd.notna(value) and str(value).strip() else default

#computer will search through PubMed abstracts
def search_pubmed(genus, keywords):
    """Fetch abstracts mentioning genus + keywords"""
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
        time.sleep(DELAY)
        return abstracts_text
    except:
        return ""

def extract_species(genus, text):
    """Extract species names (binomials) from text"""
    matches = re.findall(rf'\b{genus} [a-z]+', text)
    return ", ".join(sorted(set(matches)))

def extract_keywords(text, keyword_dict):
    """Extract relevant keywords from abstracts"""
    found = []
    for label, pattern in keyword_dict.items():
        if re.search(pattern, text, flags=re.IGNORECASE):
            found.append(label)
    return ", ".join(found)

#these are the columns you want the computer to fill
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

# Keyword dictionaries
property_keywords = {
    "Pathogenic": "pathogen OR pathogenic",
    "Antimicrobial Resistant": "antimicrobial resistant OR antibiotic resistant",
    "Culturable": "cultured OR culturable",
    "Plant Pathogen": "plant pathogen OR plant disease",
    "Endosymbiont": "endosymbiont OR symbiont"
}

environment_keywords = {
    "soil": r"\bsoil\b",
    "freshwater": r"\bfreshwater\b|\briver\b|\black\b",
    "marine": r"\bmarine\b|\bocean\b|\bsea\b",
    "plant rhizosphere": r"\brhizosphere\b",
    "human gut": r"\bgut\b|\bintestine\b|\bmicrobiome\b",
    "animal host": r"\banimal\b|\bmouse\b|\brat\b|\bhost\b",
    "extreme environments": r"\bthermophilic\b|\bhalophilic\b|\bextremophile\b"
}

disease_keywords = {
    "sepsis": r"\bsepsis\b",
    "pneumonia": r"\bpneumonia\b",
    "diarrhea": r"\bdiarrhea\b|\bgastroenteritis\b",
    "meningitis": r"\bmeningitis\b",
    "skin infections": r"\bskin infection\b|\bdermatitis\b",
    "respiratory infection": r"\binfection\b.*(lung|respiratory|airway)",
    "general infection": r"\binfection\b"
}

# New spreadsheet columns
df["Pathogenic Species"] = ""
df["Antimicrobial Resistant Species"] = ""
df["Culturable Species"] = ""
df["Plant Pathogenic Species"] = ""
df["Endosymbiont Species"] = ""
df["Environmental Sources"] = ""
df["Human Disease Impact"] = ""

#Loops through each row
for idx, row in df.iterrows():
    genus = row["Genus"]

    # Looks PubMed for genus 
    for prop, kw in property_keywords.items():
        text = search_pubmed(genus, kw)
        species_list = extract_species(genus, text)
        if species_list:
            df.at[idx, f"{prop} Species"] = species_list
            df.at[idx, f"{prop}?"] = "Yes"

    # Environmental associations
    text_env = search_pubmed(genus, "environment OR soil OR water OR gut OR rhizosphere")
    env_hits = extract_keywords(text_env, environment_keywords)
    df.at[idx, "Environmental Sources"] = env_hits if env_hits else "Unknown"

    # Human disease associations
    text_dis = search_pubmed(genus, "disease OR infection OR pathogenic")
    dis_hits = extract_keywords(text_dis, disease_keywords)
    df.at[idx, "Human Disease Impact"] = dis_hits if dis_hits else "None"

# Save to new spreadsheet
df.to_csv(OUTPUT_CSV, index=False)
print(f"CSV with environmental and disease associations saved as {OUTPUT_CSV}")
