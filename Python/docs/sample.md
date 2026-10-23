# GWAS REST API

This sample notebook shows how to answer typical GWAS questions with **pandasGWAS2**. All REST API access — URL construction, pagination, response cleaning and caching — is handled by the package: `get_associations`, `get_studies` and `get_snps` return objects that expose the data as hierarchical pandas DataFrames joined by `association_id` (for associations) or `accession_id` (for studies).

## Libraries

```python
# Install necessary library
!pip install --quiet pandasgwas
```

```python
# Import libraries
import pandas as pd

from pandasgwas import get_associations, get_studies, get_snps

print("Setup complete. pandasGWAS2 is ready to use.")
```

## Sample Questions

Below we have some sample questions and how to answer them using pandasGWAS2.

Note that in all questions below, we use the EFO traits for the diseases and traits. For more information, please visit https://www.ebi.ac.uk/efo/

### Question 1: What variants are associated with a specific disease, e.g. type 2 diabetes?

To find variants associated with "type 2 diabetes mellitus," the code queries the associations endpoint through `get_associations`, using the disease name as a filter. All result pages are fetched automatically. The main `associations` table is then merged with the `snp_effect_alleles` and `mapped_genes` tables on `association_id`, the variant's **rsID** is separated from its **risk allele**, and a clear table with the most important information (variant ID, p-value, and publication) is displayed for the first findings.

```python
# The disease we are interested in
disease_of_interest = "type 2 diabetes mellitus"

# Get all associations for the trait; pagination and caching are handled automatically
assocs = get_associations(efo_trait=disease_of_interest)

print(f"Found {len(assocs)} associations for '{disease_of_interest}'.")
print("Displaying the first 10 results:")

# Combine the main table with the variant and gene tables on association_id
display_df = (
    assocs.associations
    .merge(assocs.snp_effect_alleles, on='association_id')
    .merge(assocs.mapped_genes, on='association_id')
)

# Split 'snp_effect_allele' (e.g. 'rs123-A') into the variant ID and the allele base
display_df[['variant_rsID', 'risk_allele_base']] = display_df['snp_effect_allele'].str.split('-', n=1, expand=True)
display_df['mapped_genes'] = display_df['mapped_genes'].fillna('')

# Select and display the most relevant columns for clarity
print(display_df[[
    'variant_rsID',
    'risk_allele_base',
    'p_value',
    'risk_frequency',
    'pubmed_id',
    'mapped_genes'
]].head(10))
```

### Question 2: What studies with samples including cohort “UKB” are available for “type 2 diabetes”?

This script finds studies on "type 2 diabetes" specifically within the "UK Biobank" (UKB) cohort by using both the `disease_trait` and `cohort` filters of `get_studies`. The cohort information is stored in a separate `cohorts` table, which is merged with the main `studies` table on `accession_id` for display.

```python
# The trait and cohort we are interested in
disease_of_interest = "type 2 diabetes"
cohort_of_interest = "UKB"  # UK Biobank

studies = get_studies(disease_trait=disease_of_interest, cohort=cohort_of_interest)

print(f"Found {len(studies)} studies for '{disease_of_interest}' in the '{cohort_of_interest}' cohort.")
print("Displaying results:")

# Merge the main table with the cohort table on accession_id
studies_df = studies.studies.merge(studies.cohorts, on='accession_id')
studies_df['cohort'] = studies_df['cohort'].apply(
    lambda cohorts: ', '.join(cohorts) if isinstance(cohorts, list) else cohorts
)

# Display the DataFrame with selected, useful columns
print(studies_df[[
    'accession_id',
    'pubmed_id',
    'disease_trait',
    'initial_sample_size',
    'cohort'
]])
```

### Question 3: What are the most significant associations for a specific SNP? eg: significant associations for rs1050316

This script finds the most statistically significant trait associations for a specific SNP, `rs1050316`, using `get_associations(rs_id=...)`. The results are sorted locally by **p-value** in ascending order with pandas, so the top rows are the most significant ones. The trait names are taken from the `efo_traits` table and joined on `association_id` for a clean display.

```python
# The SNP we are interested in
snp_of_interest = "rs1050316"

assocs = get_associations(rs_id=snp_of_interest)

print(f"Found {len(assocs)} associations for SNP '{snp_of_interest}'.")
print("Displaying the 10 most significant results:")

# Sort by p-value (most significant first)
assocs_df = assocs.associations.copy()
assocs_df['p_value'] = pd.to_numeric(assocs_df['p_value'], errors='coerce')

# Join the trait names and keep the top 10 associations
top_df = (
    assocs_df.sort_values('p_value')
    .merge(assocs.efo_traits, on='association_id')
    .drop_duplicates(subset='association_id')
    .head(10)
)

print(top_df[[
    'p_value',
    'efo_trait',
    'risk_frequency',
    'pubmed_id',
    'first_author'
]])
```

### Question 4: Which genes are associated with type 2 diabetes?

This script compiles a comprehensive list of all unique genes associated with "type 2 diabetes mellitus". All associations are fetched with a single `get_associations(efo_trait=...)` call — no manual pagination is needed. The gene names are extracted from the `mapped_genes` table and gathered into a **Python set**, which automatically removes any duplicates.

```python
# The disease we want to find associated genes for
disease_of_interest = "type 2 diabetes mellitus"
print(f"Searching for genes associated with '{disease_of_interest}'...")

# Fetch all associations for the trait in one call
assocs = get_associations(efo_trait=disease_of_interest)

# Use a set to store unique gene names automatically
all_genes = set()
for genes in assocs.mapped_genes['mapped_genes']:
    # A single value can have multiple comma-separated genes
    if isinstance(genes, str):
        for gene in genes.split(','):
            if gene.strip():
                all_genes.add(gene.strip())

# Display the final results
print("\n--- Search Complete ---")
if all_genes:
    sorted_genes = sorted(all_genes)
    print(f"Found a total of {len(sorted_genes)} unique genes associated with '{disease_of_interest}'.")
    print("Example genes:", sorted_genes[:100])
else:
    print(f"Could not find any genes associated with '{disease_of_interest}'.")
```

### Question 5: What are the top significant SNPs for type 2 diabetes in the European population?

This code finds the most significant SNPs for a trait within a specific population. First, it retrieves all associations for "type 2 diabetes mellitus" and all matching studies. The studies' `discovery_ancestries` table is used to collect the publications (PubMed IDs) that involve European-ancestry studies. Associations are then filtered by those PubMed IDs, sorted by p-value, and the top matches are displayed.

```python
# --- Configuration ---
trait_of_interest = "type 2 diabetes mellitus"
population_of_interest = "European"
results_to_find = 10

pd.options.display.float_format = '{:.2e}'.format

print(f"Searching for top {results_to_find} SNPs for '{trait_of_interest}' in '{population_of_interest}' ancestry studies...")

# --- Step 1: Fetch all associations and studies for the trait ---
assocs = get_associations(efo_trait=trait_of_interest)
studies = get_studies(efo_trait=trait_of_interest)
print(f"--- Finished fetching. Found {len(assocs)} total associations. ---\n")

# --- Step 2: Collect publications with European-ancestry studies ---
studies_df = studies.studies.merge(studies.discovery_ancestries, on='accession_id')
is_european = studies_df['discovery_ancestry'].apply(
    lambda ancestry: population_of_interest.lower() in str(ancestry).lower()
)
european_pmids = set(studies_df.loc[is_european, 'pubmed_id'].dropna())

# --- Step 3: Filter and sort the associations ---
assocs_df = assocs.associations.copy()
assocs_df['p_value'] = pd.to_numeric(assocs_df['p_value'], errors='coerce')
top = (
    assocs_df[assocs_df['pubmed_id'].isin(european_pmids)]
    .sort_values('p_value')
    .merge(assocs.snp_effect_alleles, on='association_id')
    .drop_duplicates(subset='association_id')
    .head(results_to_find)
)

# --- Step 4: Display the final results ---
print(f"--- Search Complete ---")
if not top.empty:
    top[['variant_rsID', 'risk_allele_base']] = top['snp_effect_allele'].str.split('-', n=1, expand=True)
    print(f"Found the top {len(top)} matching associations.")
    print(top[[
        'p_value',
        'variant_rsID',
        'risk_allele_base',
        'pubmed_id',
        'first_author'
    ]])
else:
    print("Could not find any matching associations.")

pd.reset_option('display.float_format')
```

### Question 6: Which variants are common between "type 2 diabetes" and "obesity"?

This script identifies genetic variants common to both "type 2 diabetes" and "obesity" by comparing two distinct sets of data. The variants of each trait are collected from the `snp_effect_alleles` table of the corresponding `Association` object and stored in a **set**. Finally, the code uses the `.intersection()` method — a highly efficient operation for sets — to find and display the variants present in both groups.

```python
# Define the two traits of interest
trait_1 = "type 2 diabetes mellitus"
trait_2 = "obesity"

# --- Fetch and collect variants for trait_1 ---
print(f"--- Starting search for '{trait_1}' ---")
assocs_1 = get_associations(efo_trait=trait_1)
variants_for_trait_1 = set(
    assocs_1.snp_effect_alleles['snp_effect_allele']
    .str.split('-').str[0]
    .pipe(lambda s: s[s.str.startswith('rs', na=False)])
    .dropna()
)
print(f"--- Finished fetching for '{trait_1}'. Found {len(variants_for_trait_1)} total unique variants. ---\n")

# --- Fetch and collect variants for trait_2 ---
print(f"--- Starting search for '{trait_2}' ---")
assocs_2 = get_associations(efo_trait=trait_2)
variants_for_trait_2 = set(
    assocs_2.snp_effect_alleles['snp_effect_allele']
    .str.split('-').str[0]
    .pipe(lambda s: s[s.str.startswith('rs', na=False)])
    .dropna()
)
print(f"--- Finished fetching for '{trait_2}'. Found {len(variants_for_trait_2)} total unique variants. ---\n")

# Find the intersection of the two sets
common_variants = variants_for_trait_1.intersection(variants_for_trait_2)

# --- Display the final results ---
print("\n--- Search Complete ---")
if common_variants:
    print(f"Found {len(common_variants)} variants common to both '{trait_1}' and '{trait_2}':")
    print(sorted(common_variants))
else:
    print(f"No common variants found between '{trait_1}' and '{trait_2}'.")
```

### Question 7: What are the most commonly mapped genes for a specific disease? e.g. type 2 diabetes

This script identifies the most frequently mapped genes for "type 2 diabetes mellitus" by fetching all associations for the disease with one `get_associations(efo_trait=...)` call. Every gene name is extracted from the `mapped_genes` table (handling cases where multiple genes are listed together) and gathered into a single list. Pandas then counts the occurrences of each unique gene, and a ranked list of the top 10 most common ones is displayed.

```python
trait_of_interest = "type 2 diabetes mellitus"
results_to_find = 10

print(f"Searching for the most reported genes in associations for '{trait_of_interest}'...")

# --- Step 1: Fetch all associations for the trait ---
assocs = get_associations(efo_trait=trait_of_interest)
print(f"--- Finished fetching. Found {len(assocs)} total associations. ---\n")

# --- Step 2: Extract and count all gene mentions ---
all_gene_mentions = []
for genes in assocs.mapped_genes['mapped_genes']:
    # A single value can have multiple comma-separated genes
    if isinstance(genes, str):
        all_gene_mentions.extend(gene.strip() for gene in genes.split(',') if gene.strip())

# --- Step 3: Display the top N most common genes ---
print(f"\n--- Analysis Complete ---")
if all_gene_mentions:
    gene_counts = pd.Series(all_gene_mentions).value_counts()

    print(f"Found {len(gene_counts)} unique genes in total.")
    print(f"The top {results_to_find} most frequently reported genes for '{trait_of_interest}' are:")
    print(gene_counts.head(results_to_find))
else:
    print(f"Could not find any gene associations for '{trait_of_interest}'.")
```

### Question 8: What traits are associated with a particular SNP?

This script comprehensively finds all unique traits associated with an SNP (`rs1050316`). All associations for the SNP are fetched in one call with `get_associations(rs_id=...)`, and the trait names are read from the `efo_traits` table, whose rows are already one-per-trait. A `set` is used to remove any duplicates before displaying the final, sorted list of traits.

```python
# Configuration
snp_of_interest = "rs1050316"

print(f"Searching for all traits associated with SNP '{snp_of_interest}'...")

# Fetch all associations for the SNP and collect the unique trait names
assocs = get_associations(rs_id=snp_of_interest)
associated_traits = set(assocs.efo_traits['efo_trait'].dropna())

# --- Display the final results ---
print(f"\n--- Search Complete ---")
if associated_traits:
    print(f"Found a total of {len(associated_traits)} unique traits associated with '{snp_of_interest}':")
    for trait in sorted(associated_traits):
        print(f"- {trait}")
else:
    print(f"Could not find any traits associated with '{snp_of_interest}'.")
```

### Question 9: Which studies on a trait have full summary statistics?

This script finds all studies for "type 2 diabetes" that have complete summary statistics by using the `full_pvalue_set=True` filter of `get_studies`. The server returns only matching studies, and all pages are collected automatically into a single `Study` object, which is displayed as a table showing each study's **accession ID** and the **link** to its summary statistics file.

```python
trait_of_interest = "type 2 diabetes mellitus"

print(f"Searching for studies on '{trait_of_interest}' with full summary statistics...")

# Filter for studies with full summary statistics
studies = get_studies(efo_trait=trait_of_interest, full_pvalue_set=True)

# --- Display the final results ---
print(f"\n--- Search Complete ---")
if len(studies):
    print(f"Found a total of {len(studies)} studies for '{trait_of_interest}' with full summary statistics.")
    print(studies.studies[[
        'accession_id',
        'pubmed_id',
        'initial_sample_size',
        'full_summary_stats'  # This column contains the link to the stats
    ]])
else:
    print(f"Could not find any studies for '{trait_of_interest}' with full summary statistics.")
```

### Question 10: What variants are in a gene `HBS1L` and what are their associated traits?

This is a multi-step script that first finds all variants within a specified gene and then identifies the traits associated with those variants. To start, `get_snps(mapped_gene=...)` returns all variants mapped to the gene `HBS1L`. The script then iterates through that list of variants, calling `get_associations(rs_id=...)` for each one to find any linked traits. Thanks to the package's built-in cache, the repeated per-variant requests are only downloaded once. The collected data is structured into a pandas DataFrame mapping each variant to its traits, and "exploded" so that each row shows a single variant-to-trait relationship.

```python
# --- Configuration ---
gene_of_interest = "HBS1L"

print(f"--- Step 1: Finding variants for gene '{gene_of_interest}' ---")

# --- Find all variants for the gene ---
snps = get_snps(mapped_gene=gene_of_interest)
gene_variants_rsids = snps.snps['rs_id'].dropna().tolist()
print(f"\nFound {len(gene_variants_rsids)} total variants for '{gene_of_interest}'.")

# --- Step 2: For each variant, find its associated traits ---
print(f"\n--- Step 2: Checking for associated traits for the {len(gene_variants_rsids)} variants ---")

variant_trait_map = {}
for rs_id in gene_variants_rsids:
    traits_for_snp = set(get_associations(rs_id=rs_id).efo_traits['efo_trait'].dropna())
    variant_trait_map[rs_id] = sorted(traits_for_snp)

# --- Step 3: Display the results ---
print("\n--- Results ---")

display_data = [{'variant_rsID': key, 'associated_traits': value} for key, value in variant_trait_map.items()]
results_df = pd.DataFrame(display_data)

# Explode the list of traits into separate rows for better readability
results_df_exploded = results_df.explode('associated_traits').reset_index(drop=True)

if not results_df_exploded.empty:
    print(results_df_exploded)
else:
    print("Could not find any trait associations for the variants checked.")
```

### Question 11: Which variants are associated with an environmental factor on disease outcome? Eg: variants associated with diet on total blood protein measurement (GxE)

This script performs a complex two-step search to find Gene-by-Environment (GxE) interactions, where a factor like **"diet"** influences a trait like **"total blood protein measurement."** First, `get_studies(efo_trait=..., gxe=True)` gathers all studies flagged as GxE for the trait of interest. In the second step, the code iterates through each of those studies, fetching all their associations with `get_associations(accession_id=...)` and keeping the ones whose EFO traits explicitly mention the **"diet"** environmental factor. Finally, the filtered results are organized into a pandas DataFrame and displayed.

```python
# --- Configuration ---
trait_of_interest = "total blood protein measurement"
environmental_factor = "diet"

print(f"===== Step 1: Finding GxE studies for '{trait_of_interest}' =====")

# --- Find all GxE studies for the trait ---
gxe_studies = get_studies(efo_trait=trait_of_interest, gxe=True)
gxe_accessions = gxe_studies.studies['accession_id'].dropna()
print(f"\nFound {len(gxe_accessions)} GxE studies for '{trait_of_interest}'.")

# --- Step 2: Find associations with the environmental factor within those studies ---
print(f"\n===== Step 2: Searching for '{environmental_factor}' factor within these studies =====")
gxe_associations = []
for accession_id in gxe_accessions:
    assocs = get_associations(accession_id=accession_id)

    # Keep associations whose EFO traits mention the environmental factor
    matched = assocs.efo_traits[
        assocs.efo_traits['efo_trait'].str.contains(environmental_factor, case=False, na=False)
    ]
    for association_id, traits in matched.groupby('association_id')['efo_trait']:
        reported = assocs.reported_traits[
            assocs.reported_traits['association_id'] == association_id
        ]['reported_trait']
        pubmed_id = gxe_studies.studies.loc[
            gxe_studies.studies['accession_id'] == accession_id, 'pubmed_id'
        ]
        gxe_associations.append({
            'accession_id': accession_id,
            'pubmed_id': pubmed_id.iloc[0] if len(pubmed_id) else None,
            'reported_traits': ', '.join(reported.astype(str)),
            'efo_traits': ', '.join(traits.astype(str)),
        })

# --- Step 3: Display the results ---
print(f"\n===== Step 3: Results =====")
if gxe_associations:
    results_df = pd.DataFrame(gxe_associations)
    print(f"Found {len(results_df)} associations matching the GxE criteria.")
    print(results_df[[
        'accession_id',
        'pubmed_id',
        'reported_traits',
        'efo_traits'
    ]])
else:
    print("Could not find any GxE associations matching the specified criteria.")
```

### Question 12: Which genomic region has the most variants associated with a disease? eg: the genomic region with most variants for type 2 diabetes

This script determines which genomic regions are most associated with a disease through a multi-step process. First, it gathers the associations for "type 2 diabetes mellitus" sorted by p-value and collects the unique variants (rsIDs) of the top ones. Next, it iterates through these variants, calling `get_snps(rs_id=...)` for each to find its genomic location — the region name is directly available in the `locations` table under the `region.name` column. Finally, it counts how many times each region appeared and displays a ranked list of the most common regions.

```python
# --- Configuration ---
trait_of_interest = "type 2 diabetes mellitus"
associations_to_check = 100
top_regions_to_display = 10

print(f"--- Step 1: Finding the top {associations_to_check} most significant associations for '{trait_of_interest}' ---")

# --- Step 1: Fetch the associations and collect the unique variants of the top ones ---
assocs = get_associations(efo_trait=trait_of_interest)
assocs_df = assocs.associations.copy()
assocs_df['p_value'] = pd.to_numeric(assocs_df['p_value'], errors='coerce')
top = (
    assocs_df.sort_values('p_value')
    .merge(assocs.snp_effect_alleles, on='association_id')
    .head(associations_to_check)
)

rs_ids = top['snp_effect_allele'].str.split('-').str[0]
unique_rsids = {r for r in rs_ids.dropna() if r.startswith('rs')}
print(f"--- Finished fetching. Found {len(unique_rsids)} unique variants from the top associations. ---\n")

# --- Step 2: For each unique variant, find its genomic region ---
region_list = []
if unique_rsids:
    print(f"\n--- Step 2: Looking up genomic regions for {len(unique_rsids)} variants ---")

    for i, rs_id in enumerate(sorted(unique_rsids)):
        snp = get_snps(rs_id=rs_id)
        region_list.extend(snp.locations['region.name'].dropna())

        # Print progress update
        if (i + 1) % 25 == 0:
            print(f"Processed {i + 1}/{len(unique_rsids)} variants...")

# --- Step 3: Count and display the most common regions ---
print(f"\n--- Analysis Complete ---")
if region_list:
    region_counts = pd.Series(region_list).value_counts()

    print(f"The top {top_regions_to_display} genomic regions with the most variants for '{trait_of_interest}' are:")
    print(region_counts.head(top_regions_to_display))
else:
    print("Could not retrieve region information for the variants.")
```

### Question 13. Which SNP has the strongest effect size/OR for a disease? eg: The SNP with the strongest effect size/OR for type 2 diabetes

This script identifies the single SNP with the strongest effect size, measured by **Odds Ratio** (OR), for "type 2 diabetes mellitus". All associations for the trait are fetched with `get_associations(efo_trait=...)`, sorted locally by `or_value` in descending order, and the top result is placed into a pandas DataFrame for a clean presentation showing the SNP with the highest OR, its **p-value**, and its **rsID**.

```python
# --- Configuration ---
trait_of_interest = "type 2 diabetes mellitus"

print(f"Searching for the SNP with the strongest Odds Ratio for '{trait_of_interest}'...")

# --- Step 1: Fetch the associations and find the highest or_value ---
assocs = get_associations(efo_trait=trait_of_interest)
assocs_df = assocs.associations.copy()
assocs_df['or_value'] = pd.to_numeric(assocs_df['or_value'], errors='coerce')
assocs_df = assocs_df.dropna(subset=['or_value'])

top_association = assocs_df.sort_values('or_value', ascending=False).head(1)
top_association = top_association.merge(assocs.snp_effect_alleles, on='association_id')

# --- Step 2: Display the result ---
if not top_association.empty:
    print("\nFound the following association with the strongest effect size (Odds Ratio):")

    # Parse the risk allele field
    top_association[['variant_rsID', 'risk_allele_base']] = (
        top_association['snp_effect_allele'].str.split('-', n=1, expand=True)
    )

    pd.options.display.float_format = '{:.2e}'.format
    print(top_association[[
        'or_value',
        'p_value',
        'variant_rsID',
        'risk_frequency'
    ]].head(1))
    pd.reset_option('display.float_format')
else:
    print(f"\nCould not find an association with an Odds Ratio for '{trait_of_interest}'.")
```

### Question 14. How many studies report a specific gene associated with a specific disease? eg: number of the studies reporting the KCNQ1 gene associated with type 2 diabetes

This script finds how many studies associate the **KCNQ1** gene with "type 2 diabetes". Instead of fetching all associations and filtering locally, it passes both filters directly to `get_studies(efo_trait=..., mapped_gene=...)`, and the server returns only the matching studies. The total number of unique studies is reported, and their accession IDs are printed.

```python
# --- Configuration ---
trait_of_interest = "type 2 diabetes mellitus"
gene_of_interest = "KCNQ1"

print(f"Searching for the number of studies reporting '{gene_of_interest}' associated with '{trait_of_interest}'...")

# Combine the trait and gene filters in a single query
studies = get_studies(efo_trait=trait_of_interest, mapped_gene=gene_of_interest)

# --- Display the result ---
print(f"\n--- Search Complete ---")
if len(studies):
    print(f"A total of {len(studies)} unique studies report the '{gene_of_interest}' gene in an association with '{trait_of_interest}'.")

    print("\nStudy Accession IDs:")
    for accession_id in sorted(studies.studies['accession_id'].dropna()):
        print(f"- {accession_id}")
else:
    print(f"Could not find any studies reporting '{gene_of_interest}' for '{trait_of_interest}'.")
```

### Question 15. What are the sample sizes used for the top 10 significant variants for disease, eg: find the top 10 significant variants associated with type 2 diabetes and find their sample size

This script answers the question in a two-step process. First, it fetches all associations for "type 2 diabetes mellitus" and keeps the top 10 most significant ones sorted by **p-value**. Then, for each result, it calls `get_studies(pubmed_id=...)` using the association's publication ID to fetch the corresponding study details, including its initial sample size. Finally, a combined table is displayed, linking each top variant to the sample size of the study it was found in.

```python
# --- Configuration ---
trait_of_interest = "type 2 diabetes mellitus"
results_to_find = 10

print(f"--- Step 1: Finding the top {results_to_find} most significant associations for '{trait_of_interest}' ---")

# --- Step 1: Get the top 10 most significant associations for the trait ---
assocs = get_associations(efo_trait=trait_of_interest)
assocs_df = assocs.associations.copy()
assocs_df['p_value'] = pd.to_numeric(assocs_df['p_value'], errors='coerce')
top = assocs_df.sort_values('p_value').head(results_to_find)
print(f"Found {len(top)} associations. Now fetching their study sample sizes...")

# --- Step 2: For each association, find its study's sample size ---
final_results = []
for _, assoc in top.iterrows():
    if assoc['pubmed_id'] is None:
        sample_size = 'N/A'
    else:
        studies = get_studies(pubmed_id=assoc['pubmed_id'])
        sample_size = (
            studies.studies['initial_sample_size'].iloc[0]
            if len(studies.studies) else 'N/A'
        )

    final_results.append({
        'association_id': assoc['association_id'],
        'p_value': assoc['p_value'],
        'initial_sample_size': sample_size,
    })

results_df = pd.DataFrame(final_results).merge(assocs.snp_effect_alleles, on='association_id')
results_df[['variant_rsID', 'risk_allele_base']] = results_df['snp_effect_allele'].str.split('-', n=1, expand=True)

# --- Step 3: Display the results ---
print(f"\n--- Search Complete ---")
if not results_df.empty:
    pd.options.display.float_format = '{:.2e}'.format

    print(results_df[[
        'p_value',
        'variant_rsID',
        'initial_sample_size',  # Display the sample size
    ]].drop_duplicates(subset='association_id'))

    pd.reset_option('display.float_format')
else:
    print(f"Could not retrieve top associations for '{trait_of_interest}'.")
```
