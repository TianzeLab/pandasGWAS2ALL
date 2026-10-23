# pandasGWAS2

**pandasGWAS2** is a Python package for programmatic access to the [GWAS Catalog](https://www.ebi.ac.uk/gwas) RESTful API V2. It provides a unified Python interface for retrieving GWAS Catalog data, including associations, studies, traits, variants, publications, genes, ancestries, loci, genomic contexts and summary statistics. It is the Python counterpart of the R package [pandasGWASr](https://github.com/TianzeLab/pandasGWASr).

## Installation

```bash
pip install pandasgwas
```

## Quick start

```python
from pandasgwas import (
    is_API_available,
    get_associations, get_studies, get_traits, get_snps, get_genes,
    get_publications, get_ancestries, get_loci, get_genomic_contexts,
)

# Check that the GWAS Catalog API is reachable
is_API_available()

# Query associations for a study
associations = get_associations(accession_id="GCST000854")

# Query studies, traits, SNPs, genes, publications, ancestries, loci, genomic contexts
studies = get_studies(accession_id="GCST000854")
traits  = get_traits(efo_id="EFO_0001360")
snps    = get_snps(rs_id="rs3093017")
genes   = get_genes(gene_name="HBS1L")
pubs    = get_publications(pubmed_id="35241825")
anc     = get_ancestries(accession_id="GCST000854")
loci    = get_loci(association_id="169106883")
gc      = get_genomic_contexts(rs_id="rs3093017")
```

## Data retrieval functions

Nine high-level functions cover every GWAS Catalog entity. Each accepts filter parameters, a `cached` flag (default `True`) and a `mode` argument (`"Fat"` or `"Thin"`).

| Function | Entity | Key filters |
|---|---|---|
| `get_associations()` | Association | `association_id`, `pubmed_id`, `rs_id`, `accession_id`, `efo_trait`, `efo_id`, `mapped_gene`, … |
| `get_studies()` | Study | `accession_id`, `pubmed_id`, `disease_trait`, `efo_trait`, `cohort`, `gxe`, `ancestral_group`, … |
| `get_traits()` | Trait | `efo_trait`, `uri`, `efo_id`, `pubmed_id`, `mapped_gene`, … |
| `get_snps()` | Snp | `rs_id`, `bp_location`, `bp_start`, `bp_end`, `chromosome`, `pubmed_id`, … |
| `get_genes()` | Gene | `gene_name` (optional — omit to list all genes) |
| `get_publications()` | Publication | `pubmed_id`, `title`, `first_author` |
| `get_ancestries()` | Ancestry | `accession_id` (required), `ancestry_id` (optional sub-id) |
| `get_loci()` | Locus | `association_id` (required), `locus_id` (optional sub-id) |
| `get_genomic_contexts()` | GenomicContext | `rs_id` (required), `genomicContext_id` (optional sub-id) |

`get_associations()` is unique in supporting both an `association_id` path lookup and query filters in the same call — when both are supplied the intersection is returned.

## Data objects

Each query function returns an instance of the corresponding class: `Association`, `Study`, `Trait`, `Snp`, `Gene`, `Publication`, `Ancestry`, `Locus`, `GenomicContext`.

### Fat and Thin modes

`mode="Fat"` (default) builds the main DataFrame together with nested DataFrames (e.g. `efo_traits`, `locations`, `mapped_genes`). `mode="Thin"` keeps only the raw records in `.raw_data`:

```python
thin = get_associations(association_id="169106883", mode="Thin")
```

### Attributes

Fat-mode objects expose one or more DataFrame attributes. For example, a `Study` object contains `studies`, `efo_traits`, `genotyping_technologies`, `discovery_ancestries`, `replication_ancestries`, `bg_efo_traits`, `cohorts` and `array_manufacturers`.

```python
studies = get_studies(accession_id="GCST000854")
studies.studies
studies.efo_traits
```

### Indexing

All data objects support position and identifier indexing (Python is 0-based). Identifier keys vary by class — `accession_id` for Study, `rs_id` for Snp, `efo_id` for Trait, composite keys (`"accession_id|ancestry_id"`, `"association_id|locus_id"`, `"rs_id|genomicContext_id"`) for Ancestry, Locus and GenomicContext.

```python
associations[0]
associations["169106883"]
associations[0:2]

study = get_studies(accession_id="GCST000854")
study["GCST000854"]
```

### Set operations

Combine queried results using operators or the named wrapper functions. `bind` (`+`) preserves duplicates; `union` (`|`) does not.

| Operator | Function | Behaviour |
|---|---|---|
| `+` | `bind()` | concatenation (duplicates kept) |
| `&` | `intersect()` | intersection (no duplicates) |
| `-` | `set_diff()` | difference (no duplicates) |
| `\|` | `union()` | union (no duplicates) |
| `^` | `set_xor()` | symmetric difference (no duplicates) |
| `==` | `set_equal()` | equality test |

```python
from pandasgwas import get_snps, bind, intersect, set_diff, union, set_xor, set_equal

a = get_snps(rs_id="rs3093017")
b = get_snps(rs_id="rs7671745")

bind(a, b)       # or a + b
intersect(a, b)  # or a & b
set_diff(a, b)   # or a - b
union(a, b)      # or a | b
set_xor(a, b)    # or a ^ b
set_equal(a, b)  # or a == b
```

## Browser helpers

Open relevant database pages in the default web browser:

```python
from pandasgwas import (
    open_in_pubmed, open_in_dbsnp, open_in_gtex,
    open_study_in_gwas_catalog, open_variant_in_gwas_catalog,
    open_trait_in_gwas_catalog, open_gene_in_gwas_catalog,
    open_region_in_gwas_catalog, open_publication_in_gwas_catalog,
)

open_in_pubmed("26301688")
open_in_dbsnp("rs56261590")
open_in_gtex("rs56261590")
open_study_in_gwas_catalog("GCST000016")
open_variant_in_gwas_catalog("rs146992477")
open_trait_in_gwas_catalog("EFO_0004884")
open_gene_in_gwas_catalog("DPP6")
open_region_in_gwas_catalog("2q37.1")
open_publication_in_gwas_catalog("25533513")
```

## Summary statistics submodule

The `pandasgwas.summary_statistics` module provides FTP-based retrieval of harmonised summary statistics:

```python
from pandasgwas import summary_statistics

# Search the index by PubMed ID, study accession, or EFO trait ID
search_DF = summary_statistics.search(
    PubMed_id="27918534",
    study_accession_id="GCST003966",
)

# Browse the FTP directory in your web browser
summary_statistics.browser(search_DF)

# Download files to ~/pandasgwas_home
summary_statistics.download(search_DF)

# Load downloaded files into a DataFrame
df = summary_statistics.parse(search_DF)
```

`search()` reads a bundled harmonised index (`harmonised_list.txt`) by default. Pass `online_index=True` to download the latest index from the FTP server first.

## Cache management

All API responses are cached in memory for 24 hours. Two functions control the cache:

```python
from pandasgwas import clear_cache, reinit_cache

# Clear all caches (or pass a specific entity name)
clear_cache()
clear_cache("Association")

# Reinitialise a specific cache or all caches with custom size/ttl
reinit_cache(field="Association", maxsize=4 * 1024, ttl=60 * 60)
reinit_cache(field="All")
```

## API availability

```python
from pandasgwas import is_API_available

is_API_available()  # True if the GWAS Catalog API responds with HTTP 200
```

## Documentation

See [pandasGWAS Documentation](https://caotianze.github.io/pandasgwas/)

## Licensing information

### Source code

This package is licensed under the MIT License.

### Data from NHGRI-EBI GWAS Catalog

The NHGRI-EBI GWAS Catalog and all its contents are available under the general [Terms of Use for EMBL-EBI Services](https://www.ebi.ac.uk/about/terms-of-use/). Summary statistics are made available under [CC0](https://creativecommons.org/publicdomain/zero/1.0/) unless [otherwise stated](https://www.ebi.ac.uk/gwas/docs/faq/#faq-H7).

## Similar projects

- R package [pandasGWASr](https://github.com/TianzeLab/pandasGWASr)
- R package [gwasrapidd](https://github.com/ramiromagno/gwasrapidd) by Ramiro Magno
