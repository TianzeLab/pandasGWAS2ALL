# pandasGWAS2ALL

> A combined mirror repository for the Python package **pandasGWAS2** and its R counterpart **pandasGWASr**, automatically synced from both upstream repositories via GitHub Actions.

## About

[pandasGWAS2](https://github.com/TianzeLab/pandasGWAS2) (Python) and [pandasGWASr](https://github.com/TianzeLab/pandasGWASr) (R) are feature-identical packages that provide programmatic access to the [GWAS Catalog](https://www.ebi.ac.uk/gwas) RESTful API V2. Both cover the same entities — Association, Study, Trait, Snp, Gene, Publication, Ancestry, Locus, and GenomicContext — with matching capabilities (Fat/Thin modes, set operations, browser helpers, summary statistics submodule, cache management), and are implemented line-by-line in alignment.

This repository brings both packages together for easy side-by-side reference, cross-language testing, and joint usage.

## Repository layout

```
pandasGWAS2ALL/
├── Python/               # pandasGWAS2 — Python package source
│   ├── pandasgwas/       # Main package module
│   ├── tests/
│   ├── docs/
│   ├── setup.py
│   └── pyproject.toml
├── R/                    # pandasGWASr — R package source
│   ├── R/
│   ├── man/
│   ├── tests/
│   ├── DESCRIPTION
│   └── NAMESPACE
└── .github/workflows/
    └── sync-upstream.yml # Manual upstream sync workflow
```

## Quick start

### Python — pandasGWAS2

```bash
pip install pandasgwas
```

```python
from pandasgwas import get_associations, get_studies, get_traits

associations = get_associations(accession_id="GCST000854")
studies      = get_studies(accession_id="GCST000854")
traits       = get_traits(efo_id="EFO_0001360")
```

See [Python/README.md](Python/README.md)for full usage.

### R — pandasGWASr

```r
remotes::install_github("TianzeLab/pandasGWASr")
```

```r
library(pandasGWASr)

associations <- get_associations(accession_id = "GCST000854")
studies      <- get_studies(accession_id = "GCST000854")
traits       <- get_traits(efo_id = "EFO_0001360")
```

See [R/README.md](R/README.md) for full usage.

## Sync mechanism

This repository stays in sync with its upstreams through the **Sync Upstream Repos** workflow defined in `.github/workflows/sync-upstream.yml`:

| Directory | Upstream repo | Branch |
|---|---|---|
| `R/` | `TianzeLab/pandasGWASr` | `master` |
| `Python/` | `TianzeLab/pandasGWAS2` | `main` |

Flow: clear target directory → shallow-clone upstream into a temp dir (excluding `.git`) → `rsync` into the target → auto commit + push. Trigger it from **Actions → Sync Upstream Repos → Run workflow** in the GitHub UI.

## Feature parity

The two packages are API-for-API equivalent:

| Capability | pandasGWAS2 (Python) | pandasGWASr (R) |
|---|---|---|
| Installation | `pip install pandasgwas` | `remotes::install_github("TianzeLab/pandasGWASr")` |
| 9 `get_*` query functions | ✅ | ✅ |
| Fat / Thin modes | ✅ | ✅ |
| Set operations `+ & - \| ^ ==` | ✅ | ✅ |
| Browser helper functions | ✅ | ✅ |
| Summary statistics submodule | ✅ | ✅ |
| Cache management | ✅ | ✅ |
| Data container | Custom classes with DataFrame attributes | S4 classes with data.frame slots |
| Indexing style | 0-based, e.g. `associations[0]` | 1-based, e.g. `associations[1]` |

## License

- **Source code**: MIT License — see [LICENSE](LICENSE), [Python/LICENSE](Python/LICENSE), and [R/LICENSE](R/LICENSE).
- **GWAS Catalog data**: The NHGRI-EBI GWAS Catalog follows the general [Terms of Use for EMBL-EBI Services](https://www.ebi.ac.uk/about/terms-of-use/). Summary statistics are made available under [CC0](https://creativecommons.org/publicdomain/zero/1.0/) unless [otherwise stated](https://www.ebi.ac.uk/gwas/docs/faq/#faq-H7).
