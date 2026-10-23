from copy import deepcopy
from pandas import DataFrame, set_option

set_option('display.max_columns', None)
set_option('display.width', 1000)
set_option('display.colheader_justify', 'center')
set_option('display.precision', 3)


_LOCUS_COLUMNS = [
    'association_id',
    'locus_id',
    'haplotype_snp_count',
    'description'
]

_STRONGEST_RISK_ALLELES_COLUMNS = [
    'association_id',
    'locus_id',
    'risk_allele_name',
    'risk_frequency',
    'genome_wide',
    'limited_list',
    'source',
    'mapping_method',
    'is_closest_gene'
]
_AUTHOR_REPORTED_GENES_COLUMNS = [
    'association_id',
    'locus_id',
    'gene_name',
    'gene_description',
    'location',
    'cytogenic_region',
    'biotype'
]
_ENSEMBL_GENE_IDS_COLUMNS = ['association_id', 'locus_id', 'gene_name', 'ensembl_gene_ids']
_ENTREZ_GENE_IDS_COLUMNS = ['association_id', 'locus_id', 'gene_name', 'entrez_gene_ids']

_NESTED_FIELDS = [
    'strongest_risk_alleles',
    'author_reported_genes',
    'ensembl_gene_ids',
    'entrez_gene_ids'
]


class Locus:
    """
       An object that stores data of type Locus. To understand the significance of each column of the DataFrame. Please visit "Locus" in [Gwas Catalog Documentation](https://www.ebi.ac.uk/gwas/rest/api/v2/docs/reference) for details.

       Attributes:
            raw_data: list. Converted from obtained JSON data.
            loci: DataFrame. It only exists if the parameter mode of constructor is Fat.
            strongest_risk_alleles: DataFrame. It only exists if the parameter mode of constructor is Fat.
            author_reported_genes: DataFrame. It only exists if the parameter mode of constructor is Fat.
            ensembl_gene_ids: DataFrame. It only exists if the parameter mode of constructor is Fat.
            entrez_gene_ids: DataFrame. It only exists if the parameter mode of constructor is Fat.
            mode: Fat or Thin. Specifies the mode of the returned object.

       ```python
       from pandasgwas.get_loci import get_loci

       ch = get_loci(association_id='169106883')
       ch
       ch.raw_data
       ch.mode
       ch.loci
       ch.strongest_risk_alleles
       ch.author_reported_genes
       ch.ensembl_gene_ids
       ch.entrez_gene_ids
       ```
       Subset object s by either identifier or position
       ```python
       all_df = get_loci(association_id='169106883')
       all_df[0].loci
       all_df[0:1].loci
       all_df['169106883|169106882'].loci
       all_df[('169106883|169106882')].loci
       ```
       Objects can be manipulated like sets in the mathematical sense.
       ```python
       one = get_loci(association_id='169106883', locus_id='169106882')
       two = get_loci(association_id='218958492', locus_id='218958491')
       three = get_loci(association_id='218958441', locus_id='218958440')
       one_and_two = one+two
       two_and_three = two+three
       only_one = one_and_two - two_and_three
       only_two = one_and_two & two_and_three
       one_and_two_and_three = one_and_two | two_and_three
       one_and_three = one_and_two ^ two_and_three
       ```
       """
    def _normalize_record(self, record):
        normalized = deepcopy(record)
        aliases = [
            ('association_id', ['association_id', 'associationId']),
            ('locus_id', ['locus_id', 'locusId']),
            ('haplotype_snp_count', ['haplotype_snp_count', 'haplotypeSnpCount']),
            ('strongest_risk_alleles', ['strongest_risk_alleles', 'strongestRiskAlleles']),
            ('author_reported_genes', ['author_reported_genes', 'authorReportedGenes', 'genes', 'gene']),
            ('ensembl_gene_ids', ['ensembl_gene_ids', 'ensemblGeneIds']),
            ('entrez_gene_ids', ['entrez_gene_ids', 'entrezGeneIds'])
        ]
        for key, source_keys in aliases:
            value = None
            for source_key in source_keys:
                if source_key in record and record[source_key] is not None:
                    value = record[source_key]
                    break
            if value is not None:
                normalized[key] = value
        return normalized

    def __init__(self, data: list = None, mode: str = "Fat"):
        if mode not in ['Thin', 'Fat']:
            raise Exception("Mode must be Fat or Thin")
        if data is None:
            data = []
        else:
            normalized_data = []
            for item in data:
                if not isinstance(item, dict):
                    continue
                rec = deepcopy(item)
                rec = self._normalize_record(rec)
                if 'association_id' in rec and rec['association_id'] is not None:
                    rec['association_id'] = str(rec['association_id'])
                if 'locus_id' in rec and rec['locus_id'] is not None:
                    rec['locus_id'] = str(rec['locus_id'])
                normalized_data.append(rec)
            data = normalized_data

        self.mode = mode

        if mode == "Thin":
            raw_data = []
            for item in data:
                if not isinstance(item, dict):
                    continue
                rec = deepcopy(item)
                raw_data.append(rec)
            self.raw_data = raw_data
            return

        standardized_data = []
        for item in data:
            if not isinstance(item, dict):
                continue

            rec = deepcopy(item)
            rec = self._normalize_record(rec)

            normalized = {}
            for key, value in rec.items():
                normalized[key] = value

            standardized = {}
            for col in _LOCUS_COLUMNS:
                value = normalized.get(col, None)
                standardized[col] = value

            for field in _NESTED_FIELDS:
                standardized[field] = normalized.get(field, None)

            standardized_data.append(standardized)

        self.raw_data = standardized_data

        if len(standardized_data) == 0:
            self.loci = DataFrame(columns=_LOCUS_COLUMNS)
            self.strongest_risk_alleles = DataFrame(columns=_STRONGEST_RISK_ALLELES_COLUMNS)
            self.author_reported_genes = DataFrame(columns=_AUTHOR_REPORTED_GENES_COLUMNS)
            self.ensembl_gene_ids = DataFrame(columns=_ENSEMBL_GENE_IDS_COLUMNS)
            self.entrez_gene_ids = DataFrame(columns=_ENTREZ_GENE_IDS_COLUMNS)
            return

        locus_rows = []
        strongest_risk_alleles_rows = []
        author_reported_genes_rows = []
        ensembl_gene_ids_rows = []
        entrez_gene_ids_rows = []

        for record in standardized_data:
            association_id = record.get('association_id', None)
            locus_id = record.get('locus_id', None)

            locus_rows.append({
                col: record.get(col, None) for col in _LOCUS_COLUMNS
            })

            raw_strongest_risk_alleles = record.get('strongest_risk_alleles', None)
            if raw_strongest_risk_alleles is None:
                iterable_strongest_risk_alleles = []
            elif isinstance(raw_strongest_risk_alleles, list):
                iterable_strongest_risk_alleles = raw_strongest_risk_alleles
            else:
                iterable_strongest_risk_alleles = [raw_strongest_risk_alleles]
            for item in iterable_strongest_risk_alleles:
                if isinstance(item, dict):
                    risk_allele_name = item.get('risk_allele_name', None)
                    if risk_allele_name is None:
                        risk_allele_name = item.get('riskAlleleName', None)
                    risk_frequency = item.get('risk_frequency', None)
                    if risk_frequency is None:
                        risk_frequency = item.get('riskFrequency', None)
                    genome_wide = item.get('genome_wide', None)
                    if genome_wide is None:
                        genome_wide = item.get('genomeWide', None)
                    limited_list = item.get('limited_list', None)
                    if limited_list is None:
                        limited_list = item.get('limitedList', None)
                    mapping_method = item.get('mapping_method', None)
                    if mapping_method is None:
                        mapping_method = item.get('mappingMethod', None)
                    is_closest_gene = item.get('is_closest_gene', None)
                    if is_closest_gene is None:
                        is_closest_gene = item.get('isClosestGene', None)
                    strongest_risk_alleles_rows.append({
                        'association_id': association_id,
                        'locus_id': locus_id,
                        'risk_allele_name': risk_allele_name,
                        'risk_frequency': risk_frequency,
                        'genome_wide': genome_wide,
                        'limited_list': limited_list,
                        'source': item.get('source', None),
                        'mapping_method': mapping_method,
                        'is_closest_gene': is_closest_gene
                    })
                else:
                    strongest_risk_alleles_rows.append({
                        'association_id': association_id,
                        'locus_id': locus_id,
                        'risk_allele_name': item,
                        'risk_frequency': None,
                        'genome_wide': None,
                        'limited_list': None,
                        'source': None,
                        'mapping_method': None,
                        'is_closest_gene': None
                    })

            raw_author_reported_genes = record.get('author_reported_genes', None)
            if raw_author_reported_genes is None:
                iterable_author_reported_genes = []
            elif isinstance(raw_author_reported_genes, list):
                iterable_author_reported_genes = raw_author_reported_genes
            else:
                iterable_author_reported_genes = [raw_author_reported_genes]
            for item in iterable_author_reported_genes:
                if isinstance(item, dict):
                    gene_name = item.get('gene_name', None)
                    if gene_name is None:
                        gene_name = item.get('geneName', None)
                    gene_description = item.get('gene_description', None)
                    if gene_description is None:
                        gene_description = item.get('geneDescription', None)
                    cytogenic_region = item.get('cytogenic_region', None)
                    if cytogenic_region is None:
                        cytogenic_region = item.get('cytogenetic_region', None)
                    if cytogenic_region is None:
                        cytogenic_region = item.get('cytogenicRegion', None)
                    if cytogenic_region is None:
                        cytogenic_region = item.get('cytogeneticRegion', None)
                    author_reported_genes_rows.append({
                        'association_id': association_id,
                        'locus_id': locus_id,
                        'gene_name': gene_name,
                        'gene_description': gene_description,
                        'location': item.get('location', None),
                        'cytogenic_region': cytogenic_region,
                        'biotype': item.get('biotype', None)
                    })

                    raw_ensembl_gene_ids = item.get('ensembl_gene_ids', None)
                    if raw_ensembl_gene_ids is None:
                        raw_ensembl_gene_ids = item.get('ensemblGeneIds', None)
                    if raw_ensembl_gene_ids is None:
                        iterable_ensembl_gene_ids = []
                    elif isinstance(raw_ensembl_gene_ids, list):
                        iterable_ensembl_gene_ids = raw_ensembl_gene_ids
                    else:
                        iterable_ensembl_gene_ids = [raw_ensembl_gene_ids]
                    for ensembl_gene_id in iterable_ensembl_gene_ids:
                        ensembl_gene_ids_rows.append({
                            'association_id': association_id,
                            'locus_id': locus_id,
                            'gene_name': gene_name,
                            'ensembl_gene_ids': ensembl_gene_id
                        })

                    raw_entrez_gene_ids = item.get('entrez_gene_ids', None)
                    if raw_entrez_gene_ids is None:
                        raw_entrez_gene_ids = item.get('entrezGeneIds', None)
                    if raw_entrez_gene_ids is None:
                        iterable_entrez_gene_ids = []
                    elif isinstance(raw_entrez_gene_ids, list):
                        iterable_entrez_gene_ids = raw_entrez_gene_ids
                    else:
                        iterable_entrez_gene_ids = [raw_entrez_gene_ids]
                    for entrez_gene_id in iterable_entrez_gene_ids:
                        entrez_gene_ids_rows.append({
                            'association_id': association_id,
                            'locus_id': locus_id,
                            'gene_name': gene_name,
                            'entrez_gene_ids': entrez_gene_id
                        })
                else:
                    author_reported_genes_rows.append({
                        'association_id': association_id,
                        'locus_id': locus_id,
                        'gene_name': item,
                        'gene_description': None,
                        'location': None,
                        'cytogenic_region': None,
                        'biotype': None
                    })

            raw_ensembl_gene_ids = record.get('ensembl_gene_ids', None)
            if raw_ensembl_gene_ids is None:
                iterable_ensembl_gene_ids = []
            elif isinstance(raw_ensembl_gene_ids, list):
                iterable_ensembl_gene_ids = raw_ensembl_gene_ids
            else:
                iterable_ensembl_gene_ids = [raw_ensembl_gene_ids]
            for item in iterable_ensembl_gene_ids:
                if isinstance(item, dict):
                    gene_name = item.get('gene_name', None)
                    if gene_name is None:
                        gene_name = item.get('geneName', None)
                    ensembl_gene_ids = item.get('ensembl_gene_ids', None)
                    if ensembl_gene_ids is None:
                        ensembl_gene_ids = item.get('ensemblGeneIds', None)
                    ensembl_gene_ids_rows.append({
                        'association_id': association_id,
                        'locus_id': locus_id,
                        'gene_name': gene_name,
                        'ensembl_gene_ids': ensembl_gene_ids
                    })
                else:
                    ensembl_gene_ids_rows.append({
                        'association_id': association_id,
                        'locus_id': locus_id,
                        'gene_name': None,
                        'ensembl_gene_ids': item
                    })

            raw_entrez_gene_ids = record.get('entrez_gene_ids', None)
            if raw_entrez_gene_ids is None:
                iterable_entrez_gene_ids = []
            elif isinstance(raw_entrez_gene_ids, list):
                iterable_entrez_gene_ids = raw_entrez_gene_ids
            else:
                iterable_entrez_gene_ids = [raw_entrez_gene_ids]
            for item in iterable_entrez_gene_ids:
                if isinstance(item, dict):
                    gene_name = item.get('gene_name', None)
                    if gene_name is None:
                        gene_name = item.get('geneName', None)
                    entrez_gene_ids = item.get('entrez_gene_ids', None)
                    if entrez_gene_ids is None:
                        entrez_gene_ids = item.get('entrezGeneIds', None)
                    entrez_gene_ids_rows.append({
                        'association_id': association_id,
                        'locus_id': locus_id,
                        'gene_name': gene_name,
                        'entrez_gene_ids': entrez_gene_ids
                    })
                else:
                    entrez_gene_ids_rows.append({
                        'association_id': association_id,
                        'locus_id': locus_id,
                        'gene_name': None,
                        'entrez_gene_ids': item
                    })

        self.loci = DataFrame(locus_rows).reindex(columns=_LOCUS_COLUMNS) if len(locus_rows) else DataFrame(columns=_LOCUS_COLUMNS)
        self.strongest_risk_alleles = DataFrame(strongest_risk_alleles_rows).reindex(columns=_STRONGEST_RISK_ALLELES_COLUMNS) if len(strongest_risk_alleles_rows) else DataFrame(columns=_STRONGEST_RISK_ALLELES_COLUMNS)
        self.author_reported_genes = DataFrame(author_reported_genes_rows).reindex(columns=_AUTHOR_REPORTED_GENES_COLUMNS) if len(author_reported_genes_rows) else DataFrame(columns=_AUTHOR_REPORTED_GENES_COLUMNS)
        self.ensembl_gene_ids = DataFrame(ensembl_gene_ids_rows).reindex(columns=_ENSEMBL_GENE_IDS_COLUMNS) if len(ensembl_gene_ids_rows) else DataFrame(columns=_ENSEMBL_GENE_IDS_COLUMNS)
        self.entrez_gene_ids = DataFrame(entrez_gene_ids_rows).reindex(columns=_ENTREZ_GENE_IDS_COLUMNS) if len(entrez_gene_ids_rows) else DataFrame(columns=_ENTREZ_GENE_IDS_COLUMNS)

    def __str__(self):
        if self.mode == 'Fat':
            return (
                    "Locus is running in fat mode. It has 5 DataFrames with hierarchical dependencies.\n"
                    "loci: %d rows\n"
                    "|\n -strongest_risk_alleles: %d rows\n"
                    "|\n -author_reported_genes: %d rows\n"
                    "|\n -ensembl_gene_ids: %d rows\n"
                    "|\n -entrez_gene_ids: %d rows"
                    % (
                        len(self.loci),
                        len(self.strongest_risk_alleles),
                        len(self.author_reported_genes),
                        len(self.ensembl_gene_ids),
                        len(self.entrez_gene_ids),
                    )
            )
        if self.mode == 'Thin':
            return (
                    'Locus is running in thin mode. It has 1 list that contains the raw data.\n'
                    'raw_data: a list of size %d.'
                    % len(self.raw_data)
            )

    def __repr__(self):
        return self.__str__()

    def __getitem__(self, item):
        if isinstance(item, (str, int)):
            arr = [item]
        elif isinstance(item, (list, tuple, range)):
            arr = item
        elif isinstance(item, slice):
            start, stop, step = item.indices(len(self))
            arr = list(range(start, stop, step))
        else:
            raise TypeError('Invalid argument type: {}'.format(type(item)))

        raw_data = self.raw_data
        raw_data_dict = {}
        for j in raw_data:
            raw_data_dict['{}|{}'.format(j['association_id'], j['locus_id'])] = j

        sub_set = []
        for i in arr:
            if isinstance(i, str):
                sub_set.append(raw_data_dict[i])
            elif isinstance(i, int):
                sub_set.append(raw_data[i])
            else:
                raise TypeError('Invalid item type: {}'.format(type(i)))

        return Locus(sub_set, self.mode)

    def __add__(self, other):
        if self.mode == other.mode:
            return Locus(self.raw_data + other.raw_data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __sub__(self, other):
        if self.mode == other.mode:
            self_key_set = set()
            self_dict = {}
            other_key_set = set()

            for i in self.raw_data:
                key = '{}/{}'.format(i['association_id'], i['locus_id'])
                self_key_set.add(key)
                self_dict[key] = i

            for j in other.raw_data:
                other_key_set.add('{}/{}'.format(j['association_id'], j['locus_id']))

            sub_key = self_key_set - other_key_set
            data = []
            for k in sub_key:
                data.append(self_dict[k])

            return Locus(data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __and__(self, other):
        if self.mode == other.mode:
            self_key_set = set()
            self_dict = {}
            other_key_set = set()

            for i in self.raw_data:
                key = '{}/{}'.format(i['association_id'], i['locus_id'])
                self_key_set.add(key)
                self_dict[key] = i

            for j in other.raw_data:
                other_key_set.add('{}/{}'.format(j['association_id'], j['locus_id']))

            sub_key = self_key_set & other_key_set
            data = []
            for k in sub_key:
                data.append(self_dict[k])

            return Locus(data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __or__(self, other):
        if self.mode == other.mode:
            and_dict = {}

            for i in self.raw_data:
                and_dict['{}/{}'.format(i['association_id'], i['locus_id'])] = i

            for j in other.raw_data:
                and_dict['{}/{}'.format(j['association_id'], j['locus_id'])] = j

            data = list(and_dict.values())
            return Locus(data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __xor__(self, other):
        if self.mode == other.mode:
            self_key_set = set()
            and_dict = {}
            other_key_set = set()

            for i in self.raw_data:
                key = '{}/{}'.format(i['association_id'], i['locus_id'])
                self_key_set.add(key)
                and_dict[key] = i

            for j in other.raw_data:
                key = '{}/{}'.format(j['association_id'], j['locus_id'])
                other_key_set.add(key)
                and_dict[key] = j

            sub_key = self_key_set ^ other_key_set
            data = []
            for k in sub_key:
                data.append(and_dict[k])

            return Locus(data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __eq__(self, other):
        if self is None or other is None:
            return self is None and other is None
        return self.raw_data == other.raw_data and self.mode == other.mode

    def __len__(self):
        return len(self.raw_data)
