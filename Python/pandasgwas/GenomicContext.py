from copy import deepcopy
from pandas import DataFrame, set_option

set_option('display.max_columns', None)
set_option('display.width', 1000)
set_option('display.colheader_justify', 'center')
set_option('display.precision', 3)


_GENOMIC_CONTEXT_COLUMNS = [
    'rs_id',
    'genomicContext_id',
    'is_intergenic',
    'is_upstream',
    'is_downstream',
    'distance',
    'source',
    'mapping_method',
    'is_closest_gene'
]

_GENES_COLUMNS = [
    'rs_id',
    'genomicContext_id',
    'gene_name',
    'gene_description',
    'location',
    'cytogenic_region',
    'biotype'
]
_LOCATIONS_COLUMNS = ['rs_id', 'genomicContext_id', 'chromosome_name', 'chromosome_position', 'region.name']
_ENSEMBL_GENE_IDS_COLUMNS = ['rs_id', 'genomicContext_id', 'gene_name', 'ensembl_gene_ids']
_ENTREZ_GENE_IDS_COLUMNS = ['rs_id', 'genomicContext_id', 'gene_name', 'entrez_gene_ids']

_NESTED_FIELDS = [
    'genes',
    'locations',
    'ensembl_gene_ids',
    'entrez_gene_ids'
]


class GenomicContext:
    """
       An object that stores data of type GenomicContext. To understand the significance of each column of the DataFrame. Please visit "GenomicContext" in [Gwas Catalog Documentation](https://www.ebi.ac.uk/gwas/rest/api/v2/docs/reference) for details.

       Attributes:
            raw_data: list. Converted from obtained JSON data.
            genomic_contexts: DataFrame. It only exists if the parameter mode of constructor is Fat.
            genes: DataFrame. It only exists if the parameter mode of constructor is Fat.
            locations: DataFrame. It only exists if the parameter mode of constructor is Fat.
            ensembl_gene_ids: DataFrame. It only exists if the parameter mode of constructor is Fat.
            entrez_gene_ids: DataFrame. It only exists if the parameter mode of constructor is Fat.
            mode: Fat or Thin. Specifies the mode of the returned object.

       ```python
       from pandasgwas.get_genomic_contexts import get_genomic_contexts

       ch = get_genomic_contexts(rs_id='rs3093017')
       ch
       ch.raw_data
       ch.mode
       ch.genomic_contexts
       ch.genes
       ch.locations
       ch.ensembl_gene_ids
       ch.entrez_gene_ids
       ```
       Subset object s by either identifier or position
       ```python
       all_df = get_genomic_contexts(rs_id='rs3093017')
       all_df[0].genomic_contexts
       all_df[0:3].genomic_contexts
       all_df['rs3093017|206217817'].genomic_contexts
       all_df[('rs3093017|206217817','rs3093017|206217818','rs3093017|206217811')].genomic_contexts
       ```
       Objects can be manipulated like sets in the mathematical sense.
       ```python
       one = get_genomic_contexts(rs_id='rs3093017', genomicContext_id='206217817')
       two = get_genomic_contexts(rs_id='rs3093017', genomicContext_id='206217818')
       three = get_genomic_contexts(rs_id='rs3093017', genomicContext_id='206217811')
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
            ('rs_id', ['rs_id', 'rsId']),
            ('genomicContext_id', ['genomicContext_id', 'genomicContextId', 'genomic_context_id']),
            ('is_intergenic', ['is_intergenic', 'isIntergenic']),
            ('is_upstream', ['is_upstream', 'isUpstream']),
            ('is_downstream', ['is_downstream', 'isDownstream']),
            ('mapping_method', ['mapping_method', 'mappingMethod']),
            ('is_closest_gene', ['is_closest_gene', 'isClosestGene']),
            ('genes', ['genes', 'gene']),
            ('locations', ['locations', 'location']),
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
                if 'rs_id' in rec and rec['rs_id'] is not None:
                    rec['rs_id'] = str(rec['rs_id'])
                if 'genomicContext_id' in rec and rec['genomicContext_id'] is not None:
                    rec['genomicContext_id'] = str(rec['genomicContext_id'])
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
            for col in _GENOMIC_CONTEXT_COLUMNS:
                value = normalized.get(col, None)
                standardized[col] = value

            for field in _NESTED_FIELDS:
                standardized[field] = normalized.get(field, None)

            standardized_data.append(standardized)

        self.raw_data = standardized_data

        if len(standardized_data) == 0:
            self.genomic_contexts = DataFrame(columns=_GENOMIC_CONTEXT_COLUMNS)
            self.genes = DataFrame(columns=_GENES_COLUMNS)
            self.locations = DataFrame(columns=_LOCATIONS_COLUMNS)
            self.ensembl_gene_ids = DataFrame(columns=_ENSEMBL_GENE_IDS_COLUMNS)
            self.entrez_gene_ids = DataFrame(columns=_ENTREZ_GENE_IDS_COLUMNS)
            return

        genomic_context_rows = []
        genes_rows = []
        locations_rows = []
        ensembl_gene_ids_rows = []
        entrez_gene_ids_rows = []

        for record in standardized_data:
            rs_id = record.get('rs_id', None)
            genomicContext_id = record.get('genomicContext_id', None)

            genomic_context_rows.append({
                col: record.get(col, None) for col in _GENOMIC_CONTEXT_COLUMNS
            })

            raw_genes = record.get('genes', None)
            if raw_genes is None:
                iterable_genes = []
            elif isinstance(raw_genes, list):
                iterable_genes = raw_genes
            else:
                iterable_genes = [raw_genes]
            for item in iterable_genes:
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
                    genes_rows.append({
                        'rs_id': rs_id,
                        'genomicContext_id': genomicContext_id,
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
                            'rs_id': rs_id,
                            'genomicContext_id': genomicContext_id,
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
                            'rs_id': rs_id,
                            'genomicContext_id': genomicContext_id,
                            'gene_name': gene_name,
                            'entrez_gene_ids': entrez_gene_id
                        })
                else:
                    genes_rows.append({
                        'rs_id': rs_id,
                        'genomicContext_id': genomicContext_id,
                        'gene_name': item,
                        'gene_description': None,
                        'location': None,
                        'cytogenic_region': None,
                        'biotype': None
                    })

            raw_locations = record.get('locations', None)
            if raw_locations is None:
                iterable_locations = []
            elif isinstance(raw_locations, list):
                iterable_locations = raw_locations
            else:
                iterable_locations = [raw_locations]
            for item in iterable_locations:
                if isinstance(item, dict):
                    chromosome_name = item.get('chromosome_name', None)
                    if chromosome_name is None:
                        chromosome_name = item.get('chromosomeName', None)
                    chromosome_position = item.get('chromosome_position', None)
                    if chromosome_position is None:
                        chromosome_position = item.get('chromosomePosition', None)
                    region_name = None
                    raw_region = item.get('region', None)
                    if isinstance(raw_region, dict):
                        region_name = raw_region.get('name', None)
                    elif raw_region is not None:
                        region_name = raw_region
                    locations_rows.append({
                        'rs_id': rs_id,
                        'genomicContext_id': genomicContext_id,
                        'chromosome_name': chromosome_name,
                        'chromosome_position': chromosome_position,
                        'region.name': region_name
                    })
                else:
                    locations_rows.append({
                        'rs_id': rs_id,
                        'genomicContext_id': genomicContext_id,
                        'chromosome_name': None,
                        'chromosome_position': item,
                        'region.name': None
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
                        'rs_id': rs_id,
                        'genomicContext_id': genomicContext_id,
                        'gene_name': gene_name,
                        'ensembl_gene_ids': ensembl_gene_ids
                    })
                else:
                    ensembl_gene_ids_rows.append({
                        'rs_id': rs_id,
                        'genomicContext_id': genomicContext_id,
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
                        'rs_id': rs_id,
                        'genomicContext_id': genomicContext_id,
                        'gene_name': gene_name,
                        'entrez_gene_ids': entrez_gene_ids
                    })
                else:
                    entrez_gene_ids_rows.append({
                        'rs_id': rs_id,
                        'genomicContext_id': genomicContext_id,
                        'gene_name': None,
                        'entrez_gene_ids': item
                    })

        self.genomic_contexts = DataFrame(genomic_context_rows).reindex(columns=_GENOMIC_CONTEXT_COLUMNS) if len(genomic_context_rows) else DataFrame(columns=_GENOMIC_CONTEXT_COLUMNS)
        self.genes = DataFrame(genes_rows).reindex(columns=_GENES_COLUMNS) if len(genes_rows) else DataFrame(columns=_GENES_COLUMNS)
        self.locations = DataFrame(locations_rows).reindex(columns=_LOCATIONS_COLUMNS) if len(locations_rows) else DataFrame(columns=_LOCATIONS_COLUMNS)
        self.ensembl_gene_ids = DataFrame(ensembl_gene_ids_rows).reindex(columns=_ENSEMBL_GENE_IDS_COLUMNS) if len(ensembl_gene_ids_rows) else DataFrame(columns=_ENSEMBL_GENE_IDS_COLUMNS)
        self.entrez_gene_ids = DataFrame(entrez_gene_ids_rows).reindex(columns=_ENTREZ_GENE_IDS_COLUMNS) if len(entrez_gene_ids_rows) else DataFrame(columns=_ENTREZ_GENE_IDS_COLUMNS)

    def __str__(self):
        if self.mode == 'Fat':
            return (
                    "GenomicContext is running in fat mode. It has 5 DataFrames with hierarchical dependencies.\n"
                    "genomic_contexts: %d rows\n"
                    "|\n -genes: %d rows\n"
                    "|\n -locations: %d rows\n"
                    "|\n -ensembl_gene_ids: %d rows\n"
                    "|\n -entrez_gene_ids: %d rows"
                    % (
                        len(self.genomic_contexts),
                        len(self.genes),
                        len(self.locations),
                        len(self.ensembl_gene_ids),
                        len(self.entrez_gene_ids),
                    )
            )
        if self.mode == 'Thin':
            return (
                    'GenomicContext is running in thin mode. It has 1 list that contains the raw data.\n'
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
            raw_data_dict['{}|{}'.format(j['rs_id'], j['genomicContext_id'])] = j

        sub_set = []
        for i in arr:
            if isinstance(i, str):
                sub_set.append(raw_data_dict[i])
            elif isinstance(i, int):
                sub_set.append(raw_data[i])
            else:
                raise TypeError('Invalid item type: {}'.format(type(i)))

        return GenomicContext(sub_set, self.mode)

    def __add__(self, other):
        if self.mode == other.mode:
            return GenomicContext(self.raw_data + other.raw_data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __sub__(self, other):
        if self.mode == other.mode:
            self_key_set = set()
            self_dict = {}
            other_key_set = set()

            for i in self.raw_data:
                key = '{}/{}'.format(i['rs_id'], i['genomicContext_id'])
                self_key_set.add(key)
                self_dict[key] = i

            for j in other.raw_data:
                other_key_set.add('{}/{}'.format(j['rs_id'], j['genomicContext_id']))

            sub_key = self_key_set - other_key_set
            data = []
            for k in sub_key:
                data.append(self_dict[k])

            return GenomicContext(data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __and__(self, other):
        if self.mode == other.mode:
            self_key_set = set()
            self_dict = {}
            other_key_set = set()

            for i in self.raw_data:
                key = '{}/{}'.format(i['rs_id'], i['genomicContext_id'])
                self_key_set.add(key)
                self_dict[key] = i

            for j in other.raw_data:
                other_key_set.add('{}/{}'.format(j['rs_id'], j['genomicContext_id']))

            sub_key = self_key_set & other_key_set
            data = []
            for k in sub_key:
                data.append(self_dict[k])

            return GenomicContext(data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __or__(self, other):
        if self.mode == other.mode:
            and_dict = {}

            for i in self.raw_data:
                and_dict['{}/{}'.format(i['rs_id'], i['genomicContext_id'])] = i

            for j in other.raw_data:
                and_dict['{}/{}'.format(j['rs_id'], j['genomicContext_id'])] = j

            data = list(and_dict.values())
            return GenomicContext(data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __xor__(self, other):
        if self.mode == other.mode:
            self_key_set = set()
            and_dict = {}
            other_key_set = set()

            for i in self.raw_data:
                key = '{}/{}'.format(i['rs_id'], i['genomicContext_id'])
                self_key_set.add(key)
                and_dict[key] = i

            for j in other.raw_data:
                key = '{}/{}'.format(j['rs_id'], j['genomicContext_id'])
                other_key_set.add(key)
                and_dict[key] = j

            sub_key = self_key_set ^ other_key_set
            data = []
            for k in sub_key:
                data.append(and_dict[k])

            return GenomicContext(data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __eq__(self, other):
        if self is None or other is None:
            return self is None and other is None
        return self.raw_data == other.raw_data and self.mode == other.mode

    def __len__(self):
        return len(self.raw_data)
