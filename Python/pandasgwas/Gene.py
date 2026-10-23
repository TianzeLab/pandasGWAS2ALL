from copy import deepcopy
from pandas import DataFrame, set_option

set_option('display.max_columns', None)
set_option('display.width', 1000)
set_option('display.colheader_justify', 'center')
set_option('display.precision', 3)


_GENE_COLUMNS = [
    'gene_name',
    'gene_description',
    'location',
    'cytogenic_region',
    'biotype'
]

_ENSEMBL_GENE_IDS_COLUMNS = ['gene_name', 'ensembl_gene_ids']
_ENTREZ_GENE_IDS_COLUMNS = ['gene_name', 'entrez_gene_ids']

_NESTED_FIELDS = [
    'ensembl_gene_ids',
    'entrez_gene_ids'
]


class Gene:
    """
       An object that stores data of type Gene. To understand the significance of each column of the DataFrame. Please visit "Gene" in [Gwas Catalog Documentation](https://www.ebi.ac.uk/gwas/rest/api/v2/docs/reference) for details.

       Attributes:
            raw_data: list. Converted from obtained JSON data.
            genes: DataFrame. It only exists if the parameter mode of constructor is Fat.
            ensembl_gene_ids: DataFrame. It only exists if the parameter mode of constructor is Fat.
            entrez_gene_ids: DataFrame. It only exists if the parameter mode of constructor is Fat.
            mode: Fat or Thin. Specifies the mode of the returned object.

       ```python
       from pandasgwas.get_genes import get_genes

       ch = get_genes(gene_name='HBS1L')
       ch
       ch.raw_data
       ch.mode
       ch.genes
       ch.ensembl_gene_ids
       ch.entrez_gene_ids
       ```
       Subset object s by either identifier or position
       ```python
       all_df = get_genes()
       all_df[0].genes
       all_df[0:3].genes
       all_df['SAMD5'].genes
       all_df[('SAMD5','SASH1','PSORS1C2')].genes
       ```
       Objects can be manipulated like sets in the mathematical sense.
       ```python
       one = get_genes(gene_name='HBS1L')
       two = get_genes(gene_name='SAMD5')
       three = get_genes(gene_name='SASH1')
       one_and_two = one+two
       two_and_three = two+three
       only_one = one_and_two - two_and_three
       only_two = one_and_two & two_and_three
       one_and_two_and_three = one_and_two | two_and_three
       one_and_three = one_and_two ^ two_and_three
       ```
       """
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
                if 'gene_name' in rec and rec['gene_name'] is not None:
                    rec['gene_name'] = str(rec['gene_name'])
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

            normalized = {}
            for key, value in rec.items():
                normalized[key] = value

            standardized = {}
            for col in _GENE_COLUMNS:
                value = normalized.get(col, None)
                standardized[col] = value

            for field in _NESTED_FIELDS:
                standardized[field] = normalized.get(field, None)

            standardized_data.append(standardized)

        self.raw_data = standardized_data

        if len(standardized_data) == 0:
            self.genes = DataFrame(columns=_GENE_COLUMNS)
            self.ensembl_gene_ids = DataFrame(columns=_ENSEMBL_GENE_IDS_COLUMNS)
            self.entrez_gene_ids = DataFrame(columns=_ENTREZ_GENE_IDS_COLUMNS)
            return

        gene_rows = []
        ensembl_gene_ids_rows = []
        entrez_gene_ids_rows = []

        for record in standardized_data:
            gene_name = record.get('gene_name', None)

            gene_rows.append({
                col: record.get(col, None) for col in _GENE_COLUMNS
            })

            raw_ensembl_gene_ids = record.get('ensembl_gene_ids', None)
            if raw_ensembl_gene_ids is None:
                iterable_ensembl_gene_ids = []
            elif isinstance(raw_ensembl_gene_ids, list):
                iterable_ensembl_gene_ids = raw_ensembl_gene_ids
            else:
                iterable_ensembl_gene_ids = [raw_ensembl_gene_ids]
            for item in iterable_ensembl_gene_ids:
                ensembl_gene_ids_rows.append({
                    'gene_name': gene_name,
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
                entrez_gene_ids_rows.append({
                    'gene_name': gene_name,
                    'entrez_gene_ids': item
                })

        self.genes = DataFrame(gene_rows).reindex(columns=_GENE_COLUMNS) if len(gene_rows) else DataFrame(columns=_GENE_COLUMNS)
        self.ensembl_gene_ids = DataFrame(ensembl_gene_ids_rows).reindex(columns=_ENSEMBL_GENE_IDS_COLUMNS) if len(ensembl_gene_ids_rows) else DataFrame(columns=_ENSEMBL_GENE_IDS_COLUMNS)
        self.entrez_gene_ids = DataFrame(entrez_gene_ids_rows).reindex(columns=_ENTREZ_GENE_IDS_COLUMNS) if len(entrez_gene_ids_rows) else DataFrame(columns=_ENTREZ_GENE_IDS_COLUMNS)

    def __str__(self):
        if self.mode == 'Fat':
            return (
                    "Gene is running in fat mode. It has 3 DataFrames with hierarchical dependencies.\n"
                    "genes: %d rows\n"
                    "|\n -ensembl_gene_ids: %d rows\n"
                    "|\n -entrez_gene_ids: %d rows"
                    % (
                        len(self.genes),
                        len(self.ensembl_gene_ids),
                        len(self.entrez_gene_ids),
                    )
            )
        if self.mode == 'Thin':
            return (
                    'Gene is running in thin mode. It has 1 list that contains the raw data.\n'
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
            raw_data_dict[j['gene_name']] = j

        sub_set = []
        for i in arr:
            if isinstance(i, str):
                sub_set.append(raw_data_dict[i])
            elif isinstance(i, int):
                sub_set.append(raw_data[i])
            else:
                raise TypeError('Invalid item type: {}'.format(type(i)))

        return Gene(sub_set, self.mode)

    def __add__(self, other):
        if self.mode == other.mode:
            return Gene(self.raw_data + other.raw_data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __sub__(self, other):
        if self.mode == other.mode:
            self_key_set = set()
            self_dict = {}
            other_key_set = set()

            for i in self.raw_data:
                self_key_set.add(i['gene_name'])
                self_dict[i['gene_name']] = i

            for j in other.raw_data:
                other_key_set.add(j['gene_name'])

            sub_key = self_key_set - other_key_set
            data = []
            for k in sub_key:
                data.append(self_dict[k])

            return Gene(data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __and__(self, other):
        if self.mode == other.mode:
            self_key_set = set()
            self_dict = {}
            other_key_set = set()

            for i in self.raw_data:
                self_key_set.add(i['gene_name'])
                self_dict[i['gene_name']] = i

            for j in other.raw_data:
                other_key_set.add(j['gene_name'])

            sub_key = self_key_set & other_key_set
            data = []
            for k in sub_key:
                data.append(self_dict[k])

            return Gene(data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __or__(self, other):
        if self.mode == other.mode:
            and_dict = {}

            for i in self.raw_data:
                and_dict[i['gene_name']] = i

            for j in other.raw_data:
                and_dict[j['gene_name']] = j

            data = list(and_dict.values())
            return Gene(data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __xor__(self, other):
        if self.mode == other.mode:
            self_key_set = set()
            and_dict = {}
            other_key_set = set()

            for i in self.raw_data:
                self_key_set.add(i['gene_name'])
                and_dict[i['gene_name']] = i

            for j in other.raw_data:
                other_key_set.add(j['gene_name'])
                and_dict[j['gene_name']] = j

            sub_key = self_key_set ^ other_key_set
            data = []
            for k in sub_key:
                data.append(and_dict[k])

            return Gene(data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __eq__(self, other):
        if self is None or other is None:
            return self is None and other is None
        return self.raw_data == other.raw_data and self.mode == other.mode

    def __len__(self):
        return len(self.raw_data)
