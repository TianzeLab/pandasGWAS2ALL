from copy import deepcopy
from pandas import DataFrame, set_option

set_option('display.max_columns', None)
set_option('display.width', 1000)
set_option('display.colheader_justify', 'center')
set_option('display.precision', 3)


_SNP_COLUMNS = [
    'rs_id',
    'merged',
    'functional_class',
    'current_snp',
    'last_update_date',
    'maf',
    'minor_allele',
    'alleles',
    'most_severe_consequence'
]

_LOCATIONS_COLUMNS = ['rs_id', 'chromosome_name', 'chromosome_position', 'region.name']
_MAPPED_GENES_COLUMNS = ['rs_id', 'mapped_genes']

_NESTED_FIELDS = [
    'locations',
    'mapped_genes'
]


class Snp:
    """
       An object that stores data of type Snp. To understand the significance of each column of the DataFrame. Please visit "Snp" in [Gwas Catalog Documentation](https://www.ebi.ac.uk/gwas/rest/api/v2/docs/reference) for details.

       Attributes:
            raw_data: list. Converted from obtained JSON data.
            snps: DataFrame. It only exists if the parameter mode of constructor is Fat.
            locations: DataFrame. It only exists if the parameter mode of constructor is Fat.
            mapped_genes: DataFrame. It only exists if the parameter mode of constructor is Fat.
            mode: Fat or Thin. Specifies the mode of the returned object.

       ```python
       from pandasgwas.get_snps import get_snps

       ch = get_snps(rs_id='rs3093017')
       ch
       ch.raw_data
       ch.mode
       ch.snps
       ch.locations
       ch.mapped_genes
       ```
       Subset object s by either identifier or position
       ```python
       all_df = get_snps()
       all_df[0].snps
       all_df[0:3].snps
       all_df['rs3093017'].snps
       all_df[('rs3093017','rs7671745','rs61953412')].snps
       ```
       Objects can be manipulated like sets in the mathematical sense.
       ```python
       one = get_snps(rs_id='rs3093017')
       two = get_snps(rs_id='rs7671745')
       three = get_snps(rs_id='rs61953412')
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
                if 'rs_id' in rec and rec['rs_id'] is not None:
                    rec['rs_id'] = str(rec['rs_id'])
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
            for col in _SNP_COLUMNS:
                value = normalized.get(col, None)
                standardized[col] = value

            for field in _NESTED_FIELDS:
                standardized[field] = normalized.get(field, None)

            standardized_data.append(standardized)

        self.raw_data = standardized_data

        if len(standardized_data) == 0:
            self.snps = DataFrame(columns=_SNP_COLUMNS)
            self.locations = DataFrame(columns=_LOCATIONS_COLUMNS)
            self.mapped_genes = DataFrame(columns=_MAPPED_GENES_COLUMNS)
            return

        snp_rows = []
        locations_rows = []
        mapped_genes_rows = []

        for record in standardized_data:
            rs_id = record.get('rs_id', None)

            # main table: fixed fields, missing fields filled with None
            snp_rows.append({
                col: record.get(col, None) for col in _SNP_COLUMNS
            })

            # locations
            raw_locations = record.get('locations', None)
            if raw_locations is None:
                iterable_locations = []
            elif isinstance(raw_locations, list):
                iterable_locations = raw_locations
            else:
                iterable_locations = [raw_locations]
            for item in iterable_locations:
                if isinstance(item, dict):
                    raw_region = item.get('region', None)
                    if raw_region is None:
                        region_name = None
                    elif isinstance(raw_region, dict):
                        region_name = raw_region.get('name', None)
                    else:
                        region_name = raw_region

                    locations_rows.append({
                        'rs_id': rs_id,
                        'chromosome_name': item.get('chromosome_name', None),
                        'chromosome_position': item.get('chromosome_position', None),
                        'region.name': region_name
                    })
                else:
                    locations_rows.append({
                        'rs_id': rs_id,
                        'chromosome_name': item,
                        'chromosome_position': None,
                        'region.name': None
                    })

            # mapped_genes
            raw_mapped_genes = record.get('mapped_genes', None)
            if raw_mapped_genes is None:
                iterable_mapped_genes = []
            elif isinstance(raw_mapped_genes, list):
                iterable_mapped_genes = raw_mapped_genes
            else:
                iterable_mapped_genes = [raw_mapped_genes]
            for item in iterable_mapped_genes:
                mapped_genes_rows.append({
                    'rs_id': rs_id,
                    'mapped_genes': item
                })

        self.snps = DataFrame(snp_rows).reindex(columns=_SNP_COLUMNS) if len(snp_rows) else DataFrame(columns=_SNP_COLUMNS)
        self.locations = DataFrame(locations_rows).reindex(columns=_LOCATIONS_COLUMNS) if len(locations_rows) else DataFrame(columns=_LOCATIONS_COLUMNS)
        self.mapped_genes = DataFrame(mapped_genes_rows).reindex(columns=_MAPPED_GENES_COLUMNS) if len(mapped_genes_rows) else DataFrame(columns=_MAPPED_GENES_COLUMNS)

    def __str__(self):
        if self.mode == 'Fat':
            return (
                    "Snp is running in fat mode. It has 3 DataFrames with hierarchical dependencies.\n"
                    "snps: %d rows\n"
                    "|\n -locations: %d rows\n"
                    "|\n -mapped_genes: %d rows"
                    % (
                        len(self.snps),
                        len(self.locations),
                        len(self.mapped_genes),
                    )
            )
        if self.mode == 'Thin':
            return (
                    'Snp is running in thin mode. It has 1 list that contains the raw data.\n'
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
            raw_data_dict[j['rs_id']] = j

        sub_set = []
        for i in arr:
            if isinstance(i, str):
                sub_set.append(raw_data_dict[i])
            elif isinstance(i, int):
                sub_set.append(raw_data[i])
            else:
                raise TypeError('Invalid item type: {}'.format(type(i)))

        return Snp(sub_set, self.mode)

    def __add__(self, other):
        if self.mode == other.mode:
            return Snp(self.raw_data + other.raw_data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __sub__(self, other):
        if self.mode == other.mode:
            self_key_set = set()
            self_dict = {}
            other_key_set = set()

            for i in self.raw_data:
                self_key_set.add(i['rs_id'])
                self_dict[i['rs_id']] = i

            for j in other.raw_data:
                other_key_set.add(j['rs_id'])

            sub_key = self_key_set - other_key_set
            data = []
            for k in sub_key:
                data.append(self_dict[k])

            return Snp(data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __and__(self, other):
        if self.mode == other.mode:
            self_key_set = set()
            self_dict = {}
            other_key_set = set()

            for i in self.raw_data:
                self_key_set.add(i['rs_id'])
                self_dict[i['rs_id']] = i

            for j in other.raw_data:
                other_key_set.add(j['rs_id'])

            sub_key = self_key_set & other_key_set
            data = []
            for k in sub_key:
                data.append(self_dict[k])

            return Snp(data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __or__(self, other):
        if self.mode == other.mode:
            and_dict = {}

            for i in self.raw_data:
                and_dict[i['rs_id']] = i

            for j in other.raw_data:
                and_dict[j['rs_id']] = j

            data = list(and_dict.values())
            return Snp(data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __xor__(self, other):
        if self.mode == other.mode:
            self_key_set = set()
            and_dict = {}
            other_key_set = set()

            for i in self.raw_data:
                self_key_set.add(i['rs_id'])
                and_dict[i['rs_id']] = i

            for j in other.raw_data:
                other_key_set.add(j['rs_id'])
                and_dict[j['rs_id']] = j

            sub_key = self_key_set ^ other_key_set
            data = []
            for k in sub_key:
                data.append(and_dict[k])

            return Snp(data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __eq__(self, other):
        if self is None or other is None:
            return self is None and other is None
        return self.raw_data == other.raw_data and self.mode == other.mode

    def __len__(self):
        return len(self.raw_data)
