from copy import deepcopy
from pandas import DataFrame, set_option

set_option('display.max_columns', None)
set_option('display.width', 1000)
set_option('display.colheader_justify', 'center')
set_option('display.precision', 3)


_ASSOCIATION_COLUMNS = [
    'association_id',
    'risk_frequency',
    'pvalue_description',
    'pvalue_mantissa',
    'pvalue_exponent',
    'multi_snp_haplotype',
    'snp_interaction',
    'snp_type',
    'range',
    'description',
    'or_per_copy_num',
    'or_value',
    'beta_num',
    'beta_unit',
    'beta_direction',
    'beta',
    'last_mapping_date',
    'last_update_date',
    'p_value',
    'pubmed_id',
    'first_author',
    'ci_lower',
    'ci_upper'
]

_EFO_TRAITS_COLUMNS = ['association_id', 'efo_id', 'efo_trait']
_SNP_ALLELE_COLUMNS = ['association_id', 'rs_id', 'effect_allele']
_BG_EFO_TRAITS_COLUMNS = ['association_id', 'efo_id', 'efo_trait']
_REPORTED_TRAITS_COLUMNS = ['association_id', 'reported_trait']
_LOCATIONS_COLUMNS = ['association_id', 'locations']
_MAPPED_GENES_COLUMNS = ['association_id', 'mapped_genes']
_SNP_EFFECT_ALLELES_COLUMNS = ['association_id', 'snp_effect_allele']


_NESTED_FIELDS = [
    'efo_traits',
    'snp_allele',
    'bg_efo_traits',
    'reported_trait',
    'locations',
    'mapped_genes',
    'snp_effect_allele'
]


class Association:
    """
       An object that stores data of type Association. To understand the significance of each column of the DataFrame. Please visit "Association" in [Gwas Catalog Documentation](https://www.ebi.ac.uk/gwas/rest/api/v2/docs/reference) for details.

       Attributes:
            raw_data: list. Converted from obtained JSON data.
            associations: DataFrame. It only exists if the parameter mode of constructor is Fat.
            efo_traits: DataFrame. It only exists if the parameter mode of constructor is Fat.
            snp_allele: DataFrame. It only exists if the parameter mode of constructor is Fat.
            bg_efo_traits: DataFrame. It only exists if the parameter mode of constructor is Fat.
            reported_traits: DataFrame. It only exists if the parameter mode of constructor is Fat.
            locations: DataFrame. It only exists if the parameter mode of constructor is Fat.
            mapped_genes: DataFrame. It only exists if the parameter mode of constructor is Fat.
            snp_effect_alleles: DataFrame. It only exists if the parameter mode of constructor is Fat.
            mode: Fat or Thin. Specifies the mode of the returned object.

       ```python
       from pandasgwas.get_associations import get_associations

       ch = get_associations(association_id='169106883')
       ch
       ch.raw_data
       ch.mode
       ch.associations
       ch.efo_traits
       ch.snp_allele
       ch.bg_efo_traits
       ch.reported_traits
       ch.locations
       ch.mapped_genes
       ch.snp_effect_alleles
       ```
       Subset object s by either identifier or position
       ```python
       all_df = get_associations()
       all_df[0].associations
       all_df[0:3].associations
       all_df['169106883'].associations
       all_df[('169106883','218958492','218958441')].associations
       ```
       Objects can be manipulated like sets in the mathematical sense.
       ```python
       one = get_associations(association_id='169106883')
       two = get_associations(association_id='218958492')
       three = get_associations(association_id='218958441')
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
                if 'association_id' in rec and rec['association_id'] is not None:
                    rec['association_id'] = str(rec['association_id'])
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

        # Fat 模式：先标准化字段名、去掉 links/_links、补齐缺失字段
        standardized_data = []
        for item in data:
            if not isinstance(item, dict):
                continue

            rec = deepcopy(item)

            normalized = {}
            for key, value in rec.items():
                normalized[key] = value

            standardized = {}
            for col in _ASSOCIATION_COLUMNS:
                value = normalized.get(col, None)
                if col == 'range':
                    if value is None:
                        value = None
                    elif isinstance(value, list):
                        if len(value) == 0:
                            value = None
                        elif len(value) == 1:
                            value = None if value[0] is None else str(value[0])
                        else:
                            value = '; '.join('' if v is None else str(v) for v in value)
                    else:
                        value = str(value)
                standardized[col] = value

            for field in _NESTED_FIELDS:
                standardized[field] = normalized.get(field, None)

            standardized_data.append(standardized)

        self.raw_data = standardized_data

        if len(standardized_data) == 0:
            self.associations = DataFrame(columns=_ASSOCIATION_COLUMNS)
            self.efo_traits = DataFrame(columns=_EFO_TRAITS_COLUMNS)
            self.snp_allele = DataFrame(columns=_SNP_ALLELE_COLUMNS)
            self.bg_efo_traits = DataFrame(columns=_BG_EFO_TRAITS_COLUMNS)
            self.reported_traits = DataFrame(columns=_REPORTED_TRAITS_COLUMNS)
            self.locations = DataFrame(columns=_LOCATIONS_COLUMNS)
            self.mapped_genes = DataFrame(columns=_MAPPED_GENES_COLUMNS)
            self.snp_effect_alleles = DataFrame(columns=_SNP_EFFECT_ALLELES_COLUMNS)
            return

        association_rows = []
        efo_traits_rows = []
        snp_allele_rows = []
        bg_efo_traits_rows = []
        reported_traits_rows = []
        locations_rows = []
        mapped_genes_rows = []
        snp_effect_alleles_rows = []

        for record in standardized_data:
            association_id = record.get('association_id', None)

            # 主表：固定字段，缺失补 None
            association_rows.append({
                col: record.get(col, None) for col in _ASSOCIATION_COLUMNS
            })

            # efo_traits
            raw_efo_traits = record.get('efo_traits', None)
            if raw_efo_traits is None:
                iterable_efo_traits = []
            elif isinstance(raw_efo_traits, list):
                iterable_efo_traits = raw_efo_traits
            else:
                iterable_efo_traits = [raw_efo_traits]
            for item in iterable_efo_traits:
                if isinstance(item, dict):
                    efo_traits_rows.append({
                        'association_id': association_id,
                        'efo_id': item.get('efo_id', None),
                        'efo_trait': item.get('efo_trait', None)
                    })
                else:
                    efo_traits_rows.append({
                        'association_id': association_id,
                        'efo_id': None,
                        'efo_trait': item
                    })

            # snp_allele
            raw_snp_allele = record.get('snp_allele', None)
            if raw_snp_allele is None:
                iterable_snp_allele = []
            elif isinstance(raw_snp_allele, list):
                iterable_snp_allele = raw_snp_allele
            else:
                iterable_snp_allele = [raw_snp_allele]
            for item in iterable_snp_allele:
                if isinstance(item, dict):
                    snp_allele_rows.append({
                        'association_id': association_id,
                        'rs_id': item.get('rs_id', None),
                        'effect_allele': item.get('effect_allele', None)
                    })
                else:
                    snp_allele_rows.append({
                        'association_id': association_id,
                        'rs_id': item,
                        'effect_allele': None
                    })

            # bg_efo_traits
            raw_bg_efo_traits = record.get('bg_efo_traits', None)
            if raw_bg_efo_traits is None:
                iterable_bg_efo_traits = []
            elif isinstance(raw_bg_efo_traits, list):
                iterable_bg_efo_traits = raw_bg_efo_traits
            else:
                iterable_bg_efo_traits = [raw_bg_efo_traits]
            for item in iterable_bg_efo_traits:
                if isinstance(item, dict):
                    bg_efo_traits_rows.append({
                        'association_id': association_id,
                        'efo_id': item.get('efo_id', None),
                        'efo_trait': item.get('efo_trait', None)
                    })
                else:
                    bg_efo_traits_rows.append({
                        'association_id': association_id,
                        'efo_id': None,
                        'efo_trait': item
                    })

            # reported_traits（属性名保持不变，数据字段统一用 reported_trait）
            raw_reported_trait = record.get('reported_trait', None)
            if raw_reported_trait is None:
                iterable_reported_trait = []
            elif isinstance(raw_reported_trait, list):
                iterable_reported_trait = raw_reported_trait
            else:
                iterable_reported_trait = [raw_reported_trait]
            for item in iterable_reported_trait:
                reported_traits_rows.append({
                    'association_id': association_id,
                    'reported_trait': item
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
                locations_rows.append({
                    'association_id': association_id,
                    'locations': item
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
                    'association_id': association_id,
                    'mapped_genes': item
                })

            # snp_effect_alleles（属性名保持不变，数据字段统一用 snp_effect_allele）
            raw_snp_effect_allele = record.get('snp_effect_allele', None)
            if raw_snp_effect_allele is None:
                iterable_snp_effect_allele = []
            elif isinstance(raw_snp_effect_allele, list):
                iterable_snp_effect_allele = raw_snp_effect_allele
            else:
                iterable_snp_effect_allele = [raw_snp_effect_allele]
            for item in iterable_snp_effect_allele:
                snp_effect_alleles_rows.append({
                    'association_id': association_id,
                    'snp_effect_allele': item
                })

        self.associations = DataFrame(association_rows).reindex(columns=_ASSOCIATION_COLUMNS) if len(association_rows) else DataFrame(columns=_ASSOCIATION_COLUMNS)
        self.efo_traits = DataFrame(efo_traits_rows).reindex(columns=_EFO_TRAITS_COLUMNS) if len(efo_traits_rows) else DataFrame(columns=_EFO_TRAITS_COLUMNS)
        self.snp_alleles = DataFrame(snp_allele_rows).reindex(columns=_SNP_ALLELE_COLUMNS) if len(snp_allele_rows) else DataFrame(columns=_SNP_ALLELE_COLUMNS)
        self.bg_efo_traits = DataFrame(bg_efo_traits_rows).reindex(columns=_BG_EFO_TRAITS_COLUMNS) if len(bg_efo_traits_rows) else DataFrame(columns=_BG_EFO_TRAITS_COLUMNS)
        self.reported_traits = DataFrame(reported_traits_rows).reindex(columns=_REPORTED_TRAITS_COLUMNS) if len(reported_traits_rows) else DataFrame(columns=_REPORTED_TRAITS_COLUMNS)
        self.locations = DataFrame(locations_rows).reindex(columns=_LOCATIONS_COLUMNS) if len(locations_rows) else DataFrame(columns=_LOCATIONS_COLUMNS)
        self.mapped_genes = DataFrame(mapped_genes_rows).reindex(columns=_MAPPED_GENES_COLUMNS) if len(mapped_genes_rows) else DataFrame(columns=_MAPPED_GENES_COLUMNS)
        self.snp_effect_alleles = DataFrame(snp_effect_alleles_rows).reindex(columns=_SNP_EFFECT_ALLELES_COLUMNS) if len(snp_effect_alleles_rows) else DataFrame(columns=_SNP_EFFECT_ALLELES_COLUMNS)

    def __str__(self):
        if self.mode == 'Fat':
            return (
                    "Association is running in fat mode. It has 8 DataFrames with hierarchical dependencies.\n"
                    "associations: %d rows\n"
                    "|\n -efo_traits: %d rows\n"
                    "|\n -snp_allele: %d rows\n"
                    "|\n -bg_efo_traits: %d rows\n"
                    "|\n -reported_traits: %d rows\n"
                    "|\n -locations: %d rows\n"
                    "|\n -mapped_genes: %d rows\n"
                    "|\n -snp_effect_alleles: %d rows"
                    % (
                        len(self.associations),
                        len(self.efo_traits),
                        len(self.snp_allele),
                        len(self.bg_efo_traits),
                        len(self.reported_traits),
                        len(self.locations),
                        len(self.mapped_genes),
                        len(self.snp_effect_alleles),
                    )
            )
        if self.mode == 'Thin':
            return (
                    'Association is running in thin mode. It has 1 list that contains the raw data.\n'
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
            raw_data_dict[j['association_id']] = j

        sub_set = []
        for i in arr:
            if isinstance(i, str):
                sub_set.append(raw_data_dict[i])
            elif isinstance(i, int):
                sub_set.append(raw_data[i])
            else:
                raise TypeError('Invalid item type: {}'.format(type(i)))

        return Association(sub_set, self.mode)

    def __add__(self, other):
        if self.mode == other.mode:
            return Association(self.raw_data + other.raw_data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __sub__(self, other):
        if self.mode == other.mode:
            self_key_set = set()
            self_dict = {}
            other_key_set = set()

            for i in self.raw_data:
                self_key_set.add(i['association_id'])
                self_dict[i['association_id']] = i

            for j in other.raw_data:
                other_key_set.add(j['association_id'])

            sub_key = self_key_set - other_key_set
            data = []
            for k in sub_key:
                data.append(self_dict[k])

            return Association(data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __and__(self, other):
        if self.mode == other.mode:
            self_key_set = set()
            self_dict = {}
            other_key_set = set()

            for i in self.raw_data:
                self_key_set.add(i['association_id'])
                self_dict[i['association_id']] = i

            for j in other.raw_data:
                other_key_set.add(j['association_id'])

            sub_key = self_key_set & other_key_set
            data = []
            for k in sub_key:
                data.append(self_dict[k])

            return Association(data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __or__(self, other):
        if self.mode == other.mode:
            and_dict = {}

            for i in self.raw_data:
                and_dict[i['association_id']] = i

            for j in other.raw_data:
                and_dict[j['association_id']] = j

            data = list(and_dict.values())
            return Association(data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __xor__(self, other):
        if self.mode == other.mode:
            self_key_set = set()
            and_dict = {}
            other_key_set = set()

            for i in self.raw_data:
                self_key_set.add(i['association_id'])
                and_dict[i['association_id']] = i

            for j in other.raw_data:
                other_key_set.add(j['association_id'])
                and_dict[j['association_id']] = j

            sub_key = self_key_set ^ other_key_set
            data = []
            for k in sub_key:
                data.append(and_dict[k])

            return Association(data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __eq__(self, other):
        if self is None or other is None:
            return self is None and other is None
        return self.raw_data == other.raw_data and self.mode == other.mode

    def __len__(self):
        return len(self.raw_data)



