from copy import deepcopy
from pandas import DataFrame, set_option

set_option('display.max_columns', None)
set_option('display.width', 1000)
set_option('display.colheader_justify', 'center')
set_option('display.precision', 3)


_ANCESTRY_COLUMNS = [
    'accession_id',
    'ancestry_id',
    'type',
    'number_of_individuals'
]

_ANCESTRAL_GROUPS_COLUMNS = ['accession_id', 'ancestry_id', 'ancestral_group']
_COUNTRY_OF_ORIGINS_COLUMNS = ['accession_id', 'ancestry_id', 'major_area', 'region', 'country_name']
_COUNTRY_OF_RECRUITMENTS_COLUMNS = ['accession_id', 'ancestry_id', 'major_area', 'region', 'country_name']

_NESTED_FIELDS = [
    'ancestral_groups',
    'country_of_origin',
    'country_of_recruitment'
]


class Ancestry:
    """
       An object that stores data of type Ancestry. To understand the significance of each column of the DataFrame. Please visit "Ancestry" in [Gwas Catalog Documentation](https://www.ebi.ac.uk/gwas/rest/api/v2/docs/reference) for details.

       Attributes:
            raw_data: list. Converted from obtained JSON data.
            ancestries: DataFrame. It only exists if the parameter mode of constructor is Fat.
            ancestral_groups: DataFrame. It only exists if the parameter mode of constructor is Fat.
            country_of_origins: DataFrame. It only exists if the parameter mode of constructor is Fat.
            country_of_recruitments: DataFrame. It only exists if the parameter mode of constructor is Fat.
            mode: Fat or Thin. Specifies the mode of the returned object.

       ```python
       from pandasgwas.get_ancestries import get_ancestries

       ch = get_ancestries(accession_id='GCST000854')
       ch
       ch.raw_data
       ch.mode
       ch.ancestries
       ch.ancestral_groups
       ch.country_of_origins
       ch.country_of_recruitments
       ```
       Subset object s by either identifier or position
       ```python
       all_df = get_ancestries(accession_id='GCST000854')
       all_df[0].ancestries
       all_df[0:3].ancestries
       all_df['GCST000854|4595'].ancestries
       all_df[('GCST000854|4595','GCST000854|7528','GCST000854|7529')].ancestries
       ```
       Objects can be manipulated like sets in the mathematical sense.
       ```python
       one = get_ancestries(accession_id='GCST000854', ancestry_id='4595')
       two = get_ancestries(accession_id='GCST000854', ancestry_id='7528')
       three = get_ancestries(accession_id='GCST000854', ancestry_id='7529')
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
                if 'accession_id' in rec and rec['accession_id'] is not None:
                    rec['accession_id'] = str(rec['accession_id'])
                if 'ancestry_id' in rec and rec['ancestry_id'] is not None:
                    rec['ancestry_id'] = str(rec['ancestry_id'])
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
            for col in _ANCESTRY_COLUMNS:
                value = normalized.get(col, None)
                standardized[col] = value

            for field in _NESTED_FIELDS:
                standardized[field] = normalized.get(field, None)

            standardized_data.append(standardized)

        self.raw_data = standardized_data

        if len(standardized_data) == 0:
            self.ancestries = DataFrame(columns=_ANCESTRY_COLUMNS)
            self.ancestral_groups = DataFrame(columns=_ANCESTRAL_GROUPS_COLUMNS)
            self.country_of_origins = DataFrame(columns=_COUNTRY_OF_ORIGINS_COLUMNS)
            self.country_of_recruitments = DataFrame(columns=_COUNTRY_OF_RECRUITMENTS_COLUMNS)
            return

        ancestry_rows = []
        ancestral_groups_rows = []
        country_of_origins_rows = []
        country_of_recruitments_rows = []

        for record in standardized_data:
            accession_id = record.get('accession_id', None)
            ancestry_id = record.get('ancestry_id', None)

            # main table: fixed fields, missing fields filled with None
            ancestry_rows.append({
                col: record.get(col, None) for col in _ANCESTRY_COLUMNS
            })

            # ancestral_groups
            raw_ancestral_groups = record.get('ancestral_groups', None)
            if raw_ancestral_groups is None:
                iterable_ancestral_groups = []
            elif isinstance(raw_ancestral_groups, list):
                iterable_ancestral_groups = raw_ancestral_groups
            else:
                iterable_ancestral_groups = [raw_ancestral_groups]
            for item in iterable_ancestral_groups:
                if isinstance(item, dict):
                    ancestral_groups_rows.append({
                        'accession_id': accession_id,
                        'ancestry_id': ancestry_id,
                        'ancestral_group': item.get('ancestral_group', None)
                    })
                else:
                    ancestral_groups_rows.append({
                        'accession_id': accession_id,
                        'ancestry_id': ancestry_id,
                        'ancestral_group': item
                    })

            # country_of_origins
            raw_country_of_origins = record.get('country_of_origin', None)
            if raw_country_of_origins is None:
                iterable_country_of_origins = []
            elif isinstance(raw_country_of_origins, list):
                iterable_country_of_origins = raw_country_of_origins
            else:
                iterable_country_of_origins = [raw_country_of_origins]
            for item in iterable_country_of_origins:
                if isinstance(item, dict):
                    country_of_origins_rows.append({
                        'accession_id': accession_id,
                        'ancestry_id': ancestry_id,
                        'major_area': item.get('major_area', None),
                        'region': item.get('region', None),
                        'country_name': item.get('country_name', None)
                    })
                else:
                    country_of_origins_rows.append({
                        'accession_id': accession_id,
                        'ancestry_id': ancestry_id,
                        'major_area': None,
                        'region': None,
                        'country_name': item
                    })

            # country_of_recruitments
            raw_country_of_recruitments = record.get('country_of_recruitment', None)
            if raw_country_of_recruitments is None:
                iterable_country_of_recruitments = []
            elif isinstance(raw_country_of_recruitments, list):
                iterable_country_of_recruitments = raw_country_of_recruitments
            else:
                iterable_country_of_recruitments = [raw_country_of_recruitments]
            for item in iterable_country_of_recruitments:
                if isinstance(item, dict):
                    country_of_recruitments_rows.append({
                        'accession_id': accession_id,
                        'ancestry_id': ancestry_id,
                        'major_area': item.get('major_area', None),
                        'region': item.get('region', None),
                        'country_name': item.get('country_name', None)
                    })
                else:
                    country_of_recruitments_rows.append({
                        'accession_id': accession_id,
                        'ancestry_id': ancestry_id,
                        'major_area': None,
                        'region': None,
                        'country_name': item
                    })

        self.ancestries = DataFrame(ancestry_rows).reindex(columns=_ANCESTRY_COLUMNS) if len(ancestry_rows) else DataFrame(columns=_ANCESTRY_COLUMNS)
        self.ancestral_groups = DataFrame(ancestral_groups_rows).reindex(columns=_ANCESTRAL_GROUPS_COLUMNS) if len(ancestral_groups_rows) else DataFrame(columns=_ANCESTRAL_GROUPS_COLUMNS)
        self.country_of_origins = DataFrame(country_of_origins_rows).reindex(columns=_COUNTRY_OF_ORIGINS_COLUMNS) if len(country_of_origins_rows) else DataFrame(columns=_COUNTRY_OF_ORIGINS_COLUMNS)
        self.country_of_recruitments = DataFrame(country_of_recruitments_rows).reindex(columns=_COUNTRY_OF_RECRUITMENTS_COLUMNS) if len(country_of_recruitments_rows) else DataFrame(columns=_COUNTRY_OF_RECRUITMENTS_COLUMNS)

    def __str__(self):
        if self.mode == 'Fat':
            return (
                    "Ancestry is running in fat mode. It has 4 DataFrames with hierarchical dependencies.\n"
                    "ancestries: %d rows\n"
                    "|\n -ancestral_groups: %d rows\n"
                    "|\n -country_of_origins: %d rows\n"
                    "|\n -country_of_recruitments: %d rows"
                    % (
                        len(self.ancestries),
                        len(self.ancestral_groups),
                        len(self.country_of_origins),
                        len(self.country_of_recruitments),
                    )
            )
        if self.mode == 'Thin':
            return (
                    'Ancestry is running in thin mode. It has 1 list that contains the raw data.\n'
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
            raw_data_dict['{}|{}'.format(j['accession_id'], j['ancestry_id'])] = j

        sub_set = []
        for i in arr:
            if isinstance(i, str):
                sub_set.append(raw_data_dict[i])
            elif isinstance(i, int):
                sub_set.append(raw_data[i])
            else:
                raise TypeError('Invalid item type: {}'.format(type(i)))

        return Ancestry(sub_set, self.mode)

    def __add__(self, other):
        if self.mode == other.mode:
            return Ancestry(self.raw_data + other.raw_data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __sub__(self, other):
        if self.mode == other.mode:
            self_key_set = set()
            self_dict = {}
            other_key_set = set()

            for i in self.raw_data:
                key = '{}/{}'.format(i['accession_id'], i['ancestry_id'])
                self_key_set.add(key)
                self_dict[key] = i

            for j in other.raw_data:
                other_key_set.add('{}/{}'.format(j['accession_id'], j['ancestry_id']))

            sub_key = self_key_set - other_key_set
            data = []
            for k in sub_key:
                data.append(self_dict[k])

            return Ancestry(data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __and__(self, other):
        if self.mode == other.mode:
            self_key_set = set()
            self_dict = {}
            other_key_set = set()

            for i in self.raw_data:
                key = '{}/{}'.format(i['accession_id'], i['ancestry_id'])
                self_key_set.add(key)
                self_dict[key] = i

            for j in other.raw_data:
                other_key_set.add('{}/{}'.format(j['accession_id'], j['ancestry_id']))

            sub_key = self_key_set & other_key_set
            data = []
            for k in sub_key:
                data.append(self_dict[k])

            return Ancestry(data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __or__(self, other):
        if self.mode == other.mode:
            and_dict = {}

            for i in self.raw_data:
                and_dict['{}/{}'.format(i['accession_id'], i['ancestry_id'])] = i

            for j in other.raw_data:
                and_dict['{}/{}'.format(j['accession_id'], j['ancestry_id'])] = j

            data = list(and_dict.values())
            return Ancestry(data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __xor__(self, other):
        if self.mode == other.mode:
            self_key_set = set()
            and_dict = {}
            other_key_set = set()

            for i in self.raw_data:
                key = '{}/{}'.format(i['accession_id'], i['ancestry_id'])
                self_key_set.add(key)
                and_dict[key] = i

            for j in other.raw_data:
                key = '{}/{}'.format(j['accession_id'], j['ancestry_id'])
                other_key_set.add(key)
                and_dict[key] = j

            sub_key = self_key_set ^ other_key_set
            data = []
            for k in sub_key:
                data.append(and_dict[k])

            return Ancestry(data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __eq__(self, other):
        if self is None or other is None:
            return self is None and other is None
        return self.raw_data == other.raw_data and self.mode == other.mode

    def __len__(self):
        return len(self.raw_data)
