from copy import deepcopy
from pandas import DataFrame, set_option

set_option('display.max_columns', None)
set_option('display.width', 1000)
set_option('display.colheader_justify', 'center')
set_option('display.precision', 3)


_TRAIT_COLUMNS = [
    'efo_trait',
    'uri',
    'efo_id'
]


class Trait:
    """
       An object that stores data of type Trait. To understand the significance of each column of the DataFrame. Please visit "Trait" in [Gwas Catalog Documentation](https://www.ebi.ac.uk/gwas/rest/api/v2/docs/reference) for details.

       Attributes:
            raw_data: list. Converted from obtained JSON data.
            traits: DataFrame. It only exists if the parameter mode of constructor is Fat.
            mode: Fat or Thin. Specifies the mode of the returned object.

       ```python
       from pandasgwas.get_traits import get_traits

       ch = get_traits(efo_id='EFO_0009713')
       ch
       ch.raw_data
       ch.mode
       ch.traits
       ```
       Subset object s by either identifier or position
       ```python
       all_df = get_traits()
       all_df[0].traits
       all_df[0:3].traits
       all_df['EFO_0009713'].traits
       all_df[('EFO_0009713','OBA_VT0001321','EFO_0920093')].traits
       ```
       Objects can be manipulated like sets in the mathematical sense.
       ```python
       one = get_traits(efo_id='EFO_0009713')
       two = get_traits(efo_id='OBA_VT0001321')
       three = get_traits(efo_id='EFO_0920093')
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
                if 'efo_id' in rec and rec['efo_id'] is not None:
                    rec['efo_id'] = str(rec['efo_id'])
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
            for col in _TRAIT_COLUMNS:
                value = normalized.get(col, None)
                standardized[col] = value

            standardized_data.append(standardized)

        self.raw_data = standardized_data

        if len(standardized_data) == 0:
            self.traits = DataFrame(columns=_TRAIT_COLUMNS)
            return

        trait_rows = []

        for record in standardized_data:
            trait_rows.append({
                col: record.get(col, None) for col in _TRAIT_COLUMNS
            })

        self.traits = DataFrame(trait_rows).reindex(columns=_TRAIT_COLUMNS) if len(trait_rows) else DataFrame(columns=_TRAIT_COLUMNS)

    def __str__(self):
        if self.mode == 'Fat':
            return (
                    "Trait is running in fat mode. It has 1 DataFrame with hierarchical dependencies.\n"
                    "traits: %d rows"
                    % len(self.traits)
            )
        if self.mode == 'Thin':
            return (
                    'Trait is running in thin mode. It has 1 list that contains the raw data.\n'
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
            raw_data_dict[j['efo_id']] = j

        sub_set = []
        for i in arr:
            if isinstance(i, str):
                sub_set.append(raw_data_dict[i])
            elif isinstance(i, int):
                sub_set.append(raw_data[i])
            else:
                raise TypeError('Invalid item type: {}'.format(type(i)))

        return Trait(sub_set, self.mode)

    def __add__(self, other):
        if self.mode == other.mode:
            return Trait(self.raw_data + other.raw_data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __sub__(self, other):
        if self.mode == other.mode:
            self_key_set = set()
            self_dict = {}
            other_key_set = set()

            for i in self.raw_data:
                self_key_set.add(i['efo_id'])
                self_dict[i['efo_id']] = i

            for j in other.raw_data:
                other_key_set.add(j['efo_id'])

            sub_key = self_key_set - other_key_set
            data = []
            for k in sub_key:
                data.append(self_dict[k])

            return Trait(data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __and__(self, other):
        if self.mode == other.mode:
            self_key_set = set()
            self_dict = {}
            other_key_set = set()

            for i in self.raw_data:
                self_key_set.add(i['efo_id'])
                self_dict[i['efo_id']] = i

            for j in other.raw_data:
                other_key_set.add(j['efo_id'])

            sub_key = self_key_set & other_key_set
            data = []
            for k in sub_key:
                data.append(self_dict[k])

            return Trait(data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __or__(self, other):
        if self.mode == other.mode:
            and_dict = {}

            for i in self.raw_data:
                and_dict[i['efo_id']] = i

            for j in other.raw_data:
                and_dict[j['efo_id']] = j

            data = list(and_dict.values())
            return Trait(data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __xor__(self, other):
        if self.mode == other.mode:
            self_key_set = set()
            and_dict = {}
            other_key_set = set()

            for i in self.raw_data:
                self_key_set.add(i['efo_id'])
                and_dict[i['efo_id']] = i

            for j in other.raw_data:
                other_key_set.add(j['efo_id'])
                and_dict[j['efo_id']] = j

            sub_key = self_key_set ^ other_key_set
            data = []
            for k in sub_key:
                data.append(and_dict[k])

            return Trait(data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __eq__(self, other):
        if self is None or other is None:
            return self is None and other is None
        return self.raw_data == other.raw_data and self.mode == other.mode

    def __len__(self):
        return len(self.raw_data)
