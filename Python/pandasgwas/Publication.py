from copy import deepcopy
from pandas import DataFrame, set_option

set_option('display.max_columns', None)
set_option('display.width', 1000)
set_option('display.colheader_justify', 'center')
set_option('display.precision', 3)


_PUBLICATION_COLUMNS = [
    'pubmed_id',
    'publication_date',
    'journal',
    'title',
    'first_author.full_name',
    'first_author.full_name_standard',
    'first_author.first_name',
    'first_author.last_name',
    'first_author.initials',
    'first_author.affiliation',
    'first_author.orcid'
]

_AUTHORS_COLUMNS = [
    'pubmed_id',
    'full_name',
    'full_name_standard',
    'first_name',
    'last_name',
    'initials',
    'affiliation',
    'orcid'
]

_NESTED_FIELDS = [
    'authors'
]


class Publication:
    """
       An object that stores data of type Publication. To understand the significance of each column of the DataFrame. Please visit "Publication" in [Gwas Catalog Documentation](https://www.ebi.ac.uk/gwas/rest/api/v2/docs/reference) for details.

       Attributes:
            raw_data: list. Converted from obtained JSON data.
            publications: DataFrame. It only exists if the parameter mode of constructor is Fat.
            authors: DataFrame. It only exists if the parameter mode of constructor is Fat.
            mode: Fat or Thin. Specifies the mode of the returned object.

       ```python
       from pandasgwas.get_publications import get_publications

       ch = get_publications(pubmed_id='35241825')
       ch
       ch.raw_data
       ch.mode
       ch.publications
       ch.authors
       ```
       Subset object s by either identifier or position
       ```python
       all_df = get_publications()
       all_df[0].publications
       all_df[0:3].publications
       all_df['35241825'].publications
       all_df[('35241825','42127933','42107849')].publications
       ```
       Objects can be manipulated like sets in the mathematical sense.
       ```python
       one = get_publications(pubmed_id='35241825')
       two = get_publications(pubmed_id='42127933')
       three = get_publications(pubmed_id='42107849')
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
                if 'pubmed_id' in rec and rec['pubmed_id'] is not None:
                    rec['pubmed_id'] = str(rec['pubmed_id'])
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

            raw_first_author = normalized.get('first_author', None)
            if isinstance(raw_first_author, list):
                if len(raw_first_author) == 0:
                    first_author = None
                else:
                    first_author = raw_first_author[0]
            else:
                first_author = raw_first_author

            standardized = {}
            for col in _PUBLICATION_COLUMNS:
                if col.startswith('first_author.'):
                    author_col = col.split('.', 1)[1]
                    if col in normalized:
                        value = normalized.get(col, None)
                    elif isinstance(first_author, dict):
                        value = first_author.get(author_col, None)
                    elif author_col == 'full_name':
                        value = first_author
                    else:
                        value = None
                else:
                    value = normalized.get(col, None)
                standardized[col] = value

            for field in _NESTED_FIELDS:
                standardized[field] = normalized.get(field, None)

            standardized_data.append(standardized)

        self.raw_data = standardized_data

        if len(standardized_data) == 0:
            self.publications = DataFrame(columns=_PUBLICATION_COLUMNS)
            self.authors = DataFrame(columns=_AUTHORS_COLUMNS)
            return

        publication_rows = []
        authors_rows = []

        for record in standardized_data:
            pubmed_id = record.get('pubmed_id', None)

            publication_rows.append({
                col: record.get(col, None) for col in _PUBLICATION_COLUMNS
            })

            raw_authors = record.get('authors', None)
            if raw_authors is None:
                iterable_authors = []
            elif isinstance(raw_authors, list):
                iterable_authors = raw_authors
            else:
                iterable_authors = [raw_authors]
            for item in iterable_authors:
                if isinstance(item, dict):
                    authors_rows.append({
                        'pubmed_id': pubmed_id,
                        'full_name': item.get('full_name', None),
                        'full_name_standard': item.get('full_name_standard', None),
                        'first_name': item.get('first_name', None),
                        'last_name': item.get('last_name', None),
                        'initials': item.get('initials', None),
                        'affiliation': item.get('affiliation', None),
                        'orcid': item.get('orcid', None)
                    })
                else:
                    authors_rows.append({
                        'pubmed_id': pubmed_id,
                        'full_name': item,
                        'full_name_standard': None,
                        'first_name': None,
                        'last_name': None,
                        'initials': None,
                        'affiliation': None,
                        'orcid': None
                    })

        self.publications = DataFrame(publication_rows).reindex(columns=_PUBLICATION_COLUMNS) if len(publication_rows) else DataFrame(columns=_PUBLICATION_COLUMNS)
        self.authors = DataFrame(authors_rows).reindex(columns=_AUTHORS_COLUMNS) if len(authors_rows) else DataFrame(columns=_AUTHORS_COLUMNS)

    def __str__(self):
        if self.mode == 'Fat':
            return (
                    "Publication is running in fat mode. It has 2 DataFrames with hierarchical dependencies.\n"
                    "publications: %d rows\n"
                    "|\n -authors: %d rows"
                    % (
                        len(self.publications),
                        len(self.authors),
                    )
            )
        if self.mode == 'Thin':
            return (
                    'Publication is running in thin mode. It has 1 list that contains the raw data.\n'
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
            raw_data_dict[j['pubmed_id']] = j

        sub_set = []
        for i in arr:
            if isinstance(i, str):
                sub_set.append(raw_data_dict[i])
            elif isinstance(i, int):
                sub_set.append(raw_data[i])
            else:
                raise TypeError('Invalid item type: {}'.format(type(i)))

        return Publication(sub_set, self.mode)

    def __add__(self, other):
        if self.mode == other.mode:
            return Publication(self.raw_data + other.raw_data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __sub__(self, other):
        if self.mode == other.mode:
            self_key_set = set()
            self_dict = {}
            other_key_set = set()

            for i in self.raw_data:
                self_key_set.add(i['pubmed_id'])
                self_dict[i['pubmed_id']] = i

            for j in other.raw_data:
                other_key_set.add(j['pubmed_id'])

            sub_key = self_key_set - other_key_set
            data = []
            for k in sub_key:
                data.append(self_dict[k])

            return Publication(data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __and__(self, other):
        if self.mode == other.mode:
            self_key_set = set()
            self_dict = {}
            other_key_set = set()

            for i in self.raw_data:
                self_key_set.add(i['pubmed_id'])
                self_dict[i['pubmed_id']] = i

            for j in other.raw_data:
                other_key_set.add(j['pubmed_id'])

            sub_key = self_key_set & other_key_set
            data = []
            for k in sub_key:
                data.append(self_dict[k])

            return Publication(data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __or__(self, other):
        if self.mode == other.mode:
            and_dict = {}

            for i in self.raw_data:
                and_dict[i['pubmed_id']] = i

            for j in other.raw_data:
                and_dict[j['pubmed_id']] = j

            data = list(and_dict.values())
            return Publication(data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __xor__(self, other):
        if self.mode == other.mode:
            self_key_set = set()
            and_dict = {}
            other_key_set = set()

            for i in self.raw_data:
                self_key_set.add(i['pubmed_id'])
                and_dict[i['pubmed_id']] = i

            for j in other.raw_data:
                other_key_set.add(j['pubmed_id'])
                and_dict[j['pubmed_id']] = j

            sub_key = self_key_set ^ other_key_set
            data = []
            for k in sub_key:
                data.append(and_dict[k])

            return Publication(data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __eq__(self, other):
        if self is None or other is None:
            return self is None and other is None
        return self.raw_data == other.raw_data and self.mode == other.mode

    def __len__(self):
        return len(self.raw_data)
