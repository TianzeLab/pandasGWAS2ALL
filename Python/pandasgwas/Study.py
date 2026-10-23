from copy import deepcopy
from pandas import DataFrame, set_option

set_option('display.max_columns', None)
set_option('display.width', 1000)
set_option('display.colheader_justify', 'center')
set_option('display.precision', 3)


_STUDY_COLUMNS = [
    'initial_sample_size',
    'replication_sample_size',
    'study_id',
    'gxe',
    'gxg',
    'snp_count',
    'qualifier',
    'imputed',
    'pooled',
    'accession_id',
    'user_requested',
    'pubmed_id',
    'platforms',
    'disease_trait',
    'terms_of_license',
    'full_summary_stats',
    'full_summary_stats_available'
]

_EFO_TRAITS_COLUMNS = ['accession_id', 'efo_id', 'efo_trait']
_GENOTYPING_TECHNOLOGIES_COLUMNS = ['accession_id', 'genotyping_technologies']
_DISCOVERY_ANCESTRIES_COLUMNS = ['accession_id', 'discovery_ancestry']
_REPLICATION_ANCESTRIES_COLUMNS = ['accession_id', 'replication_ancestry']
_BG_EFO_TRAITS_COLUMNS = ['accession_id', 'efo_id', 'efo_trait']
_COHORTS_COLUMNS = ['accession_id', 'cohort']
_ARRAY_MANUFACTURERS_COLUMNS = ['accession_id', 'array_manufacturer']


_NESTED_FIELDS = [
    'efo_traits',
    'genotyping_technologies',
    'discovery_ancestry',
    'replication_ancestry',
    'bg_efo_traits',
    'cohort',
    'array_manufacturer'
]


class Study:
    """
       An object that stores data of type Study. To understand the significance of each column of the DataFrame. Please visit "Study" in [Gwas Catalog Documentation](https://www.ebi.ac.uk/gwas/rest/api/v2/docs/reference) for details.

       Attributes:
            raw_data: list. Converted from obtained JSON data.
            studies: DataFrame. It only exists if the parameter mode of constructor is Fat.
            efo_traits: DataFrame. It only exists if the parameter mode of constructor is Fat.
            genotyping_technologies: DataFrame. It only exists if the parameter mode of constructor is Fat.
            discovery_ancestries: DataFrame. It only exists if the parameter mode of constructor is Fat.
            replication_ancestries: DataFrame. It only exists if the parameter mode of constructor is Fat.
            bg_efo_traits: DataFrame. It only exists if the parameter mode of constructor is Fat.
            cohorts: DataFrame. It only exists if the parameter mode of constructor is Fat.
            array_manufacturers: DataFrame. It only exists if the parameter mode of constructor is Fat.
            mode: Fat or Thin. Specifies the mode of the returned object.

       ```python
       from pandasgwas.get_studies import get_studies

       ch = get_studies(accession_id='GCST000854')
       ch
       ch.raw_data
       ch.mode
       ch.studies
       ch.efo_traits
       ch.genotyping_technologies
       ch.discovery_ancestries
       ch.replication_ancestries
       ch.bg_efo_traits
       ch.cohorts
       ch.array_manufacturers
       ```
       Subset object s by either identifier or position
       ```python
       all_df = get_studies()
       all_df[0].studies
       all_df[0:3].studies
       all_df['GCST000854'].studies
       all_df[('GCST000854','GCST90841394','GCST90256308')].studies
       ```
       Objects can be manipulated like sets in the mathematical sense.
       ```python
       one = get_studies(accession_id='GCST000854')
       two = get_studies(accession_id='GCST90841394')
       three = get_studies(accession_id='GCST90256308')
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
            for col in _STUDY_COLUMNS:
                value = normalized.get(col, None)
                standardized[col] = value

            for field in _NESTED_FIELDS:
                standardized[field] = normalized.get(field, None)

            standardized_data.append(standardized)

        self.raw_data = standardized_data

        if len(standardized_data) == 0:
            self.studies = DataFrame(columns=_STUDY_COLUMNS)
            self.efo_traits = DataFrame(columns=_EFO_TRAITS_COLUMNS)
            self.genotyping_technologies = DataFrame(columns=_GENOTYPING_TECHNOLOGIES_COLUMNS)
            self.discovery_ancestries = DataFrame(columns=_DISCOVERY_ANCESTRIES_COLUMNS)
            self.replication_ancestries = DataFrame(columns=_REPLICATION_ANCESTRIES_COLUMNS)
            self.bg_efo_traits = DataFrame(columns=_BG_EFO_TRAITS_COLUMNS)
            self.cohorts = DataFrame(columns=_COHORTS_COLUMNS)
            self.array_manufacturers = DataFrame(columns=_ARRAY_MANUFACTURERS_COLUMNS)
            return

        study_rows = []
        efo_traits_rows = []
        genotyping_technologies_rows = []
        discovery_ancestries_rows = []
        replication_ancestries_rows = []
        bg_efo_traits_rows = []
        cohorts_rows = []
        array_manufacturers_rows = []

        for record in standardized_data:
            accession_id = record.get('accession_id', None)

            # main table: fixed fields, missing fields filled with None
            study_rows.append({
                col: record.get(col, None) for col in _STUDY_COLUMNS
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
                        'accession_id': accession_id,
                        'efo_id': item.get('efo_id', None),
                        'efo_trait': item.get('efo_trait', None)
                    })
                else:
                    efo_traits_rows.append({
                        'accession_id': accession_id,
                        'efo_id': None,
                        'efo_trait': item
                    })

            # genotyping_technologies
            raw_genotyping_technologies = record.get('genotyping_technologies', None)
            if raw_genotyping_technologies is None:
                iterable_genotyping_technologies = []
            elif isinstance(raw_genotyping_technologies, list):
                iterable_genotyping_technologies = raw_genotyping_technologies
            else:
                iterable_genotyping_technologies = [raw_genotyping_technologies]
            for item in iterable_genotyping_technologies:
                genotyping_technologies_rows.append({
                    'accession_id': accession_id,
                    'genotyping_technologies': item
                })

            # discovery_ancestry
            raw_discovery_ancestry = record.get('discovery_ancestry', None)
            if raw_discovery_ancestry is None:
                iterable_discovery_ancestry = []
            elif isinstance(raw_discovery_ancestry, list):
                iterable_discovery_ancestry = raw_discovery_ancestry
            else:
                iterable_discovery_ancestry = [raw_discovery_ancestry]
            for item in iterable_discovery_ancestry:
                discovery_ancestries_rows.append({
                    'accession_id': accession_id,
                    'discovery_ancestry': item
                })

            # replication_ancestry
            raw_replication_ancestry = record.get('replication_ancestry', None)
            if raw_replication_ancestry is None:
                iterable_replication_ancestry = []
            elif isinstance(raw_replication_ancestry, list):
                iterable_replication_ancestry = raw_replication_ancestry
            else:
                iterable_replication_ancestry = [raw_replication_ancestry]
            for item in iterable_replication_ancestry:
                replication_ancestries_rows.append({
                    'accession_id': accession_id,
                    'replication_ancestry': item
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
                        'accession_id': accession_id,
                        'efo_id': item.get('efo_id', None),
                        'efo_trait': item.get('efo_trait', None)
                    })
                else:
                    bg_efo_traits_rows.append({
                        'accession_id': accession_id,
                        'efo_id': None,
                        'efo_trait': item
                    })

            # cohort
            raw_cohort = record.get('cohort', None)
            if raw_cohort is None:
                iterable_cohort = []
            elif isinstance(raw_cohort, list):
                iterable_cohort = raw_cohort
            else:
                iterable_cohort = [raw_cohort]
            for item in iterable_cohort:
                cohorts_rows.append({
                    'accession_id': accession_id,
                    'cohort': item
                })

            # array_manufacturer
            raw_array_manufacturer = record.get('array_manufacturer', None)
            if raw_array_manufacturer is None:
                iterable_array_manufacturer = []
            elif isinstance(raw_array_manufacturer, list):
                iterable_array_manufacturer = raw_array_manufacturer
            else:
                iterable_array_manufacturer = [raw_array_manufacturer]
            for item in iterable_array_manufacturer:
                array_manufacturers_rows.append({
                    'accession_id': accession_id,
                    'array_manufacturer': item
                })

        self.studies = DataFrame(study_rows).reindex(columns=_STUDY_COLUMNS) if len(study_rows) else DataFrame(columns=_STUDY_COLUMNS)
        self.efo_traits = DataFrame(efo_traits_rows).reindex(columns=_EFO_TRAITS_COLUMNS) if len(efo_traits_rows) else DataFrame(columns=_EFO_TRAITS_COLUMNS)
        self.genotyping_technologies = DataFrame(genotyping_technologies_rows).reindex(columns=_GENOTYPING_TECHNOLOGIES_COLUMNS) if len(genotyping_technologies_rows) else DataFrame(columns=_GENOTYPING_TECHNOLOGIES_COLUMNS)
        self.discovery_ancestries = DataFrame(discovery_ancestries_rows).reindex(columns=_DISCOVERY_ANCESTRIES_COLUMNS) if len(discovery_ancestries_rows) else DataFrame(columns=_DISCOVERY_ANCESTRIES_COLUMNS)
        self.replication_ancestries = DataFrame(replication_ancestries_rows).reindex(columns=_REPLICATION_ANCESTRIES_COLUMNS) if len(replication_ancestries_rows) else DataFrame(columns=_REPLICATION_ANCESTRIES_COLUMNS)
        self.bg_efo_traits = DataFrame(bg_efo_traits_rows).reindex(columns=_BG_EFO_TRAITS_COLUMNS) if len(bg_efo_traits_rows) else DataFrame(columns=_BG_EFO_TRAITS_COLUMNS)
        self.cohorts = DataFrame(cohorts_rows).reindex(columns=_COHORTS_COLUMNS) if len(cohorts_rows) else DataFrame(columns=_COHORTS_COLUMNS)
        self.array_manufacturers = DataFrame(array_manufacturers_rows).reindex(columns=_ARRAY_MANUFACTURERS_COLUMNS) if len(array_manufacturers_rows) else DataFrame(columns=_ARRAY_MANUFACTURERS_COLUMNS)

    def __str__(self):
        if self.mode == 'Fat':
            return (
                    "Study is running in fat mode. It has 8 DataFrames with hierarchical dependencies.\n"
                    "studies: %d rows\n"
                    "|\n -efo_traits: %d rows\n"
                    "|\n -genotyping_technologies: %d rows\n"
                    "|\n -discovery_ancestries: %d rows\n"
                    "|\n -replication_ancestries: %d rows\n"
                    "|\n -bg_efo_traits: %d rows\n"
                    "|\n -cohorts: %d rows\n"
                    "|\n -array_manufacturers: %d rows"
                    % (
                        len(self.studies),
                        len(self.efo_traits),
                        len(self.genotyping_technologies),
                        len(self.discovery_ancestries),
                        len(self.replication_ancestries),
                        len(self.bg_efo_traits),
                        len(self.cohorts),
                        len(self.array_manufacturers),
                    )
            )
        if self.mode == 'Thin':
            return (
                    'Study is running in thin mode. It has 1 list that contains the raw data.\n'
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
            raw_data_dict[j['accession_id']] = j

        sub_set = []
        for i in arr:
            if isinstance(i, str):
                sub_set.append(raw_data_dict[i])
            elif isinstance(i, int):
                sub_set.append(raw_data[i])
            else:
                raise TypeError('Invalid item type: {}'.format(type(i)))

        return Study(sub_set, self.mode)

    def __add__(self, other):
        if self.mode == other.mode:
            return Study(self.raw_data + other.raw_data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __sub__(self, other):
        if self.mode == other.mode:
            self_key_set = set()
            self_dict = {}
            other_key_set = set()

            for i in self.raw_data:
                self_key_set.add(i['accession_id'])
                self_dict[i['accession_id']] = i

            for j in other.raw_data:
                other_key_set.add(j['accession_id'])

            sub_key = self_key_set - other_key_set
            data = []
            for k in sub_key:
                data.append(self_dict[k])

            return Study(data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __and__(self, other):
        if self.mode == other.mode:
            self_key_set = set()
            self_dict = {}
            other_key_set = set()

            for i in self.raw_data:
                self_key_set.add(i['accession_id'])
                self_dict[i['accession_id']] = i

            for j in other.raw_data:
                other_key_set.add(j['accession_id'])

            sub_key = self_key_set & other_key_set
            data = []
            for k in sub_key:
                data.append(self_dict[k])

            return Study(data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __or__(self, other):
        if self.mode == other.mode:
            and_dict = {}

            for i in self.raw_data:
                and_dict[i['accession_id']] = i

            for j in other.raw_data:
                and_dict[j['accession_id']] = j

            data = list(and_dict.values())
            return Study(data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __xor__(self, other):
        if self.mode == other.mode:
            self_key_set = set()
            and_dict = {}
            other_key_set = set()

            for i in self.raw_data:
                self_key_set.add(i['accession_id'])
                and_dict[i['accession_id']] = i

            for j in other.raw_data:
                other_key_set.add(j['accession_id'])
                and_dict[j['accession_id']] = j

            sub_key = self_key_set ^ other_key_set
            data = []
            for k in sub_key:
                data.append(and_dict[k])

            return Study(data, self.mode)
        else:
            raise Exception("Please input the same mode")

    def __eq__(self, other):
        if self is None or other is None:
            return self is None and other is None
        return self.raw_data == other.raw_data and self.mode == other.mode

    def __len__(self):
        return len(self.raw_data)
