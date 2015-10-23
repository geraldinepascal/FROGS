#
# Copyright (C) 2015 INRA
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

__author__ = 'Frederic Escudie - Plateforme bioinformatique Toulouse'
__copyright__ = 'Copyright (C) 2015 INRA'
__license__ = 'GNU General Public License'
__version__ = '0.8.0'
__email__ = 'frogs@toulouse.inra.fr'
__status__ = 'prod'

import sys
import time
import json
import random


class DenseData( list ):
    """
    @summary: Storage for 2D matrix.
    """
    def __init__( self, matrix=None, nb_rows=0, nb_columns=0 ):
        """
        @return: [DenseData] The matrix dense representation
        @param matrix: [list] The list of list.
        @param nb_rows: [int] The number of rows. It is automatically calculated
                        if you set the argument matrix.
        @param nb_columns: [int] The number of rows. It is automatically
                           calculated if you set the argument matrix.
        """
        if matrix is not None:
            super(DenseData, self).__init__(matrix)
            if nb_rows == 0:
                nb_rows = len(self)
            if nb_columns == 0 and nb_rows != 0:
                nb_columns = len(self[0])
        else:
            super(DenseData, self).__init__()
        self.nb_rows = nb_rows
        self.nb_columns = nb_columns

    def _to_json( self ):
        """
        @return: [list] The counts for each column, by row.
        @note: [[0, 0, 0]
                [250, 0, 0], # The column 0 of the row 1 has a count of 250.
                [0, 0, 0],
                [0, 0, 100], # The column 2 of the row 3 has a count of 100.
                [0, 521, 0]] # The column 1 of the row 4 has a count of 521.
        """
        return self

    def get_matrix_type( self ):
        """
        @return: [str] The matrix type.
        """
        return "dense"

    def remove_col( self, remove_idx ):
        """
        @summary: Remove all the count for the column provided.
        @param remove_idx: [int] The real index of the column to remove.
        """
        for row_idx in range(len(self)):
            del self[row_idx][remove_idx]
        self.nb_columns -= 1

    def remove_rows( self, remove_idx ):
        """
        @summary: Remove all the count for the rows provided.
        @param remove_idx: [list] The real index of the rows to remove.
        """
        for current_idx in sorted(remove_idx, reverse=True):
            del self[current_idx]
        self.nb_rows -= len(remove_idx)

    def merge_col( self, sum_idx, added_idx ):
        """
        @summary: Merge two columns. The count of each row of the first column (sum_idx) becomes the sum of the values of the two columns ; the second column is deleted.
        @param sum_idx: [int] The index of the first column to merge. This column is replaced by the new merged column.
        @param added_idx: [int] The index of the second column to merge. This column is deleted after the process.
        """
        for row in self:
            row[sum_idx] += row[added_idx]
        self.remove_col( added_idx )

    def nb_at( self, row_idx, col_idx ):
        """
        @return: [int] The count for the column col_idx in row row_idx.
        @param row_idx: [int] The index of the row.
        @param col_idx: [int] The index of the column.
        """
        return self[row_idx][col_idx]

    def get_col_sum( self, col_idx ):
        """
        @return: [int] The sum of count for the column col_idx.
        @param col_idx : [int] The index of the column.
        """
        total = 0
        for current_row in self:
            total += current_row[col_idx]
        return total

    def get_row_sum( self, row_idx ):
        """
        @return: [int] The sum of count for the row row_idx.
        @param row_idx: [int] The index of the row.
        """
        total = 0
        for current_col in self[row_idx]:
            total += current_col
        return total

    def get_row_idx_by_col( self, col_idx ):
        """
        @return: [list] The list of indexes for the non empty rows in the 
                 specidied column.
        @param col_idx: [int] The index of the column.
        """
        for row_idx, values in enumerate(self):
            if values[col_idx] != 0:
                yield row_idx

    def get_col_idx_by_row( self, row_idx ):
        """
        @return: [list] The list of indexes for the non empty columns in the 
                 specidied row.
        @param row_idx: [int] The index of the column.
        """
        for col_idx, value in enumerate(self[row_idx]):
            if value != 0:
                yield col_idx

    def add( self, row_idx, col_idx, value ):
        """
        @summary: Add the 'value' to the count for the column col_idx in row row_idx.
        @param row_idx: [int] The index of the row.
        @param col_idx: [int] The index of the column.
        @param value: [int] The value to add.
        """
        self[row_idx][col_idx] += value

    def subtract( self, row_idx, col_idx, value ):
        """
        @summary: Subtract the 'value' to the count for the column col_idx in row row_idx.
        @param row_idx: [int] The index of the row.
        @param col_idx: [int] The index of the column.
        @param value: [int] The value to subtract.
        """
        if self[row_idx][col_idx] >= value:
            self[row_idx][col_idx] -= value
        else:
            raise Exception( "'" + str(value) + "' cannot be subtract from row " + str(row_idx) + " column " + str(col_idx) + "." ) 

    def change( self, row_idx, col_idx, value ):
        """
        @summary: Change the 'value' to the count for the column col_idx in row row_idx.
        @param row_idx: [int] The index of the row.
        @param col_idx: [int] The index of the column.
        @param value: [int] The new value.
        """
        self[row_idx][col_idx] = value

    def random_by_col( self, col_idx ):
        """
        @return: [int] The index of the selected row by random sampling in values of specified column.
        @param col_idx: [int] The column where the random sampling is processed.
        @todo test
        """
        elt_index = random.randint(1, self.get_col_sum(col_idx))
        find = False
        row_idx = 0
        previous_elt = 0
        while not find:
            current_nb = previous_elt + self.nb_at( row_idx, col_idx )
            if elt_index <= current_nb:
                find = True
            # Next row
            previous_elt = current_nb
            row_idx += 1
        return( row_idx -1 )

    def random_extract_by_col( self, col_idx, nb_elt ):
        """
        @return: [list] The index(es) of the extracted row(s) by random sampling in values of specified column.
        @param col_idx: [int] The column where the random sampling is processed.
        @param nb_elt: [int] The number of extracted elements.
        @todo test
        """
        selected_rows = dict()
        nb_col_elt = self.get_col_sum(col_idx)
        if nb_elt > nb_col_elt:
            raise Exception("The number of element selected by random sampling is superior than the total number of elements.")
        nb_selected = 0
        while nb_selected < nb_elt:
            elt_index = random.randint(1, nb_col_elt)
            find = False
            row_idx = 0
            previous_elt = 0
            while not find:
                current_nb = previous_elt + self.nb_at( row_idx, col_idx )
                if elt_index <= current_nb:
                    find = True
                # Next row
                previous_elt = current_nb
                row_idx += 1
            self.subtract(row_idx-1, col_idx, 1)
            nb_col_elt -= 1
            nb_selected += 1
            if not selected_rows.has_key(row_idx-1):
                selected_rows[row_idx-1] = 0
            selected_rows[row_idx-1] += 1
        return( selected_rows )

    def get_row_array( self, row_idx ):
        """
        @return: [list] The count for the row for each column.
        @param row_idx: [int] The index of the row.
        @note: In return '[0, 2, 0, 0]' only the second column has an observation.
        """
        return self[row_idx]

    def get_col_array( self, col_idx ):
        """
        @return: [list] The count for the column for each row.
        @param col_idx: [int] The index of the column.
        @note: In return '[0, 2, 0, 0]' only the second row has an observation.
        """
        array = []
        for row in self:
            array.append(row[col_idx])
        return array

    def add_row( self ):
        """
        @summary: Adds an empty row.
        """
        self.append([0] * self.nb_columns)
        self.nb_rows += 1

    def add_col( self ):
        """
        @summary: Adds an empty column.
        """
        for row in self:
            row.append(0)
        self.nb_columns += 1


class SparseData( dict ):
    """
    @summary: Light storage for 2D matrix (0 are not store).
    """
    def __init__( self, matrix=None, nb_rows=0, nb_columns=0):
        """
        @return: [SparseData] The matrix dense representation
        @param matrix: [dict] The dict of dict.
        @param nb_rows: [int] The number of rows.
        @param nb_columns: [int] The number of rows.
        """
        ini_list = matrix if matrix is not None else list()
        for data in ini_list:
            if not self.has_key( data[0] ):
                self[data[0]] = dict()
            self[data[0]][data[1]] = data[2]
        self.nb_rows = nb_rows
        self.nb_columns = nb_columns

    def _to_json( self ):
        """
        @return: [list] The counts ordered by row, by column.
        @example: [[ 1, 0, 250 ]   # The column 0 of the row 1 has a count of 250.
                   [ 8, 2, 100 ]   # The column 2 of the row 8 has a count of 100.
                   [ 9, 1, 521 ]]  # The column 1 of the row 9 has a count of 521.
        """
        sparse = list()
        for rows_idx in sorted(self.keys(), key=int):
            for columns_idx in sorted(self[rows_idx].keys()):
                sparse.append([ rows_idx, columns_idx, self[rows_idx][columns_idx] ])
        return sparse

    def get_matrix_type( self ):
        """
        @return: [str] The matrix type.
        """
        return "sparse"

    def remove_col( self, remove_idx ):
        """
        @summary: Remove all the count for the column provided.
        @param remove_idx: [int] The real index of the column to remove.
        """
        for rows_idx in self.keys():
            # Remove data
            if self[rows_idx].has_key( remove_idx ):
                del self[rows_idx][remove_idx]
            # Change index
            row_columns_idx = sorted( self[rows_idx].keys(), key=int )
            for column_idx in row_columns_idx:
                if column_idx > remove_idx:
                    self[rows_idx][column_idx -1] = self[rows_idx][column_idx]
                    del self[rows_idx][column_idx]
        self.nb_columns -= 1

    def remove_rows( self, remove_idx ):
        """
        @summary: Remove all the count for the rows provided.
        @param remove_idx: [list] The real index of the rows to remove.
        """
        remove_idx.sort(key=int)
        # Remove data
        for idx in remove_idx:
            if self.has_key( idx ):
                del self[idx]
        # Change indexes
        all_rows_idx = sorted( self.keys(), key=int )
        nb_removed = len(remove_idx)
        offset = 0
        next_idx_in_removed = 0
        for row_idx in all_rows_idx:
            while next_idx_in_removed < nb_removed and row_idx > remove_idx[next_idx_in_removed]:
                offset += 1
                next_idx_in_removed += 1
            if offset > 0:
                self[row_idx - offset] = self[row_idx]
                del self[row_idx]
        self.nb_rows -= len(remove_idx)

    def merge_col( self, sum_idx, added_idx ):
        """
        @summary: Merge two columns. The count of each row of the first column (sum_idx) becomes the sum of the values of the two columns ; the second column is deleted.
        @param sum_idx: [int] The index of the first column to merge. This column is replaced by the new merged column.
        @param added_idx: [int] The index of the second column to merge. This column is deleted after the process.
        """
        # Merge counts
        added_values = dict()
        for row_idx in self.keys():
            if self[row_idx].has_key( added_idx ):
                self.add( row_idx, sum_idx, self[row_idx][added_idx] )
        # Remove column
        self.remove_col( added_idx )

    def nb_at( self, row_idx, col_idx ):
        """
        @return: [int] The count for the column col_idx in row row_idx.
        @param row_idx: [int] The index of the row.
        @param col_idx: [int] The index of the column.
        """
        nb = 0
        if self.has_key(row_idx) and self[row_idx].has_key(col_idx):
            nb = self[row_idx][col_idx]
        return nb

    def get_col_sum( self, col_idx ):
        """
        @return: [int] The sum of count for the column col_idx.
        @param col_idx : [int] The index of the column.
        """
        total = 0
        for row_idx in self.keys():
            if self[row_idx].has_key( col_idx ):
                total += self[row_idx][col_idx]
        return total

    def get_row_sum( self, row_idx ):
        """
        @return: [int] The sum of count for the row row_idx.
        @param row_idx: [int] The index of the row.
        """
        total = 0
        if self.has_key( row_idx ):
            for column_idx in self[row_idx].keys():
                total += self[row_idx][column_idx]
        return total

    def get_row_idx_by_col( self, col_idx ):
        """
        @return: [list] The list of indexes for the non empty rows in the 
                 specidied column.
        @param col_idx: [int] The index of the column.
        """
        for row_idx in self.keys():
            if self[row_idx].has_key( col_idx ):
                yield row_idx

    def get_col_idx_by_row( self, row_idx ):
        """
        @return: [list] The list of indexes for the non empty columns in the 
                 specidied row.
        @param row_idx: [int] The index of the column.
        """
        if self.has_key(row_idx):
            for col_idx in self[row_idx]:
                yield col_idx

    def get_row_array( self, row_idx ):
        """
        @return: [list] The count for the row for each column.
        @param row_idx: [int] The index of the row.
        @note: In return '[0, 2, 0, 0]' only the second column has an observation.
        """
        array = [0 for current in range(self.nb_columns)]
        if self.has_key( row_idx ):
            for column_idx in sorted( self[row_idx].keys() ):
                array[column_idx] = self[row_idx][column_idx]
        return array

    def get_col_array( self, col_idx ):
        """
        @return: [list] The count for the column for each row.
        @param col_idx: [int] The index of the column.
        @note: In return '[0, 2, 0, 0]' only the second row has an observation.
        """
        array = [0 for current in range(self.nb_rows)]
        for row_idx in self:
            if self[row_idx].has_key(col_idx):
                array[row_idx] = self[row_idx][col_idx]
        return array

    def add( self, row_idx, col_idx, value ):
        """
        @summary: Add the 'value' to the count for the column col_idx in row row_idx.
        @param row_idx : [int] The index of the row.
        @param col_idx : [int] The index of the column.
        @param value : [int] The value to add.
        """
        if not self.has_key( row_idx ):
            self[row_idx] = { col_idx : 0 }
        elif not self[row_idx].has_key( col_idx ):
            self[row_idx][col_idx] = 0
        self[row_idx][col_idx] += value

    def subtract( self, row_idx, col_idx, value ):
        """
        @summary: Subtract the 'value' to the count for the column col_idx in row row_idx.
        @param row_idx: [int] The index of the row.
        @param col_idx: [int] The index of the column.
        @param value: [int] The value to subtract.
        """
        if self.has_key( row_idx ) and self[row_idx].has_key( col_idx ) and self[row_idx][col_idx] >= value:
            self[row_idx][col_idx] -= value
        else:
            raise Exception( "'" + str(value) + "' cannot be subtract from row " + str(row_idx) + " column " + str(col_idx) + "." ) 

    def change( self, row_idx, col_idx, value ):
        """
        @summary: Change the 'value' to the count for the column col_idx in row row_idx.
        @param row_idx: [int] The index of the row.
        @param col_idx: [int] The index of the column.
        @param value: [int] The new value.
        """
        if value != 0:
            if not self.has_key( row_idx ):
                self[row_idx] = { col_idx : value }
            else:
                self[row_idx][col_idx] = value
        else:
            if self.has_key( row_idx ) and self[row_idx].has_key( col_idx ) :
                del self[row_idx][col_idx]

    def random_by_col( self, col_idx ):
        """
        @return: [int] The index of the selected row by random sampling in values of specified column.
        @param col_idx: [int] The column where the random sampling is processed.
        """
        elt_index = random.randint(1, self.get_col_sum(col_idx))
        find = False
        row_idx = 0
        previous_elt = 0
        while not find:
            current_nb = previous_elt + self.nb_at( row_idx, col_idx )
            if elt_index <= current_nb:
                find = True
            # Next row
            previous_elt = current_nb
            row_idx += 1
        return( row_idx -1 )

    def random_extract_by_col( self, col_idx, nb_elt ):
        """
        @return: [list] The index(es) of the extracted row(s) by random sampling in values of specified column.
        @param col_idx: [int] The column where the random sampling is processed.
        @param nb_elt: [int] The number of extracted elements.
        """
        selected_rows = dict()
        nb_col_elt = self.get_col_sum(col_idx)
        if nb_elt > nb_col_elt:
            raise Exception("The number of element selected by random sampling is superior than the total number of elements.")
        nb_selected = 0
        while nb_selected < nb_elt:
            elt_index = random.randint(1, nb_col_elt)
            find = False
            row_idx = 0
            previous_elt = 0
            while not find:
                current_nb = previous_elt + self.nb_at( row_idx, col_idx )
                if elt_index <= current_nb:
                    find = True
                # Next row
                previous_elt = current_nb
                row_idx += 1
            self.subtract(row_idx-1, col_idx, 1)
            nb_col_elt -= 1
            nb_selected += 1
            if not selected_rows.has_key(row_idx-1):
                selected_rows[row_idx-1] = 0
            selected_rows[row_idx-1] += 1
        return( selected_rows )

    def add_row( self ):
        """
        @summary: Adds an empty row.
        """
        self.nb_rows += 1

    def add_col( self ):
        """
        @summary: Adds an empty column.
        """
        self.nb_columns += 1


class Biom:
    """
    @summary: Store biological sample by observation contingency tables.
    @see: https://github.com/biom-format
    """
    def __init__( self, id=None, format="Biological Observation Matrix 1.0.0", 
                  format_url="http://biom-format.org/documentation/format_versions/biom-1.0.html",
                  type="OTU table", generated_by=None, date=None, rows=None,
                  columns=None, matrix_type="dense", matrix_element_type="int",
                  data=None ):
        """
        @param id: [int]
        @param format: [str]
        @param format_url: [str]
        @param type: [str]
        @param generated_by: [str]
        @param date: [str]
        @param rows: [list]
        @param columns: [list]
        @param matrix_type: [str]
        @param matrix_element_type: [str]
        @param data: [list]
        """
        self.id = id
        self.format = format
        self.format_url = format_url
        self.type = type
        self.generated_by = generated_by
        self.date = date if date is not None else time.strftime('%Y-%m-%dT%H:%M:%S',time.localtime())
        self.rows = rows if rows is not None else list()
        self.columns = columns if columns is not None else list()
        self.matrix_element_type = matrix_element_type
        ini_data = data if data is not None else list()
        if matrix_type == "dense":
            self.data = DenseData( ini_data )
        else:
            self.data = SparseData( ini_data, len(self.rows), len(self.columns) )
        self._obs_index = {} if rows is None else {key['id']:idx for idx, key in enumerate(rows)} # observation index in self.rows (increase fin_idx speed)

    def __str__(self):
        return str( self.__dict__ )

    def __repr__(self):
        return str( self.__dict__ )

    def remove_observations( self, observations_names ):
        """
        @summary: Removes the specified observations.
        @param observations_names : [list] The IDs of the observations to remove.
        """
        removed_observations_idx = list()
        # Find indexes
        for current_observation in observations_names:
            observation_idx = self.find_idx( "observation", current_observation )
            removed_observations_idx.append( observation_idx )
        removed_observations_idx.sort()
        # Remove OTU from the _obs_idx
        nb_removed = len(observations_names)
        offset = 0
        next_idx_in_removed = 0
        for current_idx, observation_name in enumerate(self.get_observations_names()):
            if next_idx_in_removed < nb_removed and current_idx == removed_observations_idx[next_idx_in_removed]: # current observation must be deleted
                next_idx_in_removed += 1
                offset += 1
                del self._obs_index[observation_name]
            elif offset > 0:
                self._obs_index[observation_name] -= offset
        # Remove OTU from the self.rows
        for observation_idx in removed_observations_idx[::-1]:
            del self.rows[observation_idx]
        # Remove OTU from the self.data
        self.data.remove_rows( removed_observations_idx )

    def reset_count_by_replicates_evidence( self, samples_names, min_evidence_nb=2 ):
        """
        @summary: Puts to 0 the counts of an observation for all samples in a
                  replicate if the number of samples with this observation is
                  lower than 'min_evidence_nb' in the replicate.
        @param samples_names: [list] The names of the replicates.
        @param min_evidence_nb: [int] The minimun number of replicates with the
                                observation.
        @note: Example with min_evidence_nb = 2 and 2 triplicates (A and B)
               Before process
                     spl_A1 spl_A2 spl_A3 spl_B1 spl_B2 spl_B3
                 obs    1      1      0      0      0      2
               After process
                     spl_A1 spl_A2 spl_A3 spl_B1 spl_B2 spl_B3
                 obs    1      1      0      0      0      0
        """
        samples_idx = [self.find_idx("sample", sample) for sample in samples_names]
        # For each observation
        for row_idx in range( len(self.rows) ):
            obs_evidence = 0
            # Process evidence number
            for col_idx in samples_idx:
                if self.data.nb_at(row_idx, col_idx) > 0: # if count > 0
                    obs_evidence += 1
            # If the evidence is insufficient
            if obs_evidence < min_evidence_nb and obs_evidence > 0:
                # for each sample
                for col_idx in samples_idx:
                    # Set observation count to 0
                    self.data.change( row_idx, col_idx, 0 )

    def filter_observations_by_count( self, min_nb, max_nb=None ):
        """
        @summary: Filter observations on count value.
        @param min_nb: [int] The observations with a count inferior to 'min_nb' is removed.
        @param max_nb: [int] The observations with a count superior to 'min_nb' is removed.
        """
        removed_obs_names = list()
        for observation_idx in range( len(self.rows) ):
            observation_count = self.data.get_row_sum( observation_idx )
            if observation_count < min_nb or (max_nb is not None and observation_count > max_nb):
                removed_obs_names.append( self.rows[observation_idx]["id"] )
        self.remove_observations( removed_obs_names )

    def merge_samples( self, samples, merged_sample_id=None ):
        """
        @summary: Merge data and metadata of a list of samples.
        @param samples: [list] Samples to merge.
        @param merged_sample_id: [str] Name for the new meta-sample.
        """
        # Preprocess final_sample
        final_idx = self.find_idx( "sample", samples[0] )
        final_sample = self.columns[final_idx]
        if final_sample['metadata'] is not None:
            metadata_names = final_sample['metadata'].keys()
            for metadata_name in metadata_names:
                final_sample['metadata'][final_sample['id'] + ":" + metadata_name] = final_sample['metadata'][metadata_name]
                del final_sample['metadata'][metadata_name]
        else:
            final_sample['metadata'] = dict()
        # Merge
        for current_name in samples[1:]:
            final_idx = self.find_idx( "sample", samples[0] )
            # Find sample
            current_idx = self.find_idx( "sample", current_name )
            current_sample = self.columns[current_idx]
            # Update metadata history
            if final_sample['metadata'].has_key( "merge_history" ):
                final_sample['metadata']['merge_history'] += " AND " + current_sample['id']
            else:
                final_sample['metadata']['merge_history'] = final_sample['id'] + " AND " + current_sample['id']
            # Merge metadata
            if current_sample['metadata'] is not None:
                for metadata_name in current_sample['metadata']:
                    final_sample['metadata'][current_sample['id'] + ":" + metadata_name] = current_sample['metadata'][metadata_name]
            # Merge data
            self.data.merge_col( final_idx, current_idx )
            # Remove sample from the self.columns
            del self.columns[current_idx]
        # If rename final sample
        if merged_sample_id is not None:
            final_sample['id'] = merged_sample_id

    def find_idx( self, subject_type, query_name ):
        """
        @summary: Returns the index of the query.
        @param subject_type: [str] The type of subject : "sample" or "observation".
        @param query_name: [str] The id of the element (ex : "OTU_0012").
        @return: [int] The index of the element.
        """
        find_idx = None
        if subject_type == "observation":
            if self._obs_index.has_key(query_name):
                find_idx = self._obs_index[query_name]
        else:
            idx = 0
            while idx < len(self.columns) and find_idx is None:
                if self.columns[idx]['id'] == query_name:
                    find_idx = idx
                idx += 1
        if find_idx is None:
            raise ValueError( "The " + subject_type + " '" + query_name + "' doesn't exist." )
        return find_idx

    def add_metadata( self, subject_name, metadata_name, metadata_value, subject_type="sample"):
        """
        @summary: Add a metadata on subject (a sample or an observation).
        @param subject_name : [str] Metadata is added to the sample/observation with this name. 
        @param metadata_name : [str] The metadata category (ex : 'taxonomy').
        @param metadata_name : [str] The value of metadata (ex : 'Bacteria').
        @param subject_type : [str] The type of subject : "sample" or "observation".
        """
        # Select subject container
        if subject_type == "sample":
            subject_list = self.columns
        elif subject_type == "observation":
            subject_list = self.rows
        else:
            raise ValueError( "'" + subject_type + "' is an invalid subject type for metadata. Metadata must be add to 'observation' or 'sample'." )
        # Find subject
        try:
            subject_idx = self.find_idx( subject_type, subject_name )
        # Subject does not exist
        except ValueError:
            sys.stderr.write("[WARNING] The metadata named '" + metadata_name + "' can't be added to sample '" + subject_name + "' because it does not exist.\n")
        # Subject exists
        else:
            if subject_list[subject_idx]['metadata'] is None:
                subject_list[subject_idx]['metadata'] = dict()
            elif subject_list[subject_idx]['metadata'].has_key( metadata_name ):
                sys.stderr.write("[WARNING] You erase previous value of the metadata named '" + metadata_name + "' in " + subject_name + " (OLD:'" + str(subject_list[subject_idx]['metadata'][metadata_name]) + "' => NEW:'" + str(metadata_value) + "').\n")
            subject_list[subject_idx]['metadata'][metadata_name] = metadata_value

    def to_json( self ):
        """
        @summary: Return a json format for the data store in the Biom object.
        @return: [str] The json.
        """
        self.shape = [
                       len(self.rows),
                       len(self.columns)
        ]
        self.matrix_type = self.data.get_matrix_type()
        save_data = self.data
        self.data = save_data._to_json()
        save_index = self._obs_index
        del self._obs_index
        json_str = json.dumps( self, default=lambda o: o.__dict__, sort_keys=False, indent=4 )
        self.data = save_data
        self._obs_index = save_index
        del self.shape
        del self.matrix_type
        return json_str

    def remove_samples( self, samples_names ):
        """
        @summary: Removes sample(s) from biom.
        @param samples_names: [str] The name of the sample to rename.
        """
        for current_sample in samples_names :
            sample_idx = self.find_idx( "sample", current_sample )
            # Remove sample from the self.columns
            del self.columns[sample_idx]
            # Remove sample from the self.data
            self.data.remove_col( sample_idx )

    def change_count( self, observation_name, sample_name, value ):
        """
        @summary: Replace the count for one observation of one sample by value.
        @param observation_name: [str] The observation name.
        @param sample_name: [str] The sample name.
        @param value: [int] The value to use.
        """
        row_idx = self.find_idx( "observation", observation_name )
        col_idx = self.find_idx( "sample", sample_name )
        self.data.change( row_idx, col_idx, value )

    def subtract_count( self, observation_name, sample_name, value ):
        """
        @summary: Subtract a value to the count for one observation of one sample.
        @param observation_name: [str] The observation name.
        @param sample_name: [str] The sample name.
        @param value: [int] The value to subtract.
        """
        row_idx = self.find_idx( "observation", observation_name )
        col_idx = self.find_idx( "sample", sample_name )
        self.data.subtract( row_idx, col_idx, value )

    def add_count( self, observation_name, sample_name, value ):
        """
        @summary: Add a value to the count for one observation of one sample.
        @param observation_name: [str] The observation name.
        @param sample_name: [str] The sample name.
        @param value: [int] The value to add.
        """
        row_idx = self.find_idx( "observation", observation_name )
        col_idx = self.find_idx( "sample", sample_name )
        self.data.add( row_idx, col_idx, value )

    def get_count(self, observation_name, sample_name):
        """
        @summary : Returns the count for one observation of one sample.
         @param observation_name : [str] The observation name.
         @param sample_name : [str] The sample name.
        """
        row_idx = self.find_idx( "observation", observation_name )
        col_idx = self.find_idx( "sample", sample_name )
        return self.data.nb_at( row_idx, col_idx )

    def add_observation( self, observation_name, metadata=None ):
        """
        @summary : Add one observation in biom.
         @param observation_name : [str] The observation name.
         @param metadata : [dict] The metadata (keys : metadata names ; values : metadata values).
        """
        ini_metadata = metadata if metadata is not None else dict()
        try:
            self.find_idx( "observation", observation_name )
        # Observation doesn't exist
        except ValueError:
            self.rows.append( {'id':observation_name, 'metadata':None } )
            self._obs_index[observation_name] = len(self.rows) - 1
            self.data.add_row()
            for metadata_name in ini_metadata.keys():
                self.add_metadata( observation_name, metadata_name, ini_metadata[metadata_name], "observation" )
        # Observation already exists
        else:
            raise ValueError( "The observation '" + observation_name + "' already exists." )

    def add_sample( self, sample_name, metadata=None ):
        """
        @summary : Add one sample in biom.
         @param sample_name : [str] The sample name.
         @param metadata : [dict] The metadata (keys : metadata names ; values : metadata values).
        """
        ini_metadata = metadata if metadata is not None else dict()
        try:
            self.find_idx( "sample", sample_name )
        # Sample doesn't exist
        except ValueError:
            self.columns.append( {'id':sample_name, 'metadata':None } )
            self.data.add_col()
            for metadata_name in ini_metadata.keys():
                self.add_metadata( sample_name, metadata_name, ini_metadata[metadata_name], "sample" )
        # Sample already exists
        else:
            raise ValueError( "The sample '" + sample_name + "' already exists." )

    def get_samples_names( self ):
        """
        @summary : Returns a generator to iterate on samples names.
        @return : [generator] the generator to iterate on samples names.
        """
        for col in self.columns:
            yield col["id"]

    def get_sample_metadata( self, sample_name ):
        """
        @summary : Returns the sample metadata.
        @return : [dict] the sample metadata.
        """
        return self.columns[self.find_idx("sample", sample_name)]["metadata"]

    def get_observations_names( self ):
        """
        @summary : Returns a generator to iterate on observations names.
        @return : [generator] the generator to iterate on observations names.
        """
        for col in self.rows:
            yield col["id"]

    def get_observation_metadata( self, observation_name ):
        """
        @summary : Returns the observation metadata.
        @return : [dict] the observation metadata.
        """
        return self.rows[self.find_idx("observation", observation_name)]["metadata"]

    def has_observation_metadata(self, metadata_title):
        has_metadata = False
        for current_observ in self.get_observations():
            if current_observ["metadata"].has_key(metadata_title):
                has_metadata = True
        return has_metadata

    # Add by Maria
    def get_observations( self ):
        """
        @summary: Returns a generator to iterate on all observations.
        @return: [generator] the generator to iterate on observations.
        """
        for observation in self.rows:
            yield observation

    def get_observations_by_sample( self, sample_name ):
        """
        @summary: Returns a generator to iterate on observations present in specified sample.
        @param sample_name: The specified sample.
        @return: [generator] the generator to iterate on observations.
        """
        sample_idx = self.find_idx( "sample", sample_name )
        for observation_idx in self.data.get_row_idx_by_col( sample_idx ):
            yield self.rows[observation_idx]

    def get_samples_by_observation( self, observation_name ):
        """
        @summary: Returns a generator to iterate on samples that contain the specified observation.
        @param observation_name: The specified observation.
        @return: [generator] the generator to iterate on samples.
        """
        observation_idx = self.find_idx( "observation", observation_name )
        for sample_idx in self.data.get_col_idx_by_row( observation_idx ):
            yield self.columns[sample_idx]

    def random_obs_by_sample( self, sample_name ):
        sample_idx = self.find_idx( "sample", sample_name )
        return self.rows[self.data.random_by_col(sample_idx)]

    def random_obs_extract_by_sample( self, sample_name, nb_selected ):
        """
        @summary: Extract and returns random sampled observations.
        @param sample_name: The sample where the observation are sampled.
        @param nb_selected: The number of sampling round(s).
        @return: [list] The list of extracted observations with the number of extraction for each.
        """
        sample_idx = self.find_idx( "sample", sample_name )
        selected_rows = self.data.random_extract_by_col(sample_idx, nb_selected)
        selected_obs = list()
        for idx in selected_rows.keys():
            selected_obs.append({ 'observation':self.rows[idx], 'count': selected_rows[idx] })
        return selected_obs

    def get_sample_count( self, sample_name ):
        return self.data.get_col_sum( self.find_idx("sample", sample_name) )

    def get_sample_obs( self, sample_name ):
        """
        @summary: for sample sample_name return observations count
        @return: sample_name otu list
        """
        return self.data.get_col_array( self.find_idx("sample", sample_name) )

    def get_total_count( self ):
        total_count = 0
        for observation_name in self.get_observations_names():
            total_count += self.get_observation_count(observation_name)
        return total_count

    def get_observation_count( self, observation_name ):
        return self.data.get_row_sum( self.find_idx("observation", observation_name) )

    def get_observations_counts(self):
        """
        @summary : Return the list of the observations counts.
        @return : [list] the observation ID and the observation count for each observation.
                  Example :
                  [
                    ["OTU_1", 128],
                    ["OTU_2", 8]
                  ]
        """
        for observation_idx in range(len(self.rows)):
            yield self.rows[observation_idx]['id'], self.data.get_row_sum(observation_idx)

    def to_count( self ):
        """
        @summary : Returns the count of observations by sample.
        @return : [generator] The generator to iterate on observations. Each observation is a list of count 
                  by sample.
                  Example : [1, 0] # Iteration 1 : sample_1 has one observation_1, sample_2 has zero observation_1
                            [1, 8] # Iteration 2 : sample_1 has one observation_2, sample_2 has eight observation_2
        """
        nb_rows = len(self.rows)
        for row_idx in range(nb_rows):
            yield self.data.get_row_array( row_idx )

    def to_count_table( self ):
        """
        @summary : Returns the count of observations by sample with titles.
        @return : [generator] The generator to iterate on observations. First line is a title.
                  Example : ['#Observation', 'Sample1', 'Sample2'] # Iteration 1 : title
                            ['GG_OTU_1', 1, 0] # Iteration 2 : Sample1 has one GG_OTU_1, Sample1 has zero GG_OTU_1
                            ['GG_OTU_2', 1, 8] # Iteration 3 : Sample2 has one GG_OTU_2, Sample2 has eight GG_OTU_2
        """
        # Return Title
        yield ["#Observation"] + [col['id'] for col in self.columns]
        # Return lines
        row_idx = 0
        for row in self.to_count():
            OTU_name = self.rows[row_idx]['id']
            row_idx += 1
            yield [OTU_name] + row


class BiomIO:
    """
    Reader/Writer for the Biom format.
    The BIOM file format is a json format designed to be a general-use format for representing biological sample by observation contingency tables.
    BIOM is a recognized standard for the Earth Microbiome Project and is a Genomics Standards Consortium candidate project.
    @see: https://github.com/biom-format
    @note: Usage example
           # Remove the samples splA and splB from the observation matrix
           if BiomIO.is_BIOM("/path/to/input/biom"):
               biom = BiomIO.from_json("/path/to/input/biom")
               biom.remove_samples(["splA", "splB"])
               BiomIO.write("/path/to/output/biom", biom)
           
           # Load metadata in BIOM
           biom = BiomIO.from_json("/path/to/input/biom")
           BiomIO.load_metadata(biom, "/path/to/samples/metadata")
           BiomIO.load_metadata(biom, "/path/to/observations/metadata")
           BiomIO.write("/path/to/output/biom", biom)
    """
    @staticmethod
    def from_count_table( count_file, generated_by=None ):
        """
        @summary: Return an object 'Biom' from a count table.
        @param count_file: [str] The path of the count file.
                             Format :
                              #Cluster_ID<TAB>sample1<TAB>sample2
                              OTU1<TAB>8<TAB>10
                              ...
        @param generated_by: [str] The method/software used to generate data.
        @return [Biom] The Biom object.
        """
        biom = Biom()
        biom.data = SparseData()
        biom.generated_by = generated_by

        count_fh = open( count_file )
        row_idx = 0
        for line in count_fh:
            line = line.strip()
            line_fields = line.split()
            # Title line
            if line.startswith('#'):
                for sample in line_fields[1:]:
                    # Load sample (biom.columns)
                    biom.add_sample( sample, None )
            # OTU line
            else:
                # Load OTU (biom.rows)
                biom.add_observation( line_fields[0], None )
                # Load count (biom.data)
                col_idx = 0
                for count in line_fields[1:]:
                    if int(count) != 0:
                        biom.data.add( row_idx, col_idx, int(count) )
                    col_idx += 1
                row_idx += 1
        count_fh.close()

        return biom

    @staticmethod
    def from_json( path ):
        """
        @summary: Return an object 'Biom' from a biom file.
        @param path: [str] The path of the biom file.
        @return: [Biom] The Biom object.
        """
        json_data = open( path )
        python_dict = json.load( json_data )
        json_data.close()

        return Biom( python_dict["id"],
                     python_dict["format"],
                     python_dict["format_url"],
                     python_dict["type"],
                     python_dict["generated_by"],
                     python_dict["date"],
                     python_dict["rows"],
                     python_dict["columns"],
                     python_dict["matrix_type"],
                     python_dict["matrix_element_type"],
                     python_dict["data"] )

    @staticmethod
    def is_BIOM( path ):
        """
        @summary: Return true if the file is a BIOM file.
        @param path: [str] The path of checked file.
        @return: [bool] True if the file is a BIOM file.
        @TODO: test
        """
        is_biom = None
        try:
            BiomIO.from_json(path)
            is_biom = True
        except:
            is_biom = False
        return is_biom

    @staticmethod
    def write( path, biom ):
        """
        @summary: Write a biom file from a 'Biom'.
        @param path: [str] The path of the biom file.
        @param biom: [Biom] The Biom object to write.
        """
        out_fh = open( path, "w" )
        out_fh.write( biom.to_json() )
        out_fh.close()

    @staticmethod
    def write_count_table( path, biom ):
        """
        @summary: Write count table from an object 'Biom'.
        @param path: [str] The path of the biom file.
        @param biom: [Biom] The Biom object to write.
        """
        out_fh = open( path, "w" )
        for line in biom.to_count_table():
            out_fh.write( "\t".join(map(str, line)) + "\n" )
        out_fh.close()

    @staticmethod
    def write_krona_table( path, biom, taxonomy_key="taxonomy" ):
        """
        @todo test
        """
        out_fh = open( path, "w" )
        for idx in range(len(biom.rows)):
            count = biom.data.get_row_sum( idx )
            tax = biom.rows[idx]["metadata"][taxonomy_key]
            if isinstance(tax, list) or isinstance(tax, tuple):
                tax = "\t".join( map(str, tax) )
            else:
                tax = str( tax )
            tax = "\t".join( map(str.strip, tax.split(";")) ) # Replace space separator between ranks by tabulation
            out_fh.write( str(count) + "\t" + tax + "\n" )
        out_fh.close()

    @staticmethod
    def load_metadata( biom, metadata_file, subject_type="sample", types=None, list_sep=None ):
        """
        @summary : Add to biom several metadata from metadata file.
         @param biom : [Biom] The Biom object to update.
         @param metadata_file : [str] The path of the metadata file.
                                Format :
                                #TITLE<TAB>Metadata_1_name<TAB>Metadata_2_name
                                Subject_name<TAB>Metadata_1_value<TAB>Metadata_2_value
                                ...
         @param subject_type : [str] The type of subject : "sample" or "observation".
         @param types : [dict] Types for of the metadata values ("str", "int", "float").
                        Example :
                        {
                          'confidence' : 'float',
                          'rank'       : 'int'
                        }
         @param list_sep : [dict] Separator if the metadata is a list.
                        Example :
                        {
                          'taxonomy'      : ';', # Bacteria;Proteobacteria
                          'environnement' : '/'  # Sea/Ocean
                        }
        """
        ini_types = types if types is not None else dict()
        ini_list_sep = list_sep if list_sep is not None else dict()
        metadata_fh = open( metadata_file )
        metadata = list()
        # Names and type of metadata
        title_line = metadata_fh.readline().strip()
        title_fields = title_line.split()
        for metadata_name in title_fields[1:]:
            metadata_type = "str"
            if ini_types.has_key( metadata_name ):
                metadata_type = ini_types[metadata_name]
            metadata_list_sep = None
            if ini_list_sep.has_key( metadata_name ):
                metadata_list_sep = ini_list_sep[metadata_name]
            metadata.append( {
                              'name'     : metadata_name,
                              'type'     : metadata_type,
                              'list_sep' : metadata_list_sep
            })
        # Values of metadata
        for line in metadata_fh:
            line = line.strip()
            if not line.startswith('#'):
                line_fields = line.split()
                metadata_subject = line_fields[0]
                title_idx = 0
                for metadata_value in line_fields[1:]:
                    # Manage cast metadata value
                    if metadata[title_idx]['type'] == "str":
                        cast = str
                    elif metadata[title_idx]['type'] == "int":
                        cast = int
                    elif metadata[title_idx]['type'] == "float":
                        cast = float
                    else:
                        raise ValueError( "'" + metadata[title_idx]['type'] + "' is an invalid type for metadata. Metadata must be 'str' or 'int' or 'float'." )
                    # Manage split
                    if metadata[title_idx]['list_sep'] is None:
                        metadata_value = cast( metadata_value )
                    else:
                        metadata_value = [cast(value) for value in metadata_value.split(metadata[title_idx]['list_sep'])]
                    # Add metadata
                    biom.add_metadata( metadata_subject, metadata[title_idx]['name'], metadata_value, subject_type)
                    # Next medata title
                    title_idx += 1