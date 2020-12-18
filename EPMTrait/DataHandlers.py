from collections import Counter, defaultdict
import io
import gzip
from typing import Dict, List, Tuple, Union

import joblib
import numpy as np
from tqdm import tqdm


def open_file(file_path) -> Union[io.BufferedReader, io.open]:
    if file_path.endswith('.gz'):
        with io.BufferedReader(gzip.open(file_path, 'rb')) as file:
            return file
    else:
        with open(file_path, 'rb') as file:
            return file


def load_methylation_matrix(file_path: str, samples: List[str] = None,
                            verbose=False, total_sites=450000) -> Tuple[np.ndarray, List[str], List[str]]:
    meth_matrix, row_ids = [], []
    header, sample_indices = None, None
    with io.BufferedReader(gzip.open(file_path, 'rb')) as matrix:
        for line in tqdm(matrix, disable=True if not verbose else False, total=total_sites):
            d_line = line.decode('utf-8').strip().split(',')
            if not header:
                header = [sample for sample in d_line]
                if samples:
                    try:
                        sample_indices = [header.index(sample) - 1 for sample in samples]
                    except ValueError as e:
                        print(f'Passed sample not in matrix header')
                        raise e
                    else:
                        header = [header[0]] + [header[1:][index] for index in sample_indices]
                else:
                    sample_indices = [sample for sample in range(len(header) - 1)]
            else:
                meth_matrix.append(np.asarray(d_line[1:], dtype=np.float)[sample_indices])
                row_ids.append(d_line[0].replace('"', ''))
    return np.array(meth_matrix), header, row_ids


def retrieve_sample_methylation(meth_samples: Dict[str, str], n_jobs=1,
                                verbosity=0, consensus_rows=None):
    batches = defaultdict(list)
    for sample, sample_file in meth_samples.items():
        batches[sample_file].append(sample)
    batch_values = joblib.Parallel(n_jobs=n_jobs,
                                   verbose=verbosity)(joblib.delayed(load_methylation_matrix)
                                                      (*[file, samples]) for file, samples in batches.items())
    batch_rows = []
    for batch in batch_values:
        batch_rows.extend(batch[2])
    if not consensus_rows:
        consensus_rows = [row for row, count in Counter(batch_rows).items() if count == len(batch_values)]
    combined_matrix, combined_header = np.zeros((len(consensus_rows), len(meth_samples))), []
    matrix_index = 0
    for matrix, header, rows in batch_values:
        row_indices = [rows.index(row) for row in consensus_rows]
        _, matrix_columns = matrix.shape
        matrix_end = matrix_index + matrix_columns
        combined_matrix[:, matrix_index:matrix_end] = matrix[row_indices, :]
        matrix_index = matrix_end
        if not combined_header:
            combined_header.extend(header)
        else:
            combined_header.extend(header[1:])
    return combined_matrix, combined_header, consensus_rows

