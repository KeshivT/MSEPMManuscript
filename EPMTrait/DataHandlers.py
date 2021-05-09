from collections import Counter, defaultdict
import io
import gzip
import os
from typing import Dict, List, Tuple, Union

import joblib
import numpy as np
from tqdm import tqdm


def get_normalization_ref(epic_ref=False):
    ref_dir = f'{os.path.dirname(os.path.realpath(__file__))}/NormalizationRefs/'
    ref_file = f'{ref_dir}450k_ref.txt.gz' if not epic_ref else f'{ref_dir}epic_ref.txt.gz'
    normalization_ref = {'type_1': set(), 'type_2': set()}
    with io.BufferedReader(gzip.open(ref_file, 'rb')) as file:
        for line in file:
            probe, probe_type = line.decode().strip().split('\t')
            normalization_ref[probe_type].add(probe)
    return normalization_ref


def open_file(file_path) -> Union[io.BufferedReader, io.open]:
    if file_path.endswith('.gz'):
        with io.BufferedReader(gzip.open(file_path, 'rb')) as file:
            return file
    else:
        with open(file_path, 'rb') as file:
            return file


def load_methylation_matrix(file_path: str, samples: List[str] = None,
                            verbose=False, total_sites=450000, row_labels: list = None,
                            quantile_norm=False) -> Tuple[np.ndarray, List[str], List[str]]:
    meth_matrix, row_ids = [], []
    header, sample_indices = None, None
    r_labels = None if not row_labels else set(row_labels)
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
                row_id = d_line[0].replace('"', '')
                if r_labels and not quantile_norm:
                    if row_id not in r_labels:
                        continue
                meth_matrix.append(np.asarray(d_line[1:], dtype=np.float)[sample_indices])
                row_ids.append(row_id)
    meth_matrix = np.array(meth_matrix)
    if quantile_norm:
        for sample in range(meth_matrix.shape[1]):
            norm_sample = quantile_norm.normalize_sample(row_ids, meth_matrix[:, sample])
            meth_matrix[:, sample] = norm_sample
        if r_labels:
            row_sites, row_ids = (zip(*[(count, label) for count, label in enumerate(row_ids) if label in r_labels]))
            meth_matrix = meth_matrix[row_sites, :]
    return meth_matrix, header, row_ids


def retrieve_sample_methylation(meth_samples: Dict[str, str], n_jobs=1,
                                verbosity=0, rows=None, quantile_norm=None):
    batches = defaultdict(list)
    for sample, sample_file in meth_samples.items():
        batches[sample_file].append(sample)
    batch_values = joblib.Parallel(n_jobs=n_jobs,
                                   verbose=verbosity)(joblib.delayed(load_methylation_matrix)
                                                      (*[file, samples, False, 1, rows, quantile_norm]) for
                                                      file, samples in batches.items())
    batch_rows = []
    for batch in batch_values:
        batch_rows.extend(batch[2])
    consensus_rows = [row for row, count in Counter(batch_rows).items() if count == len(batch_values)]
    combined_matrix, combined_header = np.zeros((len(consensus_rows), len(meth_samples))), []
    matrix_index = 0
    for matrix, header, rows in batch_values:
        row_counts = {row: count for count, row in enumerate(rows)}
        row_indices = [row_counts[row] for row in consensus_rows]
        _, matrix_columns = matrix.shape
        matrix_end = matrix_index + matrix_columns
        combined_matrix[:, matrix_index:matrix_end] = matrix[row_indices, :]
        matrix_index = matrix_end
        if not combined_header:
            combined_header.extend(header)
        else:
            combined_header.extend(header[1:])
    return combined_matrix, combined_header, consensus_rows

