import os

from MethAnalysisHelpers.DataHandlers import load_methylation_matrix, retrieve_sample_methylation

test_data_dir = os.getcwd() + '/test_data/'
test_matrix_1 = f'{test_data_dir}test_matrix.txt.gz'
test_matrix_2 = f'{test_data_dir}test_matrix_2.txt.gz'

test_samples = [f'"GSM{count +1 }"' for count in range(10)]

all_samples = {sample: test_matrix_1 for sample in test_samples[:5]}
all_samples.update({sample: test_matrix_2 for sample in test_samples[5:]})

subset_samples = {sample: test_matrix_1 for sample in test_samples[:4]}
subset_samples.update({sample: test_matrix_2 for sample in test_samples[6:]})


matrix_all, header_all, rows_all = load_methylation_matrix(test_matrix_1, test_samples[0:5],
                                                           verbose=True, total_sites=9999)
matrix_1_2, header_1_2, rows_1_2 = load_methylation_matrix(test_matrix_1, test_samples[:2],
                                                           verbose=True, total_sites=9999)
combined_matrix, combined_header, combined_rows = retrieve_sample_methylation(all_samples,
                                                                              n_jobs=2, verbosity=10)
subset_matrix, subset_header, subset_rows = retrieve_sample_methylation(subset_samples,
                                                                        n_jobs=2, verbosity=10)

# define tests


def test_matrix_all_size():
    assert matrix_all.shape == (9999, 5)


def test_matrix_all_samples():
    for x, y in zip(header_all[1:], test_samples):
        assert x == y


def test_matrix_1_2_size():
    assert matrix_1_2.shape == (9999, 2)


def test_matrix_1_2_samples():
    for x, y in zip(header_1_2[1:], test_samples):
        assert x == y


def test_combined_matrix_size():
    assert combined_matrix.shape == (8999, 10)


def test_combined_matrix_samples():
    for sample in combined_header[1:]:
        assert sample in test_samples


def test_subset_matrix_size():
    assert subset_matrix.shape == (8999, 8)


def test_subset_matrix_samples():
    subset_val_samples = test_samples[:4] + test_samples[6:]
    for sample in subset_header[1:]:
        assert sample in subset_val_samples
