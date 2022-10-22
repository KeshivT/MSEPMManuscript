from typing import Dict, List
import numpy as np


class QuantileNormalizeSample:

    def __init__(self, normalization_reference: Dict[str, np.ndarray] = None,
                 probe_ids: Dict[str, set] = None):
        self.ref = normalization_reference
        self.probe_ids = probe_ids

    def normalize_sample(self, row_labels: List[str], row_values: np.ndarray) -> np.ndarray:
        normalized_vals = {}
        for probe_type, probe_ids in self.probe_ids.items():
            probe_vals = [(pid, val) for pid, val in zip(row_labels, row_values) if pid in probe_ids]
            probe_vals.sort(key=lambda x: x[1])
            normalized_vals.update({probe[0]: val for probe, val in zip(probe_vals, self.ref[probe_type])})
        return np.array([normalized_vals[probe] for probe in row_labels])
