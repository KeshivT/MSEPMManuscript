import numpy as np


class GenerateForm:
    
    def __init__(self, form=1 / 2):
        self.form = float(form)

    def functional_form(self, x, a, b, c) -> np.ndarray:
        return a * np.asarray(x) ** self.form + c
