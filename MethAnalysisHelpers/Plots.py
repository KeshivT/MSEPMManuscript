from typing import List

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from MethAnalysisHelpers.ModelHelpers import basic_func, fit_trend, mae, r2


def plot_known_predicted(known: np.ndarray, predicted: np.ndarray,
                         ax: plt.axis, title: str = '', x_label: str = '', y_label: str = '',
                         sample_colors: List[str] = None, func=None,
                         form=(1 / 2)):
    if func:
        # fit trend line]
        popt, pcov, expected = fit_trend(known, predicted, basic_func, form)
        # get r squared
        model_mae = mae(predicted + 1, expected)
        model_r2 = r2(predicted + 1, expected)
        # format plot label
        plot_label = f'$f(x)={popt[0]:.3f}x^{{{str(form)}}} {"+" if popt[2] > 0 else ""}{popt[2]:.3f},' \
                     f' MAE={model_mae:.3f}, R^{2}={model_r2:.3f}$'
        # plot trend line
        sorted_known, sorted_expected = zip(*np.array(sorted(zip(known, expected), key=lambda x: x[0])))
        ax.plot(sorted_known, sorted_expected, 'k--', label=plot_label)
    # scatter plot
    sns.scatterplot(known, predicted, marker='o', alpha=0.8, hue=sample_colors, ax=ax, color='k', s=100)
    ax.set_title(title, fontsize=18)
    ax.set_xlabel(x_label, fontsize=28)
    ax.set_ylabel(y_label, fontsize=28)
    ax.tick_params(axis='both', which='major', labelsize=24)
    ax.legend(fontsize=24, loc='upper left')