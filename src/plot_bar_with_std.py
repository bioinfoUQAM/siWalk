'''
2024-05-30
Chao-Jung Wu
plot barplots for machine learning model performance with standard deviation
'''
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def run(infile):
    outfig = infile[:-3] + 'bar_with_std.png'
    df = pd.read_csv(infile, sep = '\t', index_col=False)
    valid_matrices = df.columns
    matrices = 'mcc, accuracy, precision, recall, specificity, f1, roc_auc, f1_weighted, tp, tn, fp, fn'.split(', ')
    matrices = [matrix for matrix in matrices if matrix in valid_matrices]
    num_matrices = max(2, len(matrices))
    fig, axes = plt.subplots(num_matrices, 1, figsize=(12, 6 * num_matrices))
    for i, matrix in enumerate(matrices):
        ax = axes[i]
        df.sort_values(by=matrix, ascending=False, inplace=True)
        models = df['#model'].tolist()
        x = np.arange(len(models))
        means = df[matrix]
        stds = df[matrix + '_std']
        bars = ax.bar(x, means, yerr=stds, capsize=5, color='skyblue', edgecolor='black')
        ax.set_xlabel('Model #')
        ax.set_ylabel(matrix, fontsize=22)
        ax.set_title('Model performance with standard deviation: ' + matrix)
        ax.set_xticks(x)
        ax.set_xticklabels(models)
        ax.tick_params(axis='both', labelsize=16)
    plt.tight_layout(pad=3.0)
    plt.subplots_adjust(hspace=0.5)
    plt.savefig(outfig)
    return outfig


if __name__ == "__main__":
  infile = '../output/random_forrest_k_report_GridSearch.tsv'
  run(infile)
  pass