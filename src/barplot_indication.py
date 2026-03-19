import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

bar_width = 0.35

def draw_6candidates_interface(infile1, infile2, sequence, ground_truths=None):
    """Draw indication bar chart with 6 candidate predictions.

    Parameters
    ----------
    infile1 : str  - full indication TSV (all positions)
    infile2 : str  - top-6 recommendation TSV
    sequence : str - precursor sequence (DNA)
    ground_truths : list of (label, start, end) or None
        e.g. [('GT1 1st tasiARF', 0, 20), ('GT2 2nd tasiARF', 21, 41)]
    """
    with open(infile2, 'r') as fh:
        fh.readline()  # header
        c,_,_,_,_,d,_ = fh.readline().split('\t')
        e,_,_,_,_,f,_ = fh.readline().split('\t')
        g,_,_,_,_,h,_ = fh.readline().split('\t')
        i,_,_,_,_,j,_ = fh.readline().split('\t')
        k,_,_,_,_,l,_ = fh.readline().split('\t')
        m,_,_,_,_,n,_ = fh.readline().split('\t')
    [c,d,e,f,g,h,i,j,k,l,m,n] = [int(x) for x in [c,d,e,f,g,h,i,j,k,l,m,n]]

    df = pd.read_csv(infile1, sep='\t', index_col=False)
    df.sort_values(by="Position p", ascending=True, inplace=True)
    p, Start, End = df['Position p'], df['Start S(p)'], df['End E(p)']
    bar_width = 0.2
    positions_start = p
    positions_end = positions_start + bar_width

    plt.figure(figsize=(20, 8))
    plt.bar(positions_start, Start, width=bar_width, label="Start S(p)")
    plt.bar(positions_end, End, width=bar_width, label="End E(p)", color='orange')

    # Determine y range for annotation lines
    all_vals = list(Start) + list(End)
    y_min = min(v for v in all_vals if not np.isnan(v))
    y_max = max(v for v in all_vals if not np.isnan(v))
    y_range = y_max - y_min if y_max != y_min else 1
    base_y = y_min - 0.12 * y_range   # just below the bars

    # Draw candidate annotation hlines (6 predictions)
    candidate_colors = ['red', 'blue', 'purple', 'brown', 'black', 'cyan']
    candidates = [(c,d), (e,f), (g,h), (i,j), (k,l), (m,n)]
    for rank, ((xs, xe), color) in enumerate(zip(candidates, candidate_colors), 1):
        y = base_y - (rank - 1) * 0.07 * y_range
        plt.hlines(y=y, xmin=xs, xmax=xe, colors=color, linewidth=4,
                   label=f'Pred {rank}: start={xs}, end={xe}')
        plt.text(xs, y, f' {rank}', va='center', ha='left',
                 color=color, fontsize=9, fontweight='bold')

    # Draw ground truth annotations (clearly distinct: dashed, green/magenta)
    if ground_truths:
        gt_colors = ['green', 'magenta']
        gt_dash = (5, 2)   # dashed style
        n_gt = len(ground_truths)
        for gi, (gt_label, gt_start, gt_end) in enumerate(ground_truths):
            y_gt = base_y - (6 + gi) * 0.07 * y_range
            color = gt_colors[gi % len(gt_colors)]
            plt.hlines(y=y_gt, xmin=gt_start, xmax=gt_end,
                       colors=color, linewidth=5, linestyles='dashed',
                       label=f'[GROUND TRUTH] {gt_label}: {gt_start}–{gt_end}')
            plt.text(gt_start, y_gt, f' GT{gi+1}', va='center', ha='left',
                     color=color, fontsize=10, fontweight='bold')
            # Vertical span to highlight GT region on x-axis
            plt.axvspan(gt_start - 0.4, gt_end + 0.4, alpha=0.08,
                        color=color, label='_nolegend_')

    substrings = {
        (c, d): sequence[c:d+1],
        (e, f): sequence[e:f+1],
        (g, h): sequence[g:h+1],
        (i, j): sequence[i:j+1],
        (k, l): sequence[k:l+1],
        (m, n): sequence[m:n+1],
    }
    extended_positions = list(range(min(p), min(p) + len(sequence)))
    main_sequence = list(sequence)
    combined_labels = main_sequence.copy()
    for (start, end), substring in substrings.items():
        for idx in range(start, end + 1):
            combined_labels[idx] = f'{main_sequence[idx]}\n{substring[idx - start]}'

    plt.ylabel("Indication", fontsize=14)
    plt.xticks(ticks=extended_positions, labels=combined_labels, fontsize=8)
    plt.tick_params(axis='y', which='both', labelsize=14)
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelsize=8)
    plt.legend(fontsize=10, loc='upper right')
    plt.tight_layout()
    plt.savefig(infile1[:-3] + 'png', dpi=120)
    print('outfig:', infile1[:-3] + 'png')


if __name__ == "__main__":
    infile1 = '../output/3_5862187_5862334.effector_localization_indication.tsv'
    infile2 = '../output/3_5862187_5862334.effector_localization_top6_recommendation.tsv'
    sequence = 'CUUGACCUUGUAAGGCCUUUUCUUGACCUUGUAAGACCCCAUCUCUUUCUAAACGUUUUAUUAUUUUCUCGUUUUACAGAUUCUAUUCUAUCUCUUCUCAAUAUAGAAUAGAUAUCUAUCU'
    draw_6candidates_interface(infile1, infile2, sequence)
