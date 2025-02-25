import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

bar_width = 0.35

def draw_6candidates_interface(infile1, infile2, sequence):
    with open (infile2, 'r') as fh:
        fh.readline()#Position p	Start S(p)	End E(p)	Sum of Indications	Best length	End position	ML predicted #[0, 5]        
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
    # Plot S(p) and E(p)
    plt.figure(figsize=(16, 6))
    plt.bar(positions_start, Start, width=bar_width, label="Start S(p)")
    plt.bar(positions_end, End, width=bar_width, label="End E(p)", color='orange')

    # Draw annotation hlines
    hlines_positions = [
        (10*2-900, c, d, '1st candidate start {c}, end {d}', 'red'),
        (20*2-900, e, f, '2nd candidate start {e}, end {f}', 'blue'),
        (30*2-900, g, h, '3rd candidate start {g}, end {h}', 'purple'),
        (40*2-900, i, j, '4th candidate start {i}, end {j}', 'brown'),
        (50*2-900, k, l, '5th candidate start {k}, end {l}', 'black'),
        (60*2-900, m, n, '6th candidate start {m}, end {n}', 'cyan')
    ]
    for y, xmin, xmax, label, color in hlines_positions:
        plt.hlines(y=y, xmin=xmin, xmax=xmax, colors=color, linewidth=3, 
                   label=f'{label.format(a=xmin, b=xmax, c=c, d=d, e=e, f=f, g=g, h=h, i=i, j=j, k=k, l=l, m=m, n=n)}')

    substrings = {
        (c, d): sequence[c:d+1],
        (e, f): sequence[e:f+1],
        (g, h): sequence[g:h+1],
        (i, j): sequence[i:j+1],
        (k, l): sequence[k:l+1],
        (m, n): sequence[m:n+1],
    }
    # Display the sequence and substrings on the x-axis
    extended_positions = list(range(min(p), min(p) + len(sequence)))  # Extend ticks for full sequence length
    main_sequence = list(sequence)
    combined_labels = main_sequence.copy()
    for (start, end), substring in substrings.items():
        for i in range(start, end + 1):
            combined_labels[i] = f'{main_sequence[i]}\n{substring[i - start]}'

    plt.ylabel("Indication", fontsize=14)
    plt.xticks(ticks=extended_positions, labels=combined_labels, fontsize=10)
    plt.tick_params(axis='y', which='both', labelsize=14)
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelsize=14)
    plt.legend(fontsize=11, loc='upper left')
    plt.tight_layout()
    plt.savefig(infile1[:-3] + 'png')
    print('outfig:', infile1[:-3] + 'png')














if __name__ == "__main__":
  infile1 = '../output/3_5862187_5862334.effector_localization_indication.tsv'
  infile2 = '../output/3_5862187_5862334.effector_localization_top6_recommendation.tsv'
  sequence = 'CUUGACCUUGUAAGGCCUUUUCUUGACCUUGUAAGACCCCAUCUCUUUCUAAACGUUUUAUUAUUUUCUCGUUUUACAGAUUCUAUUCUAUCUCUUCUCAAUAUAGAAUAGAUAUCUAUCU'
  draw_6candidates_interface(infile1, infile2, sequence)
