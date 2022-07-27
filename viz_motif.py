#!/usr/bin/env python3

# Read k-mer motif to score mappings from stdin and plot the scores of each parent motif for each allowed mutation type
# as a heat map.

import sys
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# The threshold for determining equal rank scores
EPS = 1.0e-7

def main():

    output_filename = ''

    if ( len(sys.argv) == 2 ) and ('.png' in sys.argv[1]):
        output_filename = sys.argv[1]

    # Specify the display order for mutations. Since the row labels are sequences that are sorted in
    # lexographic order, the allowed substitutions are also specified in lexographic order
    allowed_mutations = ['->A', '->C', '->G', '->T',  'del', 'ins']
    #allowed_mutations = ['->A', '->C', '->G', '->T']

    print('Reading k-mer score mappings from STDIN. Press Ctrl+d to finish', file=sys.stderr)

    # Store k-mer score mappings as nested dictionaries: mutation -> codon -> region -> kmer -> score
    motif = dict()

    mutation = ''
    codon = ''
    region = ''

    # The number of bases in the k-mer
    motif_len = -1

    # Validate the inputs, count the number of codons and regions
    allowed_codons = set()
    allowed_regions = set()
    allowed_score = set()

    # Read the kmer motifs from stdin
    for line in sys.stdin:

        line = line.strip()

        if 'motif for mutation' in line:

            match = re.search('motif for mutation (.*); codon loc (\d+); region (\d+)', line)

            if not match:
                throw_error('Could not parse motif header')
            
            mutation = match.group(1)
            codon = int(match.group(2))
            region = int(match.group(3))

            allowed_codons.add(codon)
            allowed_regions.add(region)

            if mutation not in allowed_mutations:
                throw_error('Did not find {} in the list of allowed mutations'.format(mutation))

        match = re.search('([A|T|G|C]+)\t(.*)', line)

        if match:
            
            if (mutation == '') or (codon == '') or (region == ''):
                throw_error('Invalid motif header information') 

            seq = match.group(1)

            # Convert the score from a string to a float!
            score = float( match.group(2) )

            allowed_score.add(score)

            if motif_len == -1:
                motif_len = len(seq)
            elif motif_len != len(seq):
                throw_error('Unexpected motif length variation!')
            
            # mutation -> codon -> region -> base -> list of scores
            if mutation not in motif:
                motif[mutation] = dict()
            
            if codon not in motif[mutation]:
                motif[mutation][codon] = dict()
            
            if region not in motif[mutation][codon]:
                motif[mutation][codon][region] = dict()
            
            if seq not in motif[mutation][codon][region]:
                motif[mutation][codon][region][seq] = dict()
            
            motif[mutation][codon][region][seq] = score

    # Convert from sets to lists
    allowed_regions = list(allowed_regions)
    allowed_codons = list(allowed_codons)
    
    # Associate each score with a rank
    allowed_score = list(allowed_score)
    allowed_score.sort()

    score_rank = dict()

    # Only populate the score_rank dictionary when plotting by the rank of each score
    #for i in range( len(allowed_score) ):
    #    score_rank[allowed_score[i]] = i

    if (len(allowed_regions) == 1) and (len(allowed_codons) == 1):
        plot_simple(motif, score_rank, motif_len, allowed_mutations)
    elif (len(allowed_regions) == 1) and (len(allowed_codons) > 1):
        plot_by_codon(motif, score_rank, motif_len, allowed_mutations, allowed_codons)
    elif (len(allowed_regions) > 1) and (len(allowed_codons) == 1):
        plot_by_region(motif, score_rank, motif_len, allowed_mutations, allowed_regions)

def plot_simple(m_motif, m_score_rank, m_motif_len, m_allowed_mutations):

    if m_motif_len == 1:
        plot_simple_single(m_motif, m_score_rank, m_allowed_mutations)
    else: # m_motif_len > 1
        plot_simple_by_central_base(m_motif, m_score_rank, m_allowed_mutations)

def plot_simple_single(m_motif, m_score_rank, m_allowed_mutations):

    if len(m_score_rank) > 0:
        plot_by_rank = True
    else:
        plot_by_rank = False

    # Plot a single 4 by 6 plot
    fig, heatmap = plt.subplots(nrows=1,ncols=1)

    seq = ['A', 'C', 'G', 'T']
    
    data = np.zeros( (len(seq), len(m_allowed_mutations) ) )

    for s in range( len(seq) ):
        for c in range( len(m_allowed_mutations) ):
            if s == c:
                data[s][c] = np.nan # Mask self-substitution
            else:

                if plot_by_rank:
                    data[s][c] = m_score_rank[ m_motif[ m_allowed_mutations[c] ][0][0][seq[s]] ]
                else:
                    data[s][c] = m_motif[ m_allowed_mutations[c] ][0][0][seq[s]]

    if plot_by_rank:
        im = heatmap.imshow(data, cmap='hot', interpolation='none', aspect=0.5, vmin=0, vmax=len(m_score_rank))
    else:
        im = heatmap.imshow(data, cmap='hot', interpolation='none', aspect=0.5, vmin=-1.0, vmax=1.0)

    heatmap.axis('scaled')

    heatmap.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)
    heatmap.xaxis.set_major_locator( mticker.FixedLocator( [loc for loc in range(len(m_allowed_mutations))] ) )
    heatmap.set_xticklabels(m_allowed_mutations, rotation=0)

    heatmap.yaxis.set_major_locator( mticker.FixedLocator( [loc for loc in range( len(seq) )] ) )
    heatmap.set_yticklabels(seq)

    # Since masked values are not plotted, the background color determines the apparent color of 
    # masked cells
    heatmap.set_facecolor('black')
    
    #annotate_heatmap(im, valfmt='{x:.0f}')

    fig.colorbar(im, ax=heatmap, orientation='vertical', pad=0.05, format=mticker.FuncFormatter(lambda x, pos: '{0:.2f}'.format(x)))
   
    plt.show()

def plot_simple_by_central_base(m_motif, m_score_rank, m_allowed_mutations):

    # Plot a 1 x 4 grid of heatmaps, where each heat map is 4^(m_motif_len - 1) by 5 plot    
    if len(m_score_rank) > 0:
        plot_by_rank = True
    else:
        plot_by_rank = False

    fig, heatmap = plt.subplots(nrows=1,ncols=4)

    base = ['A', 'C', 'G', 'T']
    
    for central_base_index in range( len(base) ):

        # Extract the 4^(m_motif_len - 1) kmers that will be plotted in this graph
        kmer = list()

        for k,v in m_motif[ m_allowed_mutations[0] ][0][0].items():
            if k[ int(len(k)/2) ] == base[central_base_index]:
                kmer.append(k)
        
        kmer.sort()

        # Exclude the substitution of a base to itself
        allowed_mutations = [m for m in m_allowed_mutations if '->' + base[central_base_index] not in m]

        data = np.zeros( (len(kmer), len(allowed_mutations) ) )

        for s in range( len(kmer) ):
            for c in range( len(allowed_mutations) ):
                if plot_by_rank:
                    data[s][c] = m_score_rank[ m_motif[ allowed_mutations[c] ][0][0][kmer[s]] ]
                else:
                    data[s][c] = m_motif[ allowed_mutations[c] ][0][0][kmer[s]]

        if plot_by_rank:
            im = heatmap[central_base_index].imshow(data, cmap='hot', interpolation='none', aspect=0.5, vmin=0, vmax=len(m_score_rank))
        else:
            im = heatmap[central_base_index].imshow(data, cmap='hot', interpolation='none', aspect=0.5, vmin=-1, vmax=1.0)
        
        heatmap[central_base_index].axis('scaled')

        heatmap[central_base_index].tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)
        heatmap[central_base_index].xaxis.set_major_locator( mticker.FixedLocator( [loc for loc in range(len(allowed_mutations))] ) )
        heatmap[central_base_index].set_xticklabels(allowed_mutations, rotation=0)

        heatmap[central_base_index].yaxis.set_major_locator( mticker.FixedLocator( [loc for loc in range(len(kmer))] ) )
        heatmap[central_base_index].set_yticklabels(kmer)

        # Since masked values are not plotted, the background color determines the apparent color of 
        # masked cells
        heatmap[central_base_index].set_facecolor('black')

        #annotate_heatmap(im, valfmt='{x:.0f}')

    #fig.colorbar(im, ax=heatmap, orientation='vertical', pad=0.05, format=mticker.FuncFormatter(lambda x, pos: '{0:.0f}'.format(x)))
   
    plt.show()

def plot_by_codon(m_motif, m_score_rank, m_motif_len, m_allowed_mutations, m_allowed_codons):

    if m_motif_len == 1:
        plot_simple_by_codon(m_motif, m_score_rank, m_allowed_mutations)
    else: # m_motif_len > 1
        plot_by_central_base_and_codon(m_motif, m_score_rank, m_allowed_mutations)

def plot_simple_by_codon(m_motif, m_score_rank, m_allowed_mutations):

    # Plot a 2 x 2 grid of heatmaps, where each heat map is 4 by 6 plot
    if len(m_score_rank) > 0:
        plot_by_rank = True
    else:
        plot_by_rank = False

    fig, heatmap = plt.subplots(nrows=2,ncols=2)

    seq = ['A', 'C', 'G', 'T']

    for r in [0, 1]:
        for c in [0, 1]:

            data = np.zeros( (len(seq), len(m_allowed_mutations) ) )
            codon = r*2 + c

            for i in range( len(seq) ):
                for j in range( len(m_allowed_mutations) ):

                    if plot_by_rank:
                        data[i][j] = m_score_rank[ m_motif[ m_allowed_mutations[j] ][codon][0][seq[i]] ]
                    else:
                        data[i][j] = m_motif[ m_allowed_mutations[j] ][codon][0][seq[i]]

                    # Mask self substitution
                    if i == j:
                        data[i][j] = np.nan

            if plot_by_rank:
                im = heatmap[r][c].imshow(data, cmap='hot', interpolation='none', aspect=0.5, vmin=0, vmax=len(m_score_rank))
            else:
                im = heatmap[r][c].imshow(data, cmap='hot', interpolation='none', aspect=0.5, vmin=-1, vmax=1.0)
            
            heatmap[r][c].axis('scaled')
            #heatmap[r][c].set_aspect('auto') #Allow rectangular cells (instead of square)

            heatmap[r][c].tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)
            heatmap[r][c].xaxis.set_major_locator( mticker.FixedLocator( [loc for loc in range(len(m_allowed_mutations))] ) )
            heatmap[r][c].set_xticklabels(m_allowed_mutations, rotation=0)

            heatmap[r][c].yaxis.set_major_locator( mticker.FixedLocator( [loc for loc in range(len(seq))] ) )
            heatmap[r][c].set_yticklabels(seq)

            # Since masked values are not plotted, the background color determines the apparent color of 
            # masked cells
            heatmap[r][c].set_facecolor('black')

            #annotate_heatmap(im, valfmt='{x:.0f}')

    plt.show()

def plot_by_central_base_and_codon(m_motif, m_score_rank, m_allowed_mutations):

    # Plot a 4 x 4 grid of heatmaps, where each heat map is 4^(kmer-1) by 6 plot:
    # 
    #   [codon 0] [codon 1] [codon 2] [codon 3] <-- central base A
    #   [codon 0] [codon 1] [codon 2] [codon 3] <-- central base T
    #   [codon 0] [codon 1] [codon 2] [codon 3] <-- central base G
    #   [codon 0] [codon 1] [codon 2] [codon 3] <-- central base C

    if len(m_score_rank) > 0:
        plot_by_rank = True
    else:
        plot_by_rank = False

    fig, heatmap = plt.subplots(nrows=4,ncols=4)

    base = ['A', 'C', 'G', 'T']

    for central_base_index in range( len(base) ):

        for codon in [0, 1, 2, 3]:

            # Extract the 4^(m_motif_len - 1) kmers that will be plotted in this graph
            kmer = list()

            for k,v in m_motif[ m_allowed_mutations[0] ][codon][0].items():
                if k[ int(len(k)/2) ] == base[central_base_index]:
                    kmer.append(k)
            
            kmer.sort()

            # Exclude the substitution of a base to itself
            allowed_mutations = [m for m in m_allowed_mutations if '->' + base[central_base_index] not in m]

            data = np.zeros( (len(kmer), len(allowed_mutations) ) )

            for s in range( len(kmer) ):
                for c in range( len(allowed_mutations) ):
                    if plot_by_rank:
                        data[s][c] = m_score_rank[ m_motif[ allowed_mutations[c] ][codon][0][kmer[s]] ]
                    else:
                        data[s][c] = m_motif[ allowed_mutations[c] ][codon][0][kmer[s]]

            if plot_by_rank:
                im = heatmap[central_base_index][codon].imshow(data, cmap='hot', interpolation='none', aspect=0.5, vmin=0, vmax=len(m_score_rank))
            else:
                im = heatmap[central_base_index][codon].imshow(data, cmap='hot', interpolation='none', aspect=0.5, vmin=-1, vmax=1.0)
            
            heatmap[central_base_index][codon].axis('scaled')
            heatmap[central_base_index][codon].set_aspect('auto') #Allow rectangular cells (instead of square)

            if central_base_index == 0:

                heatmap[central_base_index][codon].tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)
                heatmap[central_base_index][codon].xaxis.set_major_locator( mticker.FixedLocator( [loc for loc in range(len(allowed_mutations))] ) )
                heatmap[central_base_index][codon].set_xticklabels(allowed_mutations, rotation=0)
            else:
                heatmap[central_base_index][codon].tick_params(top=False, bottom=False, labeltop=False, labelbottom=False)
                heatmap[central_base_index][codon].xaxis.set_major_locator( mticker.FixedLocator( [] ) )
                heatmap[central_base_index][codon].set_xticklabels([], rotation=0)

            # Only label the y-axis for the left hand heat maps
            if codon == 0:
                heatmap[central_base_index][codon].yaxis.set_major_locator( mticker.FixedLocator( [loc for loc in range(len(kmer))] ) )
                heatmap[central_base_index][codon].set_yticklabels(kmer)
            else:
                heatmap[central_base_index][codon].yaxis.set_major_locator( mticker.FixedLocator( [] ) )
                heatmap[central_base_index][codon].set_yticklabels([])

            # Since masked values are not plotted, the background color determines the apparent color of 
            # masked cells
            heatmap[central_base_index][codon].set_facecolor('black')

            #annotate_heatmap(im, valfmt='{x:.0f}')

    #fig.colorbar(im, ax=heatmap, orientation='vertical', pad=0.05, format=mticker.FuncFormatter(lambda x, pos: '{0:.0f}'.format(x)))
   
    plt.show()

def throw_error(m_msg):

    print('Caught the error: {}'.format(m_msg))
    sys.exit(1)

def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=('white', 'black'),
                     threshold=None, **textkw):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A pair of colors.  The first is used for values below a threshold,
        the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **textkw
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = mticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):

            # Don't label masked elements
            if np.isnan(data[i, j]):
                continue
            
            kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts

main()