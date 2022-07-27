#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logomaker

def main():

	print('Create a PSSWM logo from command-line inputs. Press Ctrl+d after pasting data')

	# Logo maker accepts a pandas DataFrame with bases (A, T, G, C) as columns and positions as rows
	# Read in the data from stdin and store as a dict(), we can then convert this dict() into a pandas DataFrame.
	# The keys are the column headers and the values are lists.
	data = dict()

	num_col = -1

	# Read in data with bases on the rows and position along the columns.
	# 	A	1	3	-2	-1	...
	# 	T	-1	4	0	6	...
	# 	G	4	1	7	3	...
	# 	C	2	1	8	1	...
	for line in sys.stdin:
		
		col = line.strip().split()

		if num_col < 0:
			num_col = len(col)
		else:
			if num_col != len(col):

				print( 'Did not read the expected number of columns for {}'.format(line) )
				sys.exit(1)

		# The data is in string form -- convert to float for PSWM values
		data[ col[0] ] = [ float(i) for i in col[1:] ]

	# Force scaling by adding a dummy column
	force_scaling = False

	max_y = 1.5
	min_y = -0.9

	if force_scaling:

		num_col += 1
		data['A'].append(max_y)
		data['T'].append(min_y)
		data['G'].append(0)
		data['C'].append(0)
	
	data_df = pd.DataFrame(data=data)

	# Create the Logo object
	motif_logo = logomaker.Logo(data_df, width=.8, vpad=.05)

	# style using Logo methods
	motif_logo.style_spines(spines=['left', 'right'], visible=False)

	# style using Axes methods
	motif_logo.ax.set_xticks(range(len(data_df)))

	x_labels = list()

	pos = -( (num_col - 1)//2 )

	if force_scaling:
		pos += 1

	for i in range(num_col - 1):

		x_labels.append( str(pos) )
		pos += 1

	motif_logo.ax.set_xlabel('position')
	motif_logo.ax.set_xticklabels(x_labels)
	#motif_logo.ax.set_yticks([0, .5, 1])
	#motif_logo.ax.axvline(2.5, color='k', linewidth=1, linestyle=':')
	motif_logo.ax.set_ylabel('score')

	# Set the initial size of the displayed window
	motif_logo.fig.set_size_inches(5, 5)

	# When *not* running in a jupyter notebook, we need to explicitly call plt.show() to 
	# display the plot in a new window.
	plt.show()

main()
