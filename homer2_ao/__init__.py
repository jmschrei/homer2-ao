"""homer2-ao: HOMER v5.1 with faster motif scanning and motif finding.

The installed package holds a configured HOMER tree in `homer/`. Every HOMER
program is installed as a command of the same name; `home()` and `bin_dir()`
give the tree's location for code that calls the programs directly.
"""

import os

__version__ = "1.0.0"


def home():
	"""The directory of the installed HOMER tree (bin/, data/, motifs/,
	config.txt)."""

	return os.path.join(os.path.dirname(os.path.abspath(__file__)), 'homer')


def bin_dir():
	"""The directory holding HOMER's programs and scripts."""

	return os.path.join(home(), 'bin')
