"""Console entry points for HOMER's programs.

Every HOMER program and script is installed as a command of the same name
that points here. The launcher replaces itself with the real program in the
installed tree, with that tree's bin/ first on PATH, because HOMER's scripts
call each other and the compiled programs by bare name.
"""

import os
import sys

from . import __version__, bin_dir, home


def _exec(path, args):
	env = dict(os.environ)
	env['PATH'] = bin_dir() + os.pathsep + env.get('PATH', '')
	env['HOMER2_AO_HOME'] = home()

	try:
		os.execve(path, [path] + list(args), env)
	except FileNotFoundError:
		sys.exit('homer2-ao: {} is not in this installation'.format(path))


def main():
	"""Run the HOMER program this command is named after."""

	name = os.path.basename(sys.argv[0])
	if name == 'configureHomer.pl':
		path = os.path.join(home(), name)
	else:
		path = os.path.join(bin_dir(), name)

	_exec(path, sys.argv[1:])


def info():
	"""`homer2-ao`: print the version and where HOMER is installed, or run a
	HOMER program by name (`homer2-ao <program> [args]`)."""

	args = sys.argv[1:]
	if not args or args[0] in ('-h', '--help'):
		print('homer2-ao {}'.format(__version__))
		print('HOMER directory: {}'.format(home()))
		print('Programs:        {}'.format(bin_dir()))
		print()
		print('Every HOMER program is installed as a command of the same name,')
		print('e.g. findMotifsGenome.pl or homer2; `homer2-ao <program> [args]`')
		print('runs one too. If you use this software, cite the HOMER paper:')
		print('Heinz S, Benner C, et al. Mol Cell 2010;38(4):576-589.')
		return
	if args[0] == '--version':
		print(__version__)
		return
	if args[0] == '--home':
		print(home())
		return

	name = args[0]
	if name == 'configureHomer.pl':
		_exec(os.path.join(home(), name), args[1:])
	_exec(os.path.join(bin_dir(), name), args[1:])
