"""Console entry points for HOMER's programs.

Every HOMER program and script is installed as a command of the same name
that points here. The launcher replaces itself with the real program in the
installed tree, with that tree's bin/ first on PATH, because HOMER's scripts
call each other and the compiled programs by bare name.
"""

import os
import sys
import shutil

from . import __version__, bin_dir, home


def _command(name, args):
	"""The command line that runs HOMER program `name`. configureHomer.pl is
	not executable in HOMER's release, which runs it as `perl configureHomer.pl`."""

	if name == 'configureHomer.pl':
		return ['perl', os.path.join(home(), name)] + list(args)
	return [os.path.join(bin_dir(), name)] + list(args)


def ensure_config():
	"""Write config.txt from config.default.txt if it is missing. HOMER's
	scripts and programs stop without it, and configureHomer.pl -install
	records genomes in it, so it is not part of the installed files: an
	upgrade then keeps it, and the genomes it lists."""

	path = os.path.join(home(), 'config.txt')
	if os.path.exists(path):
		return
	try:
		tmp = '{}.{}.tmp'.format(path, os.getpid())
		shutil.copyfile(os.path.join(home(), 'config.default.txt'), tmp)
		os.replace(tmp, path)
	except OSError:
		pass


def _exec(argv):
	ensure_config()
	env = dict(os.environ)
	env['PATH'] = bin_dir() + os.pathsep + env.get('PATH', '')
	env['HOMER2_AO_HOME'] = home()

	try:
		os.execvpe(argv[0], argv, env)
	except FileNotFoundError:
		sys.exit('homer2-ao: cannot run {}: not found'.format(argv[0]))


def main():
	"""Run the HOMER program this command is named after."""

	_exec(_command(os.path.basename(sys.argv[0]), sys.argv[1:]))


def info():
	"""`homer2-ao`: print the version and where HOMER is installed, or run a
	HOMER program by name (`homer2-ao <program> [args]`)."""

	args = sys.argv[1:]
	ensure_config()
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

	_exec(_command(args[0], args[1:]))
