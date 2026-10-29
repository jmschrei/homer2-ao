"""Build hooks for the homer2-ao Python package.

The wheel holds a configured HOMER tree in homer2_ao/homer/. configureHomer.pl
configures a tree by writing its absolute path into the scripts and into
cpp/SeqTag.cpp before compiling, which would tie a wheel to the directory it
was built in. This hook configures a staged copy instead, with paths found at
run time: each script finds the tree from its own location, and the compiled
programs from /proc/self/exe (or HOMER2_AO_HOME, which the launcher sets).
The repository's own files are not changed.

The metadata hook installs every HOMER program and script as a console
command of the same name.
"""

import os
import platform
import re
import shutil
import subprocess
import sys
import sysconfig
import tempfile

from hatchling.builders.hooks.plugin.interface import BuildHookInterface
from hatchling.metadata.plugin.interface import MetadataHookInterface

# What configureHomer.pl writes as lines 1-4 of each script it configures,
# with the absolute path replaced by the script's own directory.
PERL_LINES = [
	'#!/usr/bin/env perl\n',
	'use warnings;\n',
	'use Cwd (); use File::Basename (); '
	'use lib Cwd::abs_path(File::Basename::dirname(__FILE__));\n',
	'my $homeDir = Cwd::abs_path(File::Basename::dirname(__FILE__)) . "/../";\n',
]

SEQTAG_HEAD = '''#include "SeqTag.h"
//Home directory found at run time (homer2-ao): <home>/bin/<program> -> <home>/
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#ifdef __APPLE__
#include <mach-o/dyld.h>
#include <stdint.h>
#endif
static const char* homer2aoHome() {
	static char path[4096];
#ifdef __APPLE__
	char raw[4096];
	uint32_t size = sizeof(raw);
	ssize_t n = -1;
	if (_NSGetExecutablePath(raw, &size) == 0 && realpath(raw, path) != NULL) n = strlen(path);
#else
	ssize_t n = readlink("/proc/self/exe", path, sizeof(path)-2);
#endif
	if (n > 0) {
		path[n] = '\\0';
		for (int k = 0; k < 2; k++) {
			char* slash = strrchr(path, '/');
			if (slash == NULL) break;
			*slash = '\\0';
		}
		strcat(path, "/");
		return path;
	}
	const char* env = getenv("HOMER2_AO_HOME");
	return env != NULL ? env : "./";
}
const char* HomerConfig::homeDirectory = homer2aoHome();

'''

# configureHomer.pl -install writes this for a fresh tree. The wheel ships it as
# config.default.txt and the launcher copies it to config.txt on first use, so
# that pip and uv do not own config.txt and an upgrade keeps installed genomes.
CONFIG = ('# Homer Configuration File (automatically generated)\n#\n\nSOFTWARE\n'
	'homer\tv5.1\tCode/Executables, ontologies, motifs for HOMER\t'
	'http://homer.ucsd.edu/homer/data/software/homer.v5.1.zip\t./\t\n\n'
	'ORGANISMS\n\nPROMOTERS\n\nGENOMES\n')


def compiled_programs(root):
	"""Names of the programs cpp/Makefile copies into bin/."""

	names = []
	with open(os.path.join(root, 'cpp', 'Makefile')) as handle:
		for line in handle:
			m = re.match(r'\s+cp (\S+) \.\./bin/(\S*)\s*$', line)
			if m:
				names.append(m.group(2) or m.group(1))
	return names


def commands(root):
	"""Every HOMER command: the non-empty scripts in bin/, the compiled
	programs and configureHomer.pl."""

	names = set(compiled_programs(root))
	bin_dir = os.path.join(root, 'bin')
	for name in os.listdir(bin_dir):
		if name.endswith('.pl') and os.path.getsize(os.path.join(bin_dir, name)) > 0:
			names.add(name)
	names.add('configureHomer.pl')
	return sorted(names)


def configure_list(root):
	"""The scripts configureHomer.pl -make rewrites, from its @files list."""

	text = open(os.path.join(root, 'configureHomer.pl')).read()
	block = re.search(r'my @files = \((.*?)\);', text, re.S).group(1)
	return [name for name in re.findall(r'"([^"]+)"', block) if '/' not in name]


def configure_scripts(bin_dir, listed):
	"""Rewrite the scripts' first four lines as configureHomer.pl does, for
	the scripts it lists, and lines 3 and 4 of the other scripts that carry
	the same two path lines."""

	for name in sorted(os.listdir(bin_dir)):
		if not name.endswith(('.pl', '.pm')):
			continue
		path = os.path.join(bin_dir, name)
		with open(path, encoding='latin-1', newline='') as handle:
			lines = handle.readlines()

		if name in listed:
			first = range(4)
		elif (len(lines) >= 4 and lines[2].startswith('use lib "')
				and lines[3].startswith('my $homeDir = "')):
			first = (2, 3)
		else:
			continue

		for i in first:
			if i < len(lines):
				lines[i] = PERL_LINES[i]
		with open(path, 'w', encoding='latin-1', newline='') as handle:
			handle.writelines(lines)
		os.chmod(path, 0o755)


def configure_seqtag(cpp_dir):
	"""Replace SeqTag.cpp's first four lines, as configureHomer.pl does, with
	a home directory found at run time."""

	path = os.path.join(cpp_dir, 'SeqTag.cpp')
	with open(path, encoding='latin-1', newline='') as handle:
		lines = handle.readlines()
	with open(path, 'w', encoding='latin-1', newline='') as handle:
		handle.write(SEQTAG_HEAD)
		handle.writelines(lines[4:])


class CustomMetadataHook(MetadataHookInterface):
	def update(self, metadata):
		scripts = {name: 'homer2_ao._launcher:main' for name in commands(self.root)}
		scripts['homer2-ao'] = 'homer2_ao._launcher:info'
		metadata['scripts'] = scripts


class CustomBuildHook(BuildHookInterface):
	def initialize(self, version, build_data):
		if self.target_name != 'wheel':
			return

		for tool in ('make', 'perl'):
			if shutil.which(tool) is None:
				raise RuntimeError('homer2-ao needs {} to build and run'.format(tool))

		self.stage = tempfile.mkdtemp(prefix='homer2-ao-build-')
		home = os.path.join(self.stage, 'homer')
		root = self.root
		ignore_cpp = shutil.ignore_patterns('backup', '*.o', *compiled_programs(root))

		os.makedirs(os.path.join(home, 'bin'))
		for name in os.listdir(os.path.join(root, 'bin')):
			if name.endswith(('.pl', '.pm')):
				shutil.copy2(os.path.join(root, 'bin', name), os.path.join(home, 'bin', name))
		shutil.copytree(os.path.join(root, 'cpp'), os.path.join(home, 'cpp'), ignore=ignore_cpp)
		shutil.copytree(os.path.join(root, 'data'), os.path.join(home, 'data'),
			ignore=shutil.ignore_patterns('genomes', 'promoters', 'organisms'))
		shutil.copytree(os.path.join(root, 'motifs'), os.path.join(home, 'motifs'))
		for name in ('configureHomer.pl', 'COPYING', 'README.txt'):
			shutil.copy2(os.path.join(root, name), os.path.join(home, name))
		with open(os.path.join(home, 'config.default.txt'), 'w') as handle:
			handle.write(CONFIG)

		configure_scripts(os.path.join(home, 'bin'), set(configure_list(root)))
		configure_seqtag(os.path.join(home, 'cpp'))

		jobs = os.environ.get('HOMER2_AO_BUILD_JOBS', str(min(4, os.cpu_count() or 1)))
		cpp = os.path.join(home, 'cpp')
		subprocess.run(['make', '-C', cpp, 'clean'], check=True, stdout=subprocess.DEVNULL)
		make = ['make', '-C', cpp, '-j', jobs]
		flags = []
		if sys.platform == 'darwin':
			# Motif2.cpp needs C++17, which Apple's compiler may not default to
			flags.append('-std=gnu++17')
		elif platform.machine() not in ('x86_64', 'AMD64', 'i686'):
			# HOMER's code assumes char is signed, as on x86; on ARM Linux it is not
			flags.append('-fsigned-char')
		if flags:
			make.append('COMPILER=g++ ' + ' '.join(flags))
		subprocess.run(make, check=True)
		missing = [name for name in compiled_programs(root)
			if not os.path.exists(os.path.join(home, 'bin', name))]
		if missing:
			raise RuntimeError('the build did not produce ' + ', '.join(missing))
		subprocess.run(['make', '-C', cpp, 'clean'], check=True, stdout=subprocess.DEVNULL)

		# -g's debug sections hold the build directory's path and most of
		# each program's size; the code itself is unchanged.
		if shutil.which('strip'):
			for name in compiled_programs(root):
				path = os.path.join(home, 'bin', name)
				if subprocess.run(['strip', '--strip-debug', path]).returncode != 0:
					subprocess.run(['strip', '-S', path])

		build_data['force_include'][home] = 'homer2_ao/homer'
		build_data['pure_python'] = False
		platform = sysconfig.get_platform().replace('-', '_').replace('.', '_')
		build_data['tag'] = 'py3-none-' + platform

	def finalize(self, version, build_data, artifact_path):
		if getattr(self, 'stage', None):
			shutil.rmtree(self.stage, ignore_errors=True)
