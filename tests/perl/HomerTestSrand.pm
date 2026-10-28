package HomerTestSrand;
# Loaded at the start of every perl process the tests run, through
# PERL5OPT=-MHomerTestSrand (tests/homer_helpers.py): HOMER's scripts call
# rand() without srand(), so Perl seeds them from the clock. A fixed seed makes
# background selection, scrambled sequences and temporary names repeatable.
# Perl 5.20 and later use their own drand48 on every platform, so the draws are
# the same everywhere.
srand(1);
1;
