# Sometime in the future, latexmk will directly support the use of
# TeX4ht to obtain html from tex.  Meanwhile, here is how to use
# latexmk with TeX4ht.  There is a script htlatex supplied by the
# TeX4ht package: It simply runs latex a fixed number of times and
# then the programs tex4ht and t4ht.  To use latexmk to get optimal
# processing use the following instructions (under UNIX-like operating
# systems, e.g., OS-X and linux):
#
#  1. Put the scripts htlatexonly and myhtlatex2 somewhere in the PATH
#     for executables (and make sure they have excutable permissions
#     set).
#  2. Set up an initialization file for latexmk like this one.
#
#  3. To process file.tex to make file.html, run
#
#             myhtlatex2 file
#

# Since these instructions use scripts that are UNIX shell scripts,
# the instructions work as written for UNIX-like operating
# systems. Users of other operating systems will have to adjust them
# and modify the scripts suitably.


warn "latexmkrc for htlatex\n";

$dvi_mode = 1;
$pdf_mode = 0;
$quote_filenames = 0;
$latex = 'htlatexonly %S';

$clean_ext .= ' 4ct 4tc idv lg tmp xref';
$clean_full_ext .= ' css html';
