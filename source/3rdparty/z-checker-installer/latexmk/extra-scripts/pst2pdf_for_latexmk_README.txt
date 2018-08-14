Pdftricks implementation for TeXShop.

The files ``pdftricksmkrc'' and ``pst2pdf_for_latexmk''  (as well as latexmk, etc.) are stored in ~/Library/TeXShop/bin/

The file pdftrciskmk.engine is stored in ~/Library/TeXShop/Engines/.

When instructed to Typeset a the file TeXShop cd's to the tex files directory and calls pdftricksmk.engine passing the file name (with extension). This executes latexmk using the pdftricksmkrc file for initialization.

The call to pdflatex in the rc file DOES use shell escape (just in case eps files are also being input using epstopdf) so you MUST use the [noshell] option for pdftricks (\usepackage[noshell]{pdftricks}) to avoid a run condition.

The processing steps I use for the -fig* files differs from that used in the standard pdftricks (and what you used). We've discovered that the original processing sometimes rotates figures and also sometimes produces a BoundingBox that cuts off descenders on letters at the edge of the figure. (I've changed my pdftricks package to use these steps too and have let the authors of that package know.)

So far it seems to work. I'll get to work on a pst-pdf package example but, now that I'm beginning to understand what is happening with your new extensions to latexmk, I also suspect it should be fairly easy to accomodate that package. I don't see a case where anyone would use both pdftricks and pst-pdf while I do see cases where eps images, included using epstopdf, will be mixed with pstricks type images.
