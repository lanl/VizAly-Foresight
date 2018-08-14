@echo off
REM Run dvipdfm from dvipdf-style command-line
REM Assume no options specified
REM If this script is called from latexmk,
REM    we have %1=source.dvi, %2=dest.pdf
REM But for safety, let's handle correctly a one argument call,
REM    i.e., %1=source, with no %2

if "%2" == "" goto onearg

:twoarg
dvipdfm -o %2 %1
goto done

:onearg
dvipdfm %1


:done