#!/usr/bin/env bash

for file in qc_trim filter_reads; do
	outlabel=$(echo ${file^^}|sed 's/_/\\_/g')
	pod2latex -h1level 3 -modify --out stages/${file}.tex ../bin/${file}
 	cat stages/${file}.tex |  perl -ne "s/label\{(SYNOPSIS|DESCRIPTION|REQUIRED_ARGUMENTS|OPTIONAL_ARGUMENTS)\}/label\{${outlabel}_\$1\}/g; print" > stages/${file}.tex.tmp
	mv stages/${file}.tex.tmp stages/${file}.tex
done

pdflatex SCaMP_User_Guide.tex

