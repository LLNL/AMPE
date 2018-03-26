SHARED_DOC = ../../shared/doc
all:
	@echo 'Usage: make doc'
	@echo '       make clean'
	@echo
doc:
	latex toms_cvodes
	bibtex toms_cvodes
	latex toms_cvodes
	latex toms_cvodes
	latex toms_cvodes
	dvips -Ppdf -o toms_cvodes.ps toms_cvodes
	ps2pdf toms_cvodes.ps toms_cvodes.pdf

clean:
	rm -f *~ *.aux *.bbl *.blg *.log *.toc *.lot *.lof *.ind *.ilg *.idx *.dvi
	rm -f toms_cvodes.ps toms_cvodes.pdf

