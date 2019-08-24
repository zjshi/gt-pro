all:  algorithms.pdf
	echo "Done making algorithms.pdf."
	open algorithms.pdf || echo "Have a nice day."

algorithms.pdf: algorithms.tex algorithms.toc
	xelatex algorithms.tex

algorithms.toc: algorithms.tex
	xelatex algorithms.tex
