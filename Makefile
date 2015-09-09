default : attrs install

attrs :
	Rscript -e "library(Rcpp) ; compileAttributes(verbose=TRUE)"

build : attrs clean
	R CMD build . --resave-data
	R CMD Rd2pdf -o doc.pdf .

check : build
	R CMD check emIRT*.tar.gz

fullinstall : build
	R CMD INSTALL emIRT*.tar.gz

install : attrs clean
	R CMD INSTALL .

remove :
	R CMD REMOVE emIRT

clean :
	rm -f emIRT*.tar.gz
	rm -fr pkg.Rcheck
	rm -f ./src/*.o
	rm -f ./src/*.so
	rm -f ./src/*.rds
	rm -f ./inst/lib/*
	rm -f ./doc.pdf
