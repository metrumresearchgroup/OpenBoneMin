SHELL := /bin/bash
PACKAGE=OpenBoneMin
VERSION=$(shell grep Version DESCRIPTION |awk '{print $$2}')
TARBALL=${PACKAGE}_${VERSION}.tar.gz
PKGDIR=.
CHKDIR=Rchecks

## Set libPaths:
#export R_LIBS=${LIBDIR}

readme:
	Rscript -e "rmarkdown::render('README.Rmd')"

ec:
	echo ${VERSION}

all:
	make doc
	make build
	make install

.PHONY: doc
doc:
	Rscript -e 'library(devtools); document()'

build:
	R CMD build --md5 $(PKGDIR)


install:
	R CMD INSTALL --install-tests ${TARBALL}

install-build:
	R CMD INSTALL --build --install-tests ${TARBALL}

check:
	make doc
	make build
	R CMD CHECK ${TARBALL} -o ${CHKDIR}



