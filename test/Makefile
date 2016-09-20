### external setting of path to ell_type.m
ELLEVAL?=../

### setting default paths of external programs
OCTAVE?=/opt/octave-4.0.3/


export PATH:= $(OCTAVE)/bin:$(PATH)



SHELL:= /bin/bash


TESTf:= t00.efit t01.efit t02.efit t03.efit

# ETres:= .ells -3Dsym.svg -3Dasym.svg -2Dpdist.svg -2Dhist.svg -2Danno.svg -2DpdistOV.svg -2Dshist.svg
# TESTfET:= $(foreach testf,$(TESTf),$(ETres:%=$(testf)%))

TESTex:= 3Dasym_VB.svg 2Dpdist_VB.svg 2Dshist_inc-shearY_VB.svg 2Dshist_inc-pdist-shearY_VB.svg 2Dshist_inc-shearY_VB-VB.svg 2Dshist_inc-pdist-shearY_VB-VB.svg
TESTexs:= $(foreach testf,$(TESTf),$(TESTex:%=$(testf)-%))

TEST:= 3Dasym.svg 2Dpdist.svg $(TESTex)
TESTs:= $(foreach testf,$(TESTf),$(TEST:%=$(testf)-%))



.PHONY: test

test : $(TESTs)

# $(TESTfET) : # static pattern rule not possible for multi-pattern rule
%.efit.ells \
%.efit-3Dsym.svg \
%.efit-3Dasym.svg \
%.efit-2Dpdist.svg \
%.efit-2Dhist.svg \
%.efit-2Danno.svg \
%.efit-2DpdistOV.svg \
%.efit-2Dshist.svg \
 : %.efit
	octave-cli $(ELLEVAL)ell_type.m $< 2.0


%.efit-3Dasym_VB.svg : %.efit-3Dasym.svg
%.efit-2Dpdist_VB.svg : %.efit-2Dpdist.svg
%.efit-2Dshist_inc-pdist-shearY_VB.svg : %.efit-2Dshist.svg
$(TESTexs) : # target list must only include extended test files i.e. based on TESTex
	$(MAKE) -f ../Makefile ELLEVAL=./  $@
