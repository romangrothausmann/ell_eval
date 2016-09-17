
### external setting of path to ell_type.m
ELLEVAL?=./


### setting default paths of external programs
OCTAVE?=/opt/octave-4.0.3/


export PATH:= $(OCTAVE)/bin:$(PATH)



SHELL:= /bin/bash


.PHONY: test


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


%_shearY.svg : %.svg
	sed '/<rect x="0" y="0" width=".*" height=".*" fill="none"\/>/d' $< > $@ # remove BG rec, not possible from within octave as gnuplot term is overwritten by print-command? (not tried) http://stackoverflow.com/questions/18169221/gnuplot-png-file-without-border-line
	## get BBox coords with inkscape ;-) https://sourceforge.net/p/inkscape/mailman/message/21202829/ # not available in rsvg-convert from librsvg nor svg_utils https://github.com/btel/svg_utils
	## SVG transform https://www.w3.org/TR/SVG/coords.html#TransformAttributeEffectOnSiblingAttributes
	x=`inkscape --without-gui --query-id=gnuplot_canvas -X $@ ` ; \
	y=`inkscape --without-gui --query-id=gnuplot_canvas -Y $@ ` ; \
	sed -i 's/<g id="gnuplot_canvas">/<g id="gnuplot_canvas" transform="translate('$$x,$$y') skewY(-30) translate('-$$x,-$$y')" >/' $@ # basically works, ToDo: do not plot text (overlay in latex) or also repos "separation curve", "ellipse arc", "oblate line", "circle\npoint" # https://developer.mozilla.org/de/docs/Web/SVG/Attribute/transform

%_VB.svg : %.svg
	sed '/<text>.*<\/text>/d' $< | sed '/<rect x="0" y="0" width=".*" height=".*" fill="none"\/>/d' > $@ # 1st run: remove text and rec to determing BBox
	VAL=`inkscape --without-gui --query-all $@ | grep gnuplot_canvas | awk -F, '{print $$2, $$3, $$4, $$5}' ` ; \
	sed "s|viewBox=.*\$$|viewBox='$$VAL' preserveAspectRatio='xMinYMin slice'|" $< > $@ # adjust SVG viewbox for typesetting with latex ;-) # keep rec which is needed in "inc-shearY"-rule later

%-VB.svg : %.svg # remove viewBox for correct view in e.g. firefox
	sed "/viewBox.*preserveAspectRatio='xMinYMin slice'/d" $< > $@

.PRECIOUS: %_shearY.svg
%2Dshist_inc-shearY_VB.svg : %2Dshist_shearY.svg %2Danno_VB.svg
	sed 's/<rect x="0" y="0" width=".*" height=".*" fill="none"\/>/<use x="0" y="0" xlink:href="$<#gnuplot_canvas" \/>/' $(lastword $^) > $@ # <use ...> needs Inkscape 0.91, SVG inclusion with <image ...> is rasterized by Inkscape  http://superuser.com/questions/255086/is-it-possible-to-embed-or-link-one-inkscape-svg-document-inside-another-one  http://stackoverflow.com/questions/5451135/embed-svg-in-svg
#	inkscape --verb=FitCanvasToDrawing --verb=FileSave --verb=FileClose $@ # auto crop with inkscape: https://shkspr.mobi/blog/2013/03/inkscape-cropping-svg-files-on-the-command-line/ # using JS: http://stackoverflow.com/questions/23560038/html5-inline-svg-autocrop-whitespace#23698133

.PRECIOUS: %.efit-2DpdistOV.svg %_shearY.svg
%2Dshist_inc-pdist-shearY_VB.svg : %2DpdistOV.svg %2Dshist_shearY.svg %2Danno_VB.svg
	sed 's/<rect x="0" y="0" width=".*" height=".*" fill="none"\/>/<use x="0" y="0" xlink:href="$(word 1,$^)#gnuplot_canvas" \/><use x="0" y="0" style="opacity:0.90" xlink:href="$(word 2,$^)#gnuplot_canvas" \/>/' $(lastword $^) > $@


TEST:= 3Dasym.svg 2Dpdist.svg 3Dasym_VB.svg 2Dpdist_VB.svg 2Dshist_inc-shearY_VB.svg 2Dshist_inc-pdist-shearY_VB.svg 2Dshist_inc-shearY_VB-VB.svg 2Dshist_inc-pdist-shearY_VB-VB.svg
TESTf:= t00.efit t01.efit t02.efit t03.efit
TESTs:= $(foreach testf,$(TESTf),$(TEST:%=test/$(testf)-%))


test : $(TESTs)

test/% :
	$(MAKE) -C $(dir $@) -f ../Makefile ELLEVAL=../ $(notdir $@)
