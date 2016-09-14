
SHELL:= /bin/bash


%.efit-01_3Dsym.svg \
%.efit-02_3Dasym.svg \
%.efit-03_2Dpdist.svg \
%.efit-04_2Dhist.svg \
%.efit-05_2Danno.svg \
%.efit-06_2Dshist.svg \
 : %.efit
	octave ell_type.m $<


%_shearY.svg : %.svg
	sed '/<rect x="0" y="0" width=".*" height=".*" fill="none"\/>/d' $< > $@ # remove BG rec, not possible from within octave as gnuplot term is overwritten by print-command? (not tried) http://stackoverflow.com/questions/18169221/gnuplot-png-file-without-border-line
	## get BBox coords with inkscape ;-) https://sourceforge.net/p/inkscape/mailman/message/21202829/ # not available in rsvg-convert from librsvg nor svg_utils https://github.com/btel/svg_utils
	## SVG transform https://www.w3.org/TR/SVG/coords.html#TransformAttributeEffectOnSiblingAttributes
	x=`inkscape --without-gui --query-id=gnuplot_canvas -X $@ ` ; \
	y=`inkscape --without-gui --query-id=gnuplot_canvas -Y $@ ` ; \
	sed -i 's/<g id="gnuplot_canvas">/<g id="gnuplot_canvas" transform="translate('$$x,$$y') skewY(-30) translate('-$$x,-$$y')" >/' $@ # basically works, ToDo: do not plot text (overlay in latex) or also repos "separation curve", "ellipse arc", "oblate line", "circle\npoint" # https://developer.mozilla.org/de/docs/Web/SVG/Attribute/transform

%_2Danno_VB.svg : %_2Danno.svg
	sed '/<text>.*<\/text>/d' $< | sed '/<rect x="0" y="0" width=".*" height=".*" fill="none"\/>/d' > $@ # 1st run: remove text and rec to determing BBox
	VAL=`inkscape --without-gui --query-all $@ | grep gnuplot_canvas | awk -F, '{print $$2, $$3, $$4, $$5}' ` ; \
	sed "s|viewBox=.*\$$|viewBox='$$VAL' preserveAspectRatio='xMinYMin slice'|" $< > $@ # adjust SVG viewbox for typesetting with latex ;-) # keep rec which is needed in "inc-shearY"-rule later

.PRECIOUS: %_shearY.svg
%06_2Dshist_inc-shearY.svg : %06_2Dshist_shearY.svg %05_2Danno_VB.svg
	sed 's/<rect x="0" y="0" width=".*" height=".*" fill="none"\/>/<use x="0" y="0" xlink:href="$<#gnuplot_canvas" \/>/' $(lastword $^) > $@ # <use ...> needs Inkscape 0.91, SVG inclusion with <image ...> is rasterized by Inkscape  http://superuser.com/questions/255086/is-it-possible-to-embed-or-link-one-inkscape-svg-document-inside-another-one  http://stackoverflow.com/questions/5451135/embed-svg-in-svg
#	inkscape --verb=FitCanvasToDrawing --verb=FileSave --verb=FileClose $@ # auto crop with inkscape: https://shkspr.mobi/blog/2013/03/inkscape-cropping-svg-files-on-the-command-line/ # using JS: http://stackoverflow.com/questions/23560038/html5-inline-svg-autocrop-whitespace#23698133
