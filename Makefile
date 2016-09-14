
SHELL:= /bin/bash


%.efit_AR_06.svg : %.efit
	octave ell_type.m $<


%_shearY.svg : %.svg
	cp $< $@
	sed -i '/<rect x="0" y="0" width=".*" height=".*" fill="none"\/>/d' $@ # remove BG rec, not possible from within octave as gnuplot term is overwritten by print-command? (not tried) http://stackoverflow.com/questions/18169221/gnuplot-png-file-without-border-line
	## get BBox coords with inkscape ;-) https://sourceforge.net/p/inkscape/mailman/message/21202829/ # not available in rsvg-convert from librsvg nor svg_utils https://github.com/btel/svg_utils
	## SVG transform https://www.w3.org/TR/SVG/coords.html#TransformAttributeEffectOnSiblingAttributes
	x=`inkscape --without-gui --query-id=gnuplot_canvas -X $@ ` ; \
	y=`inkscape --without-gui --query-id=gnuplot_canvas -Y $@ ` ; \
	sed -i 's/<g id="gnuplot_canvas">/<g id="gnuplot_canvas" transform="translate('$$x,$$y') skewY(-30) translate('-$$x,-$$y')" >/' $@ # basically works, ToDo: do not plot text (overlay in latex) or also repos "separation curve", "ellipse arc", "oblate line", "circle\npoint" # https://developer.mozilla.org/de/docs/Web/SVG/Attribute/transform
#	inkscape --verb=FitCanvasToDrawing --verb=FileSave --verb=FileClose $*.efit_AR_21.svg # auto crop with inkscape (only makes sense without text): https://shkspr.mobi/blog/2013/03/inkscape-cropping-svg-files-on-the-command-line/ # using JS: http://stackoverflow.com/questions/23560038/html5-inline-svg-autocrop-whitespace#23698133
# ToDo: insert $*.efit_AR_21.svg into $*.efit_AR_20.svg # <use x="0" y="0" width="400" height="400" xlink:href="h3d_seg_Alv.efit_AR_21.svg#Mgnuplot_canvas" xlink:type="simple" xlink:actuate="onLoad" xlink:show="embed" /> http://superuser.com/questions/255086/is-it-possible-to-embed-or-link-one-inkscape-svg-document-inside-another-one # inkscape rastersizes <image ...>, not sure about pdf export: http://stackoverflow.com/questions/5451135/embed-svg-in-svg

%_inc-shearY.svg : %.svg
	cp $< $@
	sed -i 's/<rect x="0" y="0" width=".*" height=".*" fill="none"\/>/<use x="0" y="0" xlink:href="h3d_seg_Alv.efit_AR_06_shearY.svg#gnuplot_canvas" \/>/' $@ # <use ...> needs Inkscape 0.91, SVG inclusion with <image ...> is rasterized by Inkscape 
