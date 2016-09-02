

%.efit_AR_21.svg : %.efit
	octave ell_eval.m $<
	sed -i '/<rect x="0" y="0" width=".*" height=".*" fill="none"\/>/d' $*.efit_AR_21.svg # remove BG rec, not possible from within octave as gnuplot term is overwritten by print-command? (not tried) http://stackoverflow.com/questions/18169221/gnuplot-png-file-without-border-line
	sed -i 's/<g id="gnuplot_canvas">/<g id="gnuplot_canvas" transform="skewY(-30)">/' $*.efit_AR_21.svg # basically works, ToDo: do not plot text (overlay in latex) or also repos "separation curve", "ellipse arc", "oblate line", "circle\npoint" # https://developer.mozilla.org/de/docs/Web/SVG/Attribute/transform
	inkscape --verb=FitCanvasToDrawing --verb=FileSave --verb=FileClose $*.efit_AR_21.svg # auto crop with inkscape (only makes sense without text): https://shkspr.mobi/blog/2013/03/inkscape-cropping-svg-files-on-the-command-line/ # using JS: http://stackoverflow.com/questions/23560038/html5-inline-svg-autocrop-whitespace#23698133
# ToDo: insert $*.efit_AR_21.svg into $*.efit_AR_20.svg # inkscape rastersizes <image ...>, not sure about pdf export: http://stackoverflow.com/questions/5451135/embed-svg-in-svg

