COMMON plugin_psa, positionMark

PRO oplot_linePosition
	COMMON plugin_psa

	loadct, [39]
	; OPLOT, [24], [124], color=240, THICK=2, psym=7
	OPLOT, [128], [128], color=240, THICK=5, psym=7
	OPLOT, [123], [148], color=240, THICK=5, psym=7
	OPLOT, [118], [168], color=240, THICK=5, psym=7
	OPLOT, [113], [188], color=240, THICK=5, psym=7
	loadct, [0]
	;一度使ったら次からは自動でpositionマークをつける
	positionMark = 'true'
END