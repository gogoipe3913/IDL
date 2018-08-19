;swith all-sky plot images 
PRO go_next_image, no_plot=no_plot
	COMMON raw_psa
	COMMON plugin_psa

	now_raw_next = now_raw
	now_raw_next += 100
	go_psa_raw, now_raw_next
	; if not keyword_set(no_plot) then plot_psa_raw, smo = 10
	if keyword_set(positionMark) and positionMark eq 'true' then oplot_linePosition
END

PRO go_back_image
	COMMON raw_psa
	COMMON plugin_psa
	now_raw_back = now_raw
	now_raw_back -= 100
	go_psa_raw, now_raw_back
	plot_psa_raw, smo = 10
	if keyword_set(positionMark) and positionMark eq 'true' then oplot_linePosition
END


;save all-sky plot as png file.
pro save_allSky_plot
	COMMON raw_psa
	hh_str = STRING(hh_raw(now_raw),FORMAT='(I2.2)')
	mm_str = STRING(mm_raw(now_raw),FORMAT='(I2.2)')
	ss_str = STRING(ss_raw(now_raw),FORMAT='(I2.2)')
	ms_str = STRING(ms_raw(now_raw),FORMAT='(I3.3)')
	current_hms = hh_str + mm_str + ' ' + ss_str + 's (' + ms_str + ' ms)'
	
	filename = "./ALL_SKY_SOD_COLOR_100Hz/" + day_raw[0] + ":" + current_hms + ".png"
	write_png, filename, TVRD(TRUE = 1)
	print, 'File written to ', filename
end