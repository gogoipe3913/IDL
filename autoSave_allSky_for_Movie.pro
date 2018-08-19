pro save_for_movie

	COMMON raw_psa

	startTime = ''
	read, startTime, PROMPT = 'Put start Time :' 
	read, duration, PROMPT = 'Put duration of minutes:' 
	
	go_time_psa,startTime
	duration = duration*60
	
	; 1枚目の保存
	plot_psa_raw,smo=10
	filename = "./ALLSKY_FOR_MOVIE/" + STRMID(ti_raw(now_raw),0,17) + '.png'
	write_png, filename, TVRD(TRUE = 1)
	print, 'File written to ', filename

	; 一枚目以降の保存
	for i=1, duration*10 -1 do begin	
		plot_psa_raw,smo=10,skip=10	
		filename = "./ALLSKY_FOR_MOVIE/" + STRMID(ti_raw(now_raw),0,17) + '.png'
		write_png, filename, TVRD(TRUE = 1)
		print, 'File written to ', filename

	endfor
end