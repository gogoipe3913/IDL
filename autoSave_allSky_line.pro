pro save_as_line

	startDay = ''
	read, startDay, PROMPT = 'Put start day :' 
	read, duration, PROMPT = 'Put duration of minutes:' 
	
	
	for j=1, duration -1 do begin	
		file_psa_raw, 2, startDay, j, 1, /vbs
		time_psa, '00' + STRING(j, format = '(I02)') + '00', '00' + STRING(j + 1, format = '(I02)') + '00'
		for i=0, 59 do begin
			; print, startTime
			plot_as_line
			
			filename = "./ALLSKY_&_LINE_SOD_100Hz/" + '00' + STRING(j, format = '(I02)') + STRING(i, format = '(I02)') + '.png'
			write_png, filename, TVRD(TRUE = 1)
			print, 'File written to ', filename

			if i eq 59 then break
			go_next_image,no_plot=no_plot
		endfor
	endfor
end