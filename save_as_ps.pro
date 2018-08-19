PRO save_ps

  xsize=45
  ysize=45

  ; go_time_psa,000136
  ; set_scale_psa,2100,2600
  ; set_plot,'ps'
  ; device,filename='./TEST_PLOTS/allsky_00136.eps',xsize=xsize,ysize=ysize,/color,/encapsulated,bits=8,/cmyk
  ; ; dpanel_psa,1,1,0,0,/square
  ; plot_psa_raw,1,1,smo=10
  ; device,/close 
  ; set_plot, 'x'

  ; go_time_psa,000137
  ; set_scale_psa,2100,2600
  ; set_plot,'ps'
  ; device,filename='./TEST_PLOTS/allsky_00137.eps',xsize=xsize,ysize=ysize,/color,/encapsulated,bits=8,/cmyk
  ; ; dpanel_psa,1,1,0,0,/square
  ; plot_psa_raw,1,1,smo=10
  ; device,/close 
  ; set_plot, 'x'

  ; go_time_psa,000138
  ; set_scale_psa,2100,2600
  ; set_plot,'ps'
  ; device,filename='./TEST_PLOTS/allsky_00138.eps',xsize=xsize,ysize=ysize,/color,/encapsulated,bits=8,/cmyk
  ; ; dpanel_psa,1,1,0,0,/square
  ; plot_psa_raw,1,1,smo=10
  ; device,/close 
  ; set_plot, 'x'

  go_time_psa,000139
  set_scale_psa,2100,2600
  set_plot,'ps'
  device,filename='./TEST_PLOTS/allsky_00139.eps',xsize=xsize,ysize=ysize,/color,/encapsulated,bits=8,/cmyk
  ; dpanel_psa,1,1,0,0,/square
  plot_psa_raw,1,1,smo=10
  device,/close 
  set_plot, 'x'

  ; go_time_psa,000140
  ; set_scale_psa,2100,2600
  ; set_plot,'ps'
  ; device,filename='./TEST_PLOTS/allsky_00140.eps',xsize=xsize,ysize=ysize,/color,/encapsulated,bits=8,/cmyk
  ; ; dpanel_psa,1,1,0,0,/square
  ; plot_psa_raw,1,1,smo=10
  ; device,/close 
  ; set_plot, 'x'

  ; ; go_time_psa,000141
  ; set_scale_psa,2100,2600
  ; set_plot,'ps'
  ; device,filename='./TEST_PLOTS/colour_bar.eps',xsize=xsize,ysize=ysize,/color,/encapsulated,bits=8,/cmyk
  ; ; dpanel_psa,1,1,0,0,/square
  ; ; plot_psa_raw,1,1,smo=10
  ; plot_psa_colour_bar
  ; device,/close 
  ; set_plot, 'x'

END