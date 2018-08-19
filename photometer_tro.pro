;--------------------------------------------------------------------------------
;PHOTO_INIT.PRO
;This procedure was made by Taiki Kishiyama.
;
;This is procedure for reading data of photometer in spedas.
;Data of photometer is 20Hz.
;If you wanna use any other data of photometer, you must fix this procedure.
;
;<EXAMPLE>
;
;THEMIS>photo_load_count,0205,210000,inter=5,sec=120
;
;0205 is the day of data that you wanna read.
;210000 is the start time of analyzing.
;inter is interval of time. This example read 21 to 2.
;sec is smooth(background).
;
;--------------------------------------------------------------------------------

COMMON freq,        datearr1,timearr1,d1arr1,d2arr1,d3arr1, $
                    d4arr1,d5arr1
COMMON DAY_INFO,    day0,time0,time_tplot
COMMON SMOOTH,      smd1,smd2,smd3,smd4,smd5    
COMMON FFT,         fft_data,fft_frq
COMMON NORMAL,      nor_d1,nor_d2,nor_d3,nor_d4,nor_d5
COMMON CAPTURE,     st_num,en_num
COMMON GO_NEXT,     start_2min,end_2min,set_for_go_next,now_date
COMMON HERI_COUNT,  time1,time2,time3,time4,time5,time6,time7,$
                    count1,count2,count3,count4,count5,count6,count7,$
                    abs_time1,abs_time2,abs_time3,abs_time4,abs_time5,abs_time6,abs_time7,$
                    abs_count1,abs_count2,abs_count3,abs_count4,abs_count5,abs_count6,abs_count7,$
                    time_heri,count_heri,abTi_heri,abCt_heri
COMMON HERI_PLOT,   bg_info,per_info

pro set_sdate0,sdate0
  sdate0='20170205'
end

pro init_par1,date,time1,d1,d2,d3,d4,d5
  COMMON freq
  date=0l
  time1=0d+0
  d1=0l
  d2=0l
  d3=0l
  d4=0l
  d5=0l

end

pro make_arrs,num_photo=num_photo
  COMMON freq
  datearr1=lonarr(num_photo)
  timearr1=dblarr(num_photo)
  d1arr1= dblarr(num_photo)
  d2arr1= d1arr1
  d3arr1= d1arr1
  d4arr1= d1arr1
  d5arr1= d1arr1
  time_tplot=strarr(num_photo)

end                         

pro open_files,day,time,inter=inter
;EXAMPLE open_files,0206,010000,inter=5
  COMMON freq
  COMMON DAY_INFO
  day0=day
  time0=time
  time_end=time+(inter-1)*10000L
  fday=0l
  fday=day
  c='directy '+'is'+' '

  IF time/100000 LT 1 THEN BEGIN
    IF day mod 100 EQ 1 THEN BEGIN
      IF day/100 GE 1  AND day/100 LT 2  THEN BEGIN ;Jan
        fday=day+1130                                  
      ENDIF
         
      IF day/100 GE 2  AND day/100 LT 3  THEN BEGIN ;Feb
        fday=day-70        
      ENDIF
         
      IF day/100 GE 3  AND day/100 LT 4  THEN BEGIN ;Mar
        fday=day-73  
      ENDIF
       
      IF day/100 EQ 4  AND day/100 LT 5  THEN BEGIN ;Apl
        fday=day-70
      ENDIF
       
      IF day/100 EQ 5  AND day/100 LT 6  THEN BEGIN ;May
        fday=day-71  
      ENDIF
       
      IF day/100 EQ 6  AND day/100 LT 7  THEN BEGIN ;Jun
        fday=day-70   
      ENDIF
       
      IF day/100 EQ 7  AND day/100 LT 8  THEN BEGIN ;Jul
        fday=day-71   
      ENDIF
       
      IF day/100 EQ 8  AND day/100 LT 9  THEN BEGIN ;Aug
        fday=day-70   
      ENDIF
          
      IF day/100 EQ 9  AND day/100 LT 10 THEN BEGIN ;Sep
        fday=day-70    
      ENDIF
       
      IF day/100 EQ 10 AND day/100 LT 11 THEN BEGIN ;Oct
        fday=day-71  
      ENDIF
                         
      IF day/100 EQ 11 AND day/100 LT 12 THEN BEGIN ;Nov
        fday=day-70  
      ENDIF
                           
      IF day/100 EQ 12 AND day/100 LT 13 THEN BEGIN ;Dec
        fday=day-71   
      ENDIF
                           
    ENDIF ELSE BEGIN
      fday=day-1
    ENDELSE
  ENDIF

  sdate0='20170'+STRING(fday,FORMAT='(I0)')  
  ffday=STRING(day,FORMAT='(I4.4)')
  ftime=STRING(time,FORMAT='(I6.6)')
  ftime=STRMID(ftime,0,2)
  print,ftime
  file00='/raid/data/Kiban-S/ptmt/'+sdate0+'/ascii20*.dat'
  file_start='/raid/data/Kiban-S/ptmt/'+sdate0+'/ascii2017'+ffday+ftime+'*.dat'
  file_start=file_search(file_start)
  print,time_end
   
  IF time/100000 GE 1 THEN BEGIN
    IF time_end GE 240000 THEN BEGIN
      day+=1
        IF day EQ 132  THEN day=day-1+70    ;Jan
        IF day EQ 229  THEN day=day-1+73    ;Feb
        IF day EQ 332  THEN day=day-1+70    ;Mar
        IF day EQ 431  THEN day=day-1+71    ;Apl
        IF day EQ 532  THEN day=day-1+70    ;May
        IF day EQ 631  THEN day=day-1+71    ;Jun
        IF day EQ 732  THEN day=day-1+70    ;Jul
        IF day EQ 832  THEN day=day-1+70    ;Aug
        IF day EQ 931  THEN day=day-1+71    ;Sep
        IF day EQ 1032 THEN day=day-1+70    ;Oct
        IF day EQ 1131 THEN day=day-1+71    ;Nov 
        IF day EQ 1232 THEN day=day-1-1130  ;Dec
        time_end=time_end-240000
      ENDIF
    ENDIF

  IF time/100000 LT 1 THEN BEGIN
    IF time_end GE 240000 THEN BEGIN
      day=day
      time_end=time_end-240000
    ENDIF
  ENDIF

  day=STRING(day,FORMAT='(I4.4)')
  time2=STRING(time_end,FORMAT='(I6.6)')
  time2=STRMID(time2,0,2)

  file_end='/raid/data/Kiban-S/ptmt/'+sdate0+'/ascii2017'+day+time2+'*.dat'
  file_end=file_search(file_end)

  dfilenames=file_search(file00,count=dfiles)
  iloopend=n_elements(dfilenames)

  print,dfilenames
  num_search1=0l
  num_search2=0l
  num_search1=where(dfilenames EQ file_start(0))
  num_search2=where(dfilenames EQ file_end(0))

  num_photo=72000l*(num_search2(0)-num_search1(0)+1)
  make_arrs,num_photo=num_photo
  num_photo=0l

  FOR i=num_search1(0),num_search2(0) DO BEGIN
    j=1
    IF i EQ num_search1(0) THEN j=0

    file=dfilenames(i)
    openr,1,file
    print,'--- Open(r) = ',file
    moji="moji"
    num_photo=0l
    readf,1,format='(I8)',num_photo

    IF j EQ 0 THEN BEGIN
      back_num_photo=0l
      read_files,0,num_photo=num_photo,back_num_photo=back_num_photo
      back_num_photo=num_photo-1
    ENDIF
    IF j EQ 1 THEN BEGIN
      read_files,1,num_photo=num_photo,back_num_photo=back_num_photo
      back_num_photo+=num_photo
    ENDIF

    close,1 
  ENDFOR
  close,1
  print,file_start
  print,file_end
END

PRO read_files,key,num_photo=num_photo,back_num_photo=back_num_photo
  COMMON freq

  init_par1,date,time1,d1,d2,d3,d4,d5
  str0='(I7,F10.3,5I6)'
  i=0l
  IF key EQ 0 THEN BEGIN
    FOR i=0,num_photo-1 DO BEGIN ;20Hz(0.05s)*60sec*60min=72000
      readf,1,format=str0,date,time1,d1,d2,d3,d4,d5
      print,date,time1,d1,d2,d3,d4,d5
      zerodata=(2l)^15
      datearr1(i)=date
      timearr1(i)=time1
      d1arr1(i)=d1-zerodata
      d2arr1(i)=d2-zerodata
      d3arr1(i)=d3-zerodata
      d4arr1(i)=d4-zerodata
      d5arr1(i)=d5-zerodata       
    ENDFOR
  ENDIF
       
  IF key EQ 1 THEN BEGIN
    FOR i=back_num_photo+1,(back_num_photo+num_photo) DO BEGIN ;20Hz(0.05s)*60sec*60min=72000
      readf,1,format=str0,date,time1,d1,d2,d3,d4,d5
      print,date,time1,d1,d2,d3,d4,d5
      zerodata=(2l)^15
      datearr1(i)=date
      timearr1(i)=time1
      d1arr1(i)=d1-zerodata
      d2arr1(i)=d2-zerodata
      d3arr1(i)=d3-zerodata
      d4arr1(i)=d4-zerodata
      d5arr1(i)=d5-zerodata       
    ENDFOR
  ENDIF


END

pro mkarr_for_tplot

  COMMON freq
  COMMON DAY_INFO       

  time_tplot=STRARR(n_elements(datearr1))
  datearr2=STRARR(n_elements(datearr1))
  timearr2=DBLARR(n_elements(timearr1))
  timearr3=STRARR(n_elements(timearr1))
  decimal1=LONARR(n_elements(timearr1))
  decimal2=STRARR(n_elements(timearr1))

;----------time change from through sec to hhmmss------------
  FOR k=0,N_ELEMENTS(timearr1)-1 DO BEGIN
    timearr2(k)=LONG(timearr1(k)/3600)*10000+$
    LONG((timearr1(k)mod 3600)/60)*100+$
    (timearr1(k) mod 3600) mod 60
    decimal1(k)=(((timearr1(k) mod 3600) mod 60) mod 1)*100
  ENDFOR

;---------------date & time change for STRING----------------
  FOR j=0,n_elements(datearr1)-1 DO BEGIN
    datearr2(j)=STRING(datearr1(j))
    timearr3(j)=STRING(timearr2(j),format='(I6.6)')
    decimal2(j)=STRING(decimal1(j),format='(I02)')
  ENDFOR
     
  FOR i=0,n_elements(datearr1)-2 DO BEGIN      
    ;this_year=FLOAT(STRMID(datearr2(i),0,4))
    Judge=datearr1(i)-2017000
;-------------------------make month-------------------------
    MONTH='moji'
    DAY='moji'
    IF Judge LE 31 THEN BEGIN
      MONTH='01'
      Judge0=Judge
      DAY=STRING(Judge,format='(I0)')
    ENDIF
    IF Judge GE 32  && Judge LE 59  THEN BEGIN
      MONTH ='02'
      Judge0=Judge-31
      DAY=STRING(Judge0,format='(I0)')
    ENDIF
    IF Judge GE 60  && Judge LE 90  THEN BEGIN
      MONTH ='03'
      Judge0=Judge-59
      DAY=STRING(Judge0,format='(I0)')
    ENDIF   
    IF Judge GE 91  && Judge LE 120 THEN BEGIN
      MONTH ='04'
      Judge0=Judge-90
      DAY=STRING(Judge0,format='(I0)')
    ENDIF
    IF Judge GE 121 && Judge LE 151 THEN BEGIN
      MONTH ='05'
      Judge0=Judge-120
      DAY=STRING(Judge0,format='(I0)')
    ENDIF
    IF Judge GE 152 && Judge LE 181 THEN BEGIN
      MONTH ='06'
      Judge0=Judge-151
      DAY=STRING(Judge0,format='(I0)')
    ENDIF
    IF Judge GE 182 && Judge LE 212 THEN BEGIN
      MONTH ='07'
      Judge0=Judge-181
      DAY=STRING(Judge0,format='(I0)')
    ENDIF
    IF Judge GE 213 && Judge LE 243 THEN BEGIN
      MONTH ='08'
      Judge0=Judge-212
      DAY=STRING(Judge0,format='(I0)')
    ENDIF
    IF Judge GE 244 && Judge LE 273 THEN BEGIN
      MONTH ='09'
      Judge0=Judge-243
      DAY=STRING(Judge0,format='(I0)')
    ENDIF
    IF Judge GE 274 && Judge LE 304 THEN BEGIN
      MONTH ='10'
      Judge0=Judge-273
      DAY=STRING(Judge0,format='(I0)')
    ENDIF
    IF Judge GE 305 && Judge LE 334 THEN BEGIN
      MONTH ='11'
      Judge0=Judge-304
      DAY=STRING(Judge0,format='(I0)')
    ENDIF
    IF Judge GE 335 && Judge LE 365 THEN BEGIN
      MONTH ='12'
      Judge0=Judge-334
      DAY=STRING(Judge0,format='(I0)')
    ENDIF
;-----------------Correct degit------------------------------            
    IF Judge0 LT 10 THEN DAY='0'+DAY
;--------------------make time array for tplot---------------        
    time_tplot(i)='2017-'+MONTH+'-'+DAY+'/'+$
                STRMID(timearr3(i),0,2)+':'+STRMID(timearr3(i),2,2)+':'+$
                STRMID(timearr3(i),4,2)+'.'+decimal2(i)
  ENDFOR

END

PRO smooth_go,sec

  COMMON freq
  COMMON DAY_INFO
  COMMON SMOOTH

  smd1=smooth(d1arr1,8)-smooth(smooth(d1arr1,8),sec*20)
  smd2=smooth(d2arr1,15)-smooth(smooth(d2arr1,15),sec*20)
  smd3=smooth(d3arr1,8)-smooth(smooth(d3arr1,8),sec*20)
  smd4=smooth(d4arr1,15)-smooth(smooth(d4arr1,15),sec*20)
  smd5=smooth(d5arr1,15)-smooth(smooth(d5arr1,15),sec*20)


END   

PRO ptmt_load_plot
;Read photometer data for spedas without smooth
  COMMON freq
  COMMON DAY_INFO
  COMMON SMOOTH

  store_data,'ptmt1',data={x:time_double(time_tplot),y:d1arr1},dlimits=dl,lim=lim
  store_data,'ptmt2',data={x:time_double(time_tplot),y:d2arr1},dlimits=dl,lim=lim 
  store_data,'ptmt3',data={x:time_double(time_tplot),y:d3arr1},dlimits=dl,lim=lim
  store_data,'ptmt4',data={x:time_double(time_tplot),y:d4arr1},dlimits=dl,lim=lim
  store_data,'ptmt5',data={x:time_double(time_tplot),y:d5arr1},dlimits=dl,lim=lim

END

FUNCTION normalizer,array         
  ret_array=DBLARR(N_ELEMENTS(array))
  ret_array=(array-min(array))/(max(array)-min(array))
  RETURN,ret_array
END

PRO photo_sm,sec
;Read photometer data for spedas with smooth
;Using keyword 'sec', you can process background smooth
  COMMON freq
  COMMON DAY_INFO
  COMMON SMOOTH
  COMMON NORMAL
       
  IF NOT KEYWORD_SET(sec) THEN BEGIN
    ns_back_d1=smooth(d1arr1,8)
    ns_back_d2=smooth(d2arr1,10)
    ns_back_d3=smooth(d3arr1,8)
    ns_back_d4=smooth(d4arr1,10)
    ns_back_d5=smooth(d5arr1,15)

    cut_nbg_d1=ns_back_d1(0:N_ELEMENTS(ns_back_d1)-2)
    cut_nbg_d2=ns_back_d2(0:N_ELEMENTS(ns_back_d2)-2)
    cut_nbg_d3=ns_back_d3(0:N_ELEMENTS(ns_back_d3)-2)
    cut_nbg_d4=ns_back_d4(0:N_ELEMENTS(ns_back_d4)-2)
    cut_nbg_d5=ns_back_d5(0:N_ELEMENTS(ns_back_d5)-2)

    ;norm_d1=DBLARR(N_ELEMENTS(cut_nbg_d1))
    ;norm_d2=DBLARR(N_ELEMENTS(cut_nbg_d2))
    ;norm_d3=DBLARR(N_ELEMENTS(cut_nbg_d3))
    ;norm_d4=DBLARR(N_ELEMENTS(cut_nbg_d4))
    ;norm_d5=DBLARR(N_ELEMENTS(cut_nbg_d5))
      
    ;norm_d1=normalizer(cut_nbg_d1)
    ;norm_d2=normalizer(cut_nbg_d2)
    ;norm_d3=normalizer(cut_nbg_d3)
    ;norm_d4=normalizer(cut_nbg_d4)
    ;norm_d5=normalizer(cut_nbg_d5)

    nor_d1=cut_nbg_d1
    nor_d2=cut_nbg_d2
    nor_d3=cut_nbg_d3
    nor_d4=cut_nbg_d4
    nor_d5=cut_nbg_d5
    
    store_data,'ptmt1',data={x:time_double(time_tplot(0:N_ELEMENTS(time_tplot)-2)),y:cut_nbg_d1},dlimits=dl,lim=lim
    store_data,'ptmt2',data={x:time_double(time_tplot(0:N_ELEMENTS(time_tplot)-2)),y:cut_nbg_d2},dlimits=dl,lim=lim 
    store_data,'ptmt3',data={x:time_double(time_tplot(0:N_ELEMENTS(time_tplot)-2)),y:cut_nbg_d3},dlimits=dl,lim=lim
    store_data,'ptmt4',data={x:time_double(time_tplot(0:N_ELEMENTS(time_tplot)-2)),y:cut_nbg_d4},dlimits=dl,lim=lim
    store_data,'ptmt5',data={x:time_double(time_tplot(0:N_ELEMENTS(time_tplot)-2)),y:cut_nbg_d5},dlimits=dl,lim=lim
  ENDIF
  IF KEYWORD_SET(sec) THEN BEGIN
    smooth_go,sec
    nor_d1=smd1(0:N_ELEMENTS(smd1)-2)
    nor_d2=smd2(0:N_ELEMENTS(smd2)-2)
    nor_d3=smd3(0:N_ELEMENTS(smd3)-2)
    nor_d4=smd4(0:N_ELEMENTS(smd4)-2)
    nor_d5=smd5(0:N_ELEMENTS(smd5)-2)

    store_data,'pmsm1',data={x:time_double(time_tplot(0:N_ELEMENTS(time_tplot)-2)),y:nor_d1},dlimits=dl,lim=lim
    store_data,'pmsm2',data={x:time_double(time_tplot(0:N_ELEMENTS(time_tplot)-2)),y:nor_d2},dlimits=dl,lim=lim 
    store_data,'pmsm3',data={x:time_double(time_tplot(0:N_ELEMENTS(time_tplot)-2)),y:nor_d3},dlimits=dl,lim=lim
    store_data,'pmsm4',data={x:time_double(time_tplot(0:N_ELEMENTS(time_tplot)-2)),y:nor_d4},dlimits=dl,lim=lim
    store_data,'pmsm5',data={x:time_double(time_tplot(0:N_ELEMENTS(time_tplot)-2)),y:nor_d5},dlimits=dl,lim=lim
  ENDIF
      
END

PRO photo_load_count,day,time,sec=sec
;Store data photometer for spedas   
  COMMON freq
  COMMON DAY_INFO
  COMMON SMOOTH

  READ,inter,PROMPT='How long hour ? : '
  open_files,day,time,inter=inter
  mkarr_for_tplot
  timespan,time_tplot(0),inter,/h
  IF NOT KEYWORD_SET(sec) THEN photo_sm
  IF KEYWORD_SET(sec) THEN photo_sm,sec
  message='MOJI'
  message='Default is no smooth'
  print,message

END

PRO plot_data,time,time2,ch=ch,no_x=no_x,no_y=no_y
;EXAMPLE plot_data1,010000,013000,ch=1
  COMMON freq
  COMMON DAY_INFO

  IF KEYWORD_SET(position) THEN !P.POSITION=position
;------set timerange for Plot-------
    stime=time                         ;始まりの時間
    etime=time2                        ;終わりの時間

    stime_n=FLOAT(stime)               ;Change to FLOAT type
    etime_n=FLOAT(etime)
    ;PRINT,stime_n,etime_n

    sh=FIX(stime / 10000)              ;start hour
    sm=FIX((stime_n mod 10000)/100)    ;start min
    ss=stime_n mod 100                 ;start sec

    eh=FIX(etime_n / 10000)            ;end hour
    em=FIX((etime_n mod 10000)/100)    ;end min
    es=etime_n mod 100                 ;end sec

    d_hour=eh-sh

    IF sh EQ 24 AND eh EQ 24 THEN BEGIN
    sm=sm+60
    em=em+60
  ENDIF

  IF sm LT 1 THEN BEGIN
    start_time=ss
  ENDIF ELSE BEGIN  
    start_time=60*sm+ss                ;min&sec change to sec 
  ENDELSE

  IF d_hour GT 0 THEN BEGIN
    em=em+60 ;時間をまたぐ仕組み
  ENDIF

  IF em eq 0 and es eq 0 THEN BEGIN
    end_time=3600
  ENDIF ELSE BEGIN
    end_time=60*em+es
  ENDELSE

  mass_trim=(end_time-start_time)*20
  num_stime=start_time*20
  num_etime=end_time*20

  xtitle='UT'
  ytitle='COUNT'
  IF ch EQ 1 THEN BEGIN
    plotarr1_smo=smooth(d1arr1,8);-smooth(d1arr1,3000)
    Freq_info='427.8 nm'
  ENDIF
  IF ch EQ 2 THEN BEGIN
    plotarr1_smo=smooth(d2arr1,10)
    Freq_info='557.7 nm'
  ENDIF 
  IF ch EQ 3 THEN BEGIN
    plotarr1_smo=smooth(d3arr1,10)
    Freq_info='670.0 w'
  ENDIF
  IF ch EQ 4 THEN BEGIN
    plotarr1_smo=smooth(d4arr1,30)
    Freq_info='777.9 nm'
  ENDIF
  IF ch EQ 5 THEN BEGIN
    plotarr1_smo=smooth(d5arr1,40)
    Freq_info='844.6 nm'
  ENDIF 

  hh_stime=time/10000L
  mm_stime=(time-hh_stime*10000L)/100L
  ss_stime=time-hh_stime*10000L-mm_stime*100L
  stm_photo=LONG(3600.0*hh_stime+60.0*mm_stime+ss_stime)

  hh_etime=time2/10000L
  mm_etime=(time2-hh_etime*10000L)/100L
  ss_etime=time2-hh_etime*10000L-mm_etime*100L
  etm_photo=LONG(3600.0*hh_etime+60.0*mm_etime+ss_etime)
  ts=stm_photo/3600.0 & te=etm_photo/3600.0

  time_axis_intervals,ts,te,minor_ticks,where_tick,no_of_ticks

  search_st=where(LONG(timearr1) EQ stm_photo)
  print,'stm_photo',stm_photo
  start_num=min(search_st)

  search_ed=where(LONG(timearr1) EQ etm_photo)
  print,'etm_photo',etm_photo
  end_num=min(search_ed)
     
  IF KEYWORD_SET(no_x) THEN BEGIN
    xtit=''
    xtickformat='no_ticks'
    IF no_x EQ 1 THEN ytickname=[' '] ELSE ytickname=['']
  ENDIF ELSE BEGIN
    xtit='UT'
    xtickformat='axis_time'
    ytickname=['']
  ENDELSE
  IF KEYWORD_SET(no_y) THEN BEGIN
    ytit=''
    ytickformat='no_ticks'
  ENDIF
        
  tim_plot=timearr1[start_num :end_num]/3600.0
  freq_array=plotarr1_smo[start_num :end_num]
  freq_range=[min(freq_array),max(freq_array)+100]
     
  print,'tim_plot',N_ELEMENTS(tim_plot)
  print,'freq_array',N_ELEMENTS(freq_array)
     
  RED='0F00F0'XL
  PLOT,[0],[0],title=title,xrange=[ts,te],yrange=freq_range, $
    xtitle=xtitle,ytitle=ytitle,position=position, $
    XTICKFORMAT=xtickformat,YTICKNAME=ytickname,YTICKS=yticks,YMINOR=yminor,$
    XTICKS=no_of_ticks,XTICKV=where_tick,XMINOR=minor_ticks,$
    YTICKFORMAT=ytickformat,CHARSIZE=!P.CHARSIZE*1.5
  OPLOT,tim_plot,freq_array,color=color,LINESTYLE=0,THICK=2
;------------time--info--------------

  Timim1=STRING(time,FORMAT='(I6.6)')
  Timim2=STRING(time2,FORMAT='(I6.6)')
  HH1=STRMID(Timim1,0,2)
  HH2=STRMID(Timim2,0,2)
  MM1=STRMID(Timim1,2,2)
  MM2=STRMID(Timim2,2,2)
  SS1=STRMID(Timim1,4,2)
  SS2=STRMID(Timim2,4,2)

  YDT='2017-'+'0'+STRING(day0,FORMAT='(I0)')+' '$
     +HH1+':'+MM1+':'+SS1+' to '+HH2+':'+MM2+':'+SS2+' UT'
  time_info=YDT
  XYOUTS,!P.POSITION(0)+0.15,!P.POSITION(3)+0.004,time_info,/NORMAL,ALI=1.0
  XYOUTS,!P.POSITION(2),!P.POSITION(3)+0.004,Freq_info,/NORMAL,ALI=1.0
END

PRO heri_count_photo
  COMMON freq
  COMMON DAY_INFO
  COMMON NORMAL
  COMMON CAPTURE
  COMMON GO_NEXT

  READ,heri_type,PROMPT='Which is decreasing types? : '
  
  IF heri_type EQ 0 THEN BEGIN               ;After ON TIME
     openu,/APPEND,1,'../heri_count2/count_after.txt'
  ENDIF

  IF heri_type EQ 1 THEN BEGIN               ;Before ON TIME
     openu,/APPEND,1,'../heri_count2/count_before.txt'
  ENDIF

  psym=1
  fin_time=DBLARR(7)
  luminacity=DBLARR(7)
  tim_sub=DBLARR(7)
  lumi_sub=DBLARR(7)
  ;FOR j=0,6 DO BEGIN
  ctime,fin_time,luminacity,npoints=7,prompt="Select point of plot",$
       hours=hours,minutes=minutes,seconds=seconds
  tim_sub(6)=0
  lumi_sub(6)=0

  tim_sub(0)=fin_time(1)-fin_time(0)
  tim_sub(1)=fin_time(2)-fin_time(0)
  tim_sub(2)=fin_time(3)-fin_time(0)
  tim_sub(3)=fin_time(4)-fin_time(0)
  tim_sub(4)=fin_time(5)-fin_time(0)
  tim_sub(5)=fin_time(6)-fin_time(0)

  lumi_sub(0)=luminacity(1)-luminacity(0)
  lumi_sub(1)=luminacity(2)-luminacity(0)
  lumi_sub(2)=luminacity(3)-luminacity(0)
  lumi_sub(3)=luminacity(4)-luminacity(0)
  lumi_sub(4)=luminacity(5)-luminacity(0)
  lumi_sub(5)=luminacity(6)-luminacity(0)

  print,luminacity
  print,time_string(fin_time)
  print,tim_sub
  print,lumi_sub
         
  FOR i=0,6 DO BEGIN
    printf,1,fin_time(i),luminacity(i),tim_sub(i),lumi_sub(i)
  ENDFOR

  READ,numberThis,PROMPT='What number is this png in this time ??? : '
  numberThis=STRING(numberThis,FORMAT='(I01)')
  makepng,'heri_line'+now_date+'_'+start_2min+'_'+end_2min+'_'+numberThis
  
  print,'heri_line'+now_date+'_'+start_2min+'_'+end_2min+'_'+numberThis
  print,''
  print,'             書き込み終わったよ~!!'         
  print,''
  
  close,1
END

PRO SET_GO_NEXT,nnnext_time

  COMMON freq
  COMMON DAY_INFO
  COMMON GO_NEXT
  COMMON NORMAL
  COMMON CAPTURE
  ;nnnext_time is start time of your analyzing       
  IF NOT KEYWORD_SET(nnnext_time) THEN BEGIN
    print,''
    print,'PLEASE SET START TIME!'
    print,''
  ENDIF

  set_for_go_next=nnnext_time
  chr_go_time=STRING(set_for_go_next,FORMAT='(I06)')
  next_time_after=nnnext_time+200
  chr_go_after=STRING(next_time_after,FORMAT='(I06)')

  today_is=STRMID(time_tplot(0),0,11)
  IF set_for_go_next/100000 LT 1 THEN today_is=STRMID(time_tplot(N_ELEMENTS(time_tplot)-3),0,11)
  start_2min=STRMID(chr_go_time,0,2)+':'+STRMID(chr_go_time,2,2)+':'+STRMID(chr_go_time,4,2)
  end_2min=STRMID(chr_go_after,0,2)+':'+STRMID(chr_go_after,2,2)+':'+STRMID(chr_go_after,4,2)

  print,'NOW, TIME = '+start_2min+' to '+end_2min
  tlimit,today_is+start_2min,today_is+end_2min

  search_s2min=min(where(STRMID(time_tplot,0,19) EQ today_is+start_2min))
  search_e2min=min(where(STRMID(time_tplot,0,19) EQ today_is+end_2min))
  print,search_s2min,search_e2min

  max_lumi=max(nor_d1[search_s2min :search_e2min])
  min_lumi=min(nor_d1[search_s2min :search_e2min])
  print,'ylim '+STRING(min_lumi-200)+' to '+STRING(max_lumi+200)
  ylim,1,min_lumi-200,max_lumi+200
  tplot,['ptmt1']
  st_num=search_s2min
  en_num=search_e2min
  now_date=STRMID(today_is,5,2)+STRMID(today_is,8,2)
END

PRO GO_NEXT_2MIN

  COMMON freq
  COMMON DAY_INFO
  COMMON go_next
  COMMON NORMAL
  COMMON CAPTURE

  set_for_go_next+=200L
  IF (set_for_go_next mod 10000) GE 6000 THEN set_for_go_next=set_for_go_next-6000+10000
  IF set_for_go_next GE 240000 THEN set_for_go_next=000000

  next_time_after=set_for_go_next+200L
  IF (next_time_after mod 10000) GE 6000 THEN next_time_after=next_time_after-6000+10000
  IF next_time_after GT 240000 THEN next_time_after=000000

  chr_go_time=STRING(set_for_go_next,FORMAT='(I06)')
  chr_go_after=STRING(next_time_after,FORMAT='(I06)')

  print,set_for_go_next
  print,chr_go_after

  today_is=STRMID(time_tplot(0),0,11)
  IF set_for_go_next/100000 LT 1 THEN BEGIN
    today_is=STRMID(time_tplot(N_ELEMENTS(time_tplot)-71999),0,11)   
  ENDIF

  start_2min=STRMID(chr_go_time,0,2)+':'+STRMID(chr_go_time,2,2)+':'+STRMID(chr_go_time,4,2)
  end_2min=STRMID(chr_go_after,0,2)+':'+STRMID(chr_go_after,2,2)+':'+STRMID(chr_go_after,4,2)

  print,'----------------------------------------'
  print,'NOW, DAY  = '+today_is
  print,'NOW, TIME = '+start_2min+' to '+end_2min
  print,'----------------------------------------'

  tlimit,today_is+start_2min,today_is+end_2min

  search_s2min=min(where(STRMID(time_tplot,0,19) EQ today_is+start_2min))
  search_e2min=min(where(STRMID(time_tplot,0,19) EQ today_is+end_2min))
  print,search_s2min,search_e2min

  max_lumi=max(nor_d1[search_s2min :search_e2min])
  min_lumi=min(nor_d1[search_s2min :search_e2min])
  print,'ylim '+STRING(min_lumi-200)+' to '+STRING(max_lumi+200)
  ylim,1,min_lumi-200,max_lumi+200
  tplot,['ptmt1']

  st_num=search_s2min
  en_num=search_e2min

  now_date=STRMID(today_is,5,2)+STRMID(today_is,8,2)
END

PRO MAKE_WAVE_RATIO

  COMMON freq
  COMMON DAY_INFO
  COMMON NORMAL
        
  ns_back_d1=smooth(d1arr1,8)
  ns_back_d2=smooth(d2arr1,8)
  ns_back_d3=smooth(d3arr1,8)
  ns_back_d4=smooth(d4arr1,8)
  ns_back_d5=smooth(d5arr1,15)

  cut_nbg_d1=ns_back_d1(0:N_ELEMENTS(ns_back_d1)-2)
  cut_nbg_d2=ns_back_d2(0:N_ELEMENTS(ns_back_d2)-2)
  cut_nbg_d3=ns_back_d3(0:N_ELEMENTS(ns_back_d3)-2)
  cut_nbg_d4=ns_back_d4(0:N_ELEMENTS(ns_back_d4)-2)
  cut_nbg_d5=ns_back_d5(0:N_ELEMENTS(ns_back_d5)-2)

  RATIO_5577=cut_nbg_d2/cut_nbg_d1
  RATIO_670w=cut_nbg_d3/cut_nbg_d1
  RATIO_7779=cut_nbg_d4/cut_nbg_d1
  RATIO_8446=cut_nbg_d5/cut_nbg_d1


  nor_4278=normalizer(cut_nbg_d1)
  nor_rat_5577=normalizer(cut_nbg_d2)
  nor_rat_670w=normalizer(cut_nbg_d3)
  nor_rat_7779=normalizer(cut_nbg_d4)
  nor_rat_8446=normalizer(cut_nbg_d5)
       
  store_data,'ptmt1_normal',data={x:time_double(time_tplot[0:N_ELEMENTS(time_tplot)-2]),y:nor_4278}

  options,'ptmt1_normal','ytitle','Count'
  store_data,'pt_rat_557.7',data={x:time_double(time_tplot[0:N_ELEMENTS(time_tplot)-2]),y:nor_rat_5577},$
            dlim={color:6}
  ;options,'pt_rat_557.7',axis={yaxis:1, yrange:[950, 2000], ystyle:1}
         
  store_data,'pt_rat_670w',data={x:time_double(time_tplot[0:N_ELEMENTS(time_tplot)-2]),y:nor_rat_670w},$
            dlim={color:6}
  ;options,'pt_rat_670w',axis={yaxis:1, yrange:[1000, 1300], ystyle:1}
  options,'pt_rat_670w','ytitle','Count'

  store_data,'pt_rat_777.9',data={x:time_double(time_tplot[0:N_ELEMENTS(time_tplot)-2]),y:nor_rat_7779},$
            dlim={color:6}
  ;options,'pt_rat_777.9',axis={yaxis:1, yrange:[1000, 1200], ystyle:1}
  options,'pt_rat_777.9','ytitle','Count'
  store_data,'pt_rat_844.6',data={x:time_double(time_tplot[0:N_ELEMENTS(time_tplot)-2]),y:nor_rat_8446},$
            dlim={color:6}
  ;options,'pt_rat_844.6',axis={yaxis:1, yrange:[0, 250], ystyle:1}
  options,'pt_rat_844.6','ytitle','Count'

  options,'ptmt1_normal','labels',['427.8 nm']
  options,'pt_rat_557.7','labels',['557.7 nm']
  options,'pt_rat_670w','labels',['670w nm']
  options,'pt_rat_777.9','labels',['777.9 nm']
  options,'pt_rat_844.6','labels',['844.6 nm']

  store_data,'557.7',data=['ptmt1_normal','pt_rat_557.7']
  store_data,'670w',data=['ptmt1_normal','pt_rat_670w']
  store_data,'777.9',data=['ptmt1_normal','pt_rat_777.9']
  store_data,'844.6',data=['ptmt1_normal','pt_rat_844.6']

  options,'557.7','ytitle','Count'
  options,'670w','ytitle','Count'
  options,'777.9','ytitle','Count'
  options,'844.6','ytitle','Count'
  options,'557.7','labels',['','']
  options,'670w','labels',['','']
  options,'777.9','labels',['','']
  options,'844.6','labels',['','' ]
END

PRO compair_heri

  COMMON freq
  COMMON DAY_INFO
  COMMON NORMAL

  join_vec,['ptmt1_normal','pt_rat_557.7'],'vs5577'
  options,'vs5577','labels',['427.8 nm','557.7 nm']


  join_vec,['ptmt1_normal','pt_rat_670w'],'vs670w'
  options,'vs670w','labels',['427.8 nm','670w']

  join_vec,['ptmt1_normal','pt_rat_777.9'],'vs7779'
  options,'vs7779','labels',['427.8 nm','777.9 nm']

  join_vec,['ptmt1_normal','pt_rat_844.6'],'vs8446'
  options,'vs8446','labels',['427.8 nm','844.6 nm']
END

PRO read_h_info

  COMMON HERI_COUNT
  
  READ,types,PROMPT='Which types do you wanna check ?? : '
  IF types EQ 0 THEN file='/home/tkishiyama/spedas/heri_count2/count_after.txt'
  IF types EQ 1 THEN file='/home/tkishiyama/spedas/heri_count2/count_before.txt'
  
  file_lines=file_lines(file)
  openr,1,file
  print,'--- Open(r) = ',file

  count1=DBLARR(file_lines/7)
  count2=DBLARR(file_lines/7)
  count3=DBLARR(file_lines/7)
  count4=DBLARR(file_lines/7)
  count5=DBLARR(file_lines/7)
  count6=DBLARR(file_lines/7)
  count7=DBLARR(file_lines/7)

  time1=DBLARR(file_lines/7)
  time2=DBLARR(file_lines/7)
  time3=DBLARR(file_lines/7)
  time4=DBLARR(file_lines/7)
  time5=DBLARR(file_lines/7)
  time6=DBLARR(file_lines/7)
  time7=DBLARR(file_lines/7)

  abs_time1=DBLARR(file_lines/7)
  abs_time2=DBLARR(file_lines/7)
  abs_time3=DBLARR(file_lines/7)
  abs_time4=DBLARR(file_lines/7)
  abs_time5=DBLARR(file_lines/7)
  abs_time6=DBLARR(file_lines/7)
  abs_time7=DBLARR(file_lines/7)

  abs_count1=DBLARR(file_lines/7)
  abs_count2=DBLARR(file_lines/7)
  abs_count3=DBLARR(file_lines/7)
  abs_count4=DBLARR(file_lines/7)
  abs_count5=DBLARR(file_lines/7)
  abs_count6=DBLARR(file_lines/7)
  abs_count7=DBLARR(file_lines/7)
  
  tmp_read=DBLARR(4)
  time_heri=DBLARR(file_lines)
  count_heri=DBLARR(file_lines)
  abTi_heri=DBLARR(file_lines)
  abCt_heri=DBLARR(file_lines)
  
  i=0
  WHILE NOT EOF(1) DO BEGIN    
    readf,1,tmp_read
    time_heri(i)=tmp_read(0)
    count_heri(i)=tmp_read(1)
    abTi_heri(i)=tmp_read(2)
    abCt_heri(i)=tmp_read(3)

    print,tmp_read
    i++
  ENDWHILE


  FOR i=0, file_lines/7-1 DO BEGIN
    ;count data is devided
    count1(i)=count_heri(7*i)
    count2(i)=count_heri(7*i+1)
    count3(i)=count_heri(7*i+2)
    count4(i)=count_heri(7*i+3)
    count5(i)=count_heri(7*i+4)
    count6(i)=count_heri(7*i+5)
    count7(i)=count_heri(7*i+6)

    ;time data is devided
    time1(i)=time_heri(7*i)
    time2(i)=time_heri(7*i+1)
    time3(i)=time_heri(7*i+2)
    time4(i)=time_heri(7*i+3)
    time5(i)=time_heri(7*i+4)
    time6(i)=time_heri(7*i+5)
    time7(i)=time_heri(7*i+6)

    ;ab-time data is devided
    abs_time1(i)=abTi_heri(7*i)
    abs_time2(i)=abTi_heri(7*i+1)
    abs_time3(i)=abTi_heri(7*i+2)
    abs_time4(i)=abTi_heri(7*i+3)
    abs_time5(i)=abTi_heri(7*i+4)
    abs_time6(i)=abTi_heri(7*i+5)
    abs_time7(i)=abTi_heri(7*i+6)


    abs_count1(i)=abCt_heri(7*i)
    abs_count2(i)=abCt_heri(7*i+1)
    abs_count3(i)=abCt_heri(7*i+2)
    abs_count4(i)=abCt_heri(7*i+3)
    abs_count5(i)=abCt_heri(7*i+4)
    abs_count6(i)=abCt_heri(7*i+5)
    abs_count7(i)=abCt_heri(7*i+6)
  ENDFOR
  close,1
END

PRO cal_heri_per

  COMMON HERI_COUNT
  COMMON HERI_PLOT

  READ,types,PROMPT='Which types do you wanna calculate ?? : '
  percentage=DBLARR(N_ELEMENTS(time1))
  
  FOR i=0, N_ELEMENTS(time1)-1 DO BEGIN
    IF types EQ 0 THEN BEGIN
      IF count2(i) GT count3(i) THEN BEGIN
        denominator=count2(i)-count1(i)
        ENDIF ELSE BEGIN
          denominator=count3(i)-count1(i)
        ENDELSE
        numerator=count4(i)-count5(i)
        percentage(i)=numerator/denominator
    ENDIF
    IF types EQ 1 THEN BEGIN
      IF count5(i) GT count6(i) THEN BEGIN
        denominator=count5(i)-count4(i)
      ENDIF ELSE BEGIN
        denominator=count6(i)-count4(i)
      ENDELSE
      numerator=count2(i)-count3(i)
      percentage(i)=numerator/denominator
    ENDIF
  ENDFOR

  print,percentage

  per_ave=0d

  FOR i=0, N_ELEMENTS(percentage)-1 DO BEGIN
     per_ave+=percentage(i)
  ENDFOR

  per_info=percentage
  
  print,''
  per_ave=per_ave/N_ELEMENTS(percentage)
  print,'Number of herisugi is          ',N_ELEMENTS(percentage)
  print,'Average of percentage (OFF/ON) is ',per_ave
  print,''
        
END

PRO cal_heri_interval

  COMMON HERI_COUNT

  READ,types,PROMPT='Which types do you wanna check ?? : '
  
  interval=DBLARR(N_ELEMENTS(time1))
  ave_ON=DBLARR(N_ELEMENTS(time1))
  ave_OFF=DBLARR(N_ELEMENTS(time1))
  
  FOR i=0, N_ELEMENTS(time1)-1 DO BEGIN
    IF types EQ 0 THEN BEGIN
      interval(i)=abs_time5(i)-abs_time3(i)
      ave_ON(i)=abs_time3(i)
      ave_OFF(i)=abs_time6(i)-abs_time3(i)
    ENDIF
    IF types EQ 1 THEN BEGIN
      interval(i)=abs_time2(i)
    ENDIF
  ENDFOR
  print,interval

  per_interval=0d
  per_ON=0d
  per_OFF=0d
  FOR i=0, N_ELEMENTS(interval)-1 DO BEGIN
    per_interval+=interval(i)
    per_ON+=ave_ON(i)
    per_OFF+=ave_OFF(i)
  ENDFOR

  print,''
  per_interval=per_interval/N_ELEMENTS(interval)
  per_ON=per_ON/N_ELEMENTS(ave_ON)
  per_OFF=per_OFF/N_ELEMENTS(ave_OFF)
  print,'Number of herisugi is          ',N_ELEMENTS(interval)
  print,'Average of herisugi time is ',per_interval
  print,''
  print,'average ON is ',per_ON
  print,'average OFF is ',per_OFF
END

PRO cal_heri_bg

  COMMON HERI_COUNT
  COMMON HERI_PLOT

  READ,types,PROMPT='Which types do you wanna check ?? : '
  background_ave=DBLARR(N_ELEMENTS(time1))

  FOR i=0, N_ELEMENTS(time1)-1 DO BEGIN
    IF types EQ 0 THEN BEGIN
      background_ave(i)=(count1(i)+count4(i)+count6(i)+count7(i))/4
    ENDIF
    IF types EQ 1 THEN BEGIN
      background_ave(i)=(count1(i)+count2(i)+count4(i)+count7(i))/4
    ENDIF   
  ENDFOR

  print,background_ave
  bg_info=background_ave
  bg_ave=0d
  FOR i=0, N_ELEMENTS(count1)-1 DO BEGIN
    bg_ave+=background_ave(i)
  ENDFOR
  bg_ave=bg_ave/N_ELEMENTS(count1)

  print,''
  print,'Number of herisugi is  ',N_ELEMENTS(background_ave)
  print,'Average of background is ',bg_ave
  print,''
END

PRO barplot_background

  COMMON HERI_COUNT
  COMMON HERI_PLOT

  number_heri=DBLARR(10)
  FOR i=0, N_ELEMENTS(bg_info)-1 DO BEGIN
    bg_judge=bg_info(i)
    IF bg_judge LT 500 THEN number_heri(0)+=1
    IF bg_judge GE 500  && bg_judge LT 1000 THEN number_heri(1)+=1
    IF bg_judge GE 1000 && bg_judge LT 1500 THEN number_heri(2)+=1
    IF bg_judge GE 1500 && bg_judge LT 2000 THEN number_heri(3)+=1
    IF bg_judge GE 2000 && bg_judge LT 2500 THEN number_heri(4)+=1
    IF bg_judge GE 2500 && bg_judge LT 3000 THEN number_heri(5)+=1
    IF bg_judge GE 3000 && bg_judge LT 3500 THEN number_heri(6)+=1
    IF bg_judge GE 3500 && bg_judge LT 4000 THEN number_heri(7)+=1
    IF bg_judge GE 4000 && bg_judge LT 4500 THEN number_heri(8)+=1
    IF bg_judge GE 4500 && bg_judge LT 5000 THEN number_heri(9)+=1
  ENDFOR

  heri_background=DBLARR(10)
  FOR i=0, 9 DO BEGIN
    heri_background(i)=500*i
  ENDFOR

  READ,title1,PROMPT='after? or before? type which ! : '
  IF title1 EQ 0 THEN title2='Background level and Number of Large-decreasing after ON'
  IF title1 EQ 1 THEN title2='Background and Number of Large-decreasing before ON'
  b1=BARPLOT(heri_background,number_heri,INDEX=0,NBARS=1,FILL_COLOR='gold',YRANGE=[0.,100.],YMINOR=0.,$
         YTITLE='Number of Over-darkening',XTITLE='Count',TITLE=title2)
END

PRO barplot_time
  
  COMMON HERI_COUNT
  COMMON HERI_PLOT

  for_TiJudge=STRMID(time_string(time1),11,4)     
  number_heri=DBLARR(18)
  FOR i=0, N_ELEMENTS(time1)-1 DO BEGIN
    Ti_judge=for_TiJudge(i)
    FOR j=19, 23 DO BEGIN
      k=2*j-38
      l=2*j-37

      n=STRING(j,FORMAT='(I2)')
      m0=STRMID(n,0,2)+':0'
      m1=STRMID(n,0,2)+':1'
      m2=STRMID(n,0,2)+':2'
      IF Ti_judge EQ m0 THEN number_heri(k)+=1
      IF Ti_judge EQ m1 THEN number_heri(k)+=1
      IF Ti_judge EQ m2 THEN number_heri(k)+=1
      print,m0,m1,m2

      m0=STRMID(n,0,2)+':3'
      m1=STRMID(n,0,2)+':4'
      m2=STRMID(n,0,2)+':5'
      print,m0,m1,m2
      IF Ti_judge EQ m0 THEN number_heri(l)+=1
      IF Ti_judge EQ m1 THEN number_heri(l)+=1
      IF Ti_judge EQ m2 THEN number_heri(l)+=1

    ENDFOR
    FOR j=0,4 DO BEGIN
      k=10+2*j
      l=11+2*j

      n=STRING(j,FORMAT='(I02)')
      m0=STRMID(n,0,2)+':0'
      m1=STRMID(n,0,2)+':1'
      m2=STRMID(n,0,2)+':2'
      IF Ti_judge EQ m0 THEN number_heri(k)+=1
      IF Ti_judge EQ m1 THEN number_heri(k)+=1
      IF Ti_judge EQ m2 THEN number_heri(k)+=1

      m0=STRMID(n,0,2)+':3'
      m1=STRMID(n,0,2)+':4'
      m2=STRMID(n,0,2)+':5'
      IF Ti_judge EQ m0 THEN number_heri(l)+=1
      IF Ti_judge EQ m1 THEN number_heri(l)+=1
      IF Ti_judge EQ m2 THEN number_heri(l)+=1
    ENDFOR
  ENDFOR

  heri_jikantai=DBLARR(18)
  FOR i=0, 17 DO BEGIN
    heri_jikantai(i)=19+0.5*i
  ENDFOR
  heri_jikantai+=2.5
  READ,title1,PROMPT='after? or before? type which ! : '
  IF title1 EQ 0 THEN title2='Occurrence of PsA and Number of over-darkening after ON'
  IF title1 EQ 1 THEN title2='Occurrence of PsA and Number of over_darkening before ON'
  number_PsA=DBLARR(20)
  number_PsA=[2.5,5,7.5,10,10,5,7.5,5,10,13,20,17.5,17.5,27,25,30,32.5,32.5]
     
  b2=BARPLOT(heri_jikantai,number_PsA,INDEX=1,NBARS=2,FILL_COLOR='green',YRANGE=[0.,15])
  b1=BARPLOT(heri_jikantai,number_heri,INDEX=0,NBARS=2,FILL_COLOR='orange',YRANGE=[0.,40.],YMINOR=0.,$
            YTITLE='Number of over-darkening',XTITLE='MLT',$
            TITLE=title2,XTHICK=3,YTHICK=3,FONT_SIZE=23,/overplot)

  yaxis = AXIS('Y', LOCATION='right',TARGET=number_PsA, $
              TITLE='Total time of PsA occurence!C[Hour]', $
              COORD_TRANSFORM=[0,0.4], $
              AXIS_RANGE=[0,15],TICKFONT_SIZE=20)

  text_pine = TEXT(22, 35, 'Over-darkening', /CURRENT, $
                  COLOR='orange', FONT_SIZE=23, /DATA)
  text_yel = TEXT(22, 32, 'Pulsating aurora', /CURRENT, $
                 COLOR = 'green',FONT_SIZE=23, /DATA) 
END

PRO scatplot_background

  COMMON HERI_COUNT
  COMMON HERI_PLOT
       
  cal_heri_per
  cal_heri_bg

  distSize=SCATTERPLOT(bg_info,per_info,SYMBOL='x',/SYM_FILLED, RGB_TABLE=0,$
                      YTITLE=' Percentage of OFF/ON ',XTITLE=' Background ',$
                      TITLE='Background and Percentage of OFF/ON')
END

PRO plot_heri_shape

  COMMON HERI_COUNT
  COMMON HERI_PLOT

  new_count_heri=DBLARR(7)
  new_time_heri=DBLARR(7)

  nor_count=DBLARR(N_ELEMENTS(count_heri))
  nor_time=DBLARR(N_ELEMENTS(abTi_heri))

  c_ave1=DBLARR(N_ELEMENTS(time1))
  c_ave2=DBLARR(N_ELEMENTS(time1))
  c_ave3=DBLARR(N_ELEMENTS(time1))
  c_ave4=DBLARR(N_ELEMENTS(time1))
  c_ave5=DBLARR(N_ELEMENTS(time1))
  c_ave6=DBLARR(N_ELEMENTS(time1))
  c_ave7=DBLARR(N_ELEMENTS(time1))

  t_ave1=DBLARR(N_ELEMENTS(time1))
  t_ave2=DBLARR(N_ELEMENTS(time1))
  t_ave3=DBLARR(N_ELEMENTS(time1))
  t_ave4=DBLARR(N_ELEMENTS(time1))
  t_ave5=DBLARR(N_ELEMENTS(time1))
  t_ave6=DBLARR(N_ELEMENTS(time1))
  t_ave7=DBLARR(N_ELEMENTS(time1))

  c_ave=DBLARR(7)
  t_ave=DBLARR(7)
  abTi_norm=DBLARR(7)

  PLOT,[0],[0],max_value=1.0,min_value=0.0,thick=3,charsize=1.7,$
      XRANGE=[-10.0,25.0],YRANGE=[0,1.0],$
      TITLE='Average shape of over-darkening after ON',$
      XTITLE='Sec',YTITLE='Normalized count',$
      XTHICK=3,YTHICK=3
   
  FOR i=0, N_ELEMENTS(time1)-1 DO BEGIN
    new_time_heri(0)=0.0
    ;abTi_norm=normalizer(abTi_heri(7*i :7*i+6))
    new_count_heri=normalizer(count_heri(7*i :7*i+6))
    new_time_heri(1:6)=abTi_heri(7*i :7*i+5)

    nor_count(7*i :7*i+6)=new_count_heri
    nor_time(7*i :7*i+6)=new_time_heri
    nor_time(7*i :7*i+6)-=nor_time(7*i+3)
    
    c_ave1(i)=new_count_heri(0)
    c_ave2(i)=new_count_heri(1)
    c_ave3(i)=new_count_heri(2)
    c_ave4(i)=new_count_heri(3)
    c_ave5(i)=new_count_heri(4)
    c_ave6(i)=new_count_heri(5)
    c_ave7(i)=new_count_heri(6)

    t_ave1(i)=nor_time(7*i)
    t_ave2(i)=nor_time(7*i+1)
    t_ave3(i)=nor_time(7*i+2)
    t_ave4(i)=nor_time(7*i+3)
    t_ave5(i)=nor_time(7*i+4)
    t_ave6(i)=nor_time(7*i+5)
    t_ave7(i)=nor_time(7*i+6)

    oplot,nor_time(7*i :7*i+6),nor_count(7*i :7*i+6),thick=0.9  
  ENDFOR
  FOR j=0, N_ELEMENTS(time1)-1 DO BEGIN
    c_ave(0)+=c_ave1(j)
    c_ave(1)+=c_ave2(j)
    c_ave(2)+=c_ave3(j)
    c_ave(3)+=c_ave4(j)
    c_ave(4)+=c_ave5(j)
    c_ave(5)+=c_ave6(j)
    c_ave(6)+=c_ave7(j)

    t_ave(0)+=t_ave1(j)
    t_ave(1)+=t_ave2(j)
    t_ave(2)+=t_ave3(j)
    t_ave(3)+=t_ave4(j)
    t_ave(4)+=t_ave5(j)
    t_ave(5)+=t_ave6(j)
    t_ave(6)+=t_ave7(j)
  ENDFOR
  FOR k=0, 6 DO BEGIN
    c_ave(k)=c_ave(k)/N_ELEMENTS(c_ave1)
    t_ave(k)=t_ave(k)/N_ELEMENTS(t_ave1)
  ENDFOR
  oplot,t_ave,c_ave,COLOR=230,THICK=14.0
  print,t_ave
END

PRO plot_mlt_per

	COMMON HERI_COUNT
	COMMON HERI_PLOT

	x_ut=DBLARR(10)
  FOR i=0, 9 DO BEGIN
    x_ut(i)=19+i
	ENDFOR
	x_ut+=2.5      	
	y_offon=DBLARR(10)
	count=DBLARR(10)	
	timeJudge=STRMID(time_string(time1),11,2)
	FOR i=0, N_ELEMENTS(time1)-1 DO BEGIN
    FOR j=0, 4 DO BEGIN
      nextDay=j
      J_nextDay=STRING(nextDay,FORMAT='(I02)')
      today=19+j
      J_today=STRING(today,FORMAT='(I02)')
      IF timeJudge(i) EQ J_nextDay THEN BEGIN
        y_offon(5+j)+=per_info(i)
        count(5+j)++   
      ENDIF
      IF timeJudge(i) EQ J_today THEN BEGIN
        y_offon(j)+= per_info(i)
        count(j)++
      ENDIF    
    ENDFOR         
  ENDFOR

  FOR i=0, 8 DO BEGIN
    y_offon(i) = y_offon(i)/count(i)
  ENDFOR
  print,y_offon

	PLOT,[0],[0],max_value=1.0,min_value=0.0,thick=3,charsize=1.7,$
            XRANGE=[21,30],YRANGE=[20,80],$
            TITLE='MLT and percentage of over-darkening after ON',$
            XTITLE='MLT',YTITLE='Percentage of over-darkening!C [%]',$
            XTHICK=3,YTHICK=3
	oplot,x_ut(0:N_ELEMENTS(x_ut)-2),y_offon(0:N_ELEMENTS(y_offon)-2)*100,thick=7.0
END


      
