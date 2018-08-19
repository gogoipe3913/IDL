;--------------------------------------------------------------------------------------------------------------------------
; COMMON BLOCKS
;

	COMMON prf_psa, stm_psa,etm_psa,max_psa,min_psa,fmt_psa,dir_psa,tr_sat,tr_img

	COMMON fov_psa,	azm_psa,ele_psa,gla_psa,glo_psa

	COMMON raw_psa,	cid_raw,loc_raw,flt_raw,fps_raw,day_raw,num_raw,img_raw,typ_raw,now_raw, $
                        yr_raw,mo_raw,dy_raw,hh_raw,mm_raw,ss_raw,ms_raw,ts_raw,ti_raw,YYYY_orb, $
                        MM_orb,DD_orb,hour_orb,min_orb,dura_orb,fft_peak,smo_data,ms_changed,lun,fft_hist_img,fft_hist_fbk

	COMMON tif_psa,	cid_tif,loc_tif,flt_tif,fps_tif,day_tif,num_tif,img_tif,typ_tif,now_tif, $
                        yr_tif,mo_tif,dy_tif,hh_tif,mm_tif,ss_tif,ms_tif,ts_tif,ti_tif

        COMMON cor_psa, max_cor,min_cor,max_pixel,max_pixel_cor,cor,img_changed,img_changed_tim, $
                        img_changed_tim_ms,fbk_data_changed,fbk_tim_changed,cal_span

        COMMON struct_psa, geogra_map,cor_contour,fft_contour

        COMMON apa_psa, day_apa,num_apa,img_apa,now_apa,yr_apa,mo_apa,dy_apa,hh_apa,mm_apa,ss_apa,hms_apa,tsc_apa,siz_apa_x,siz_apa_y,time_str_apa,img_apa_cut,time_str_apa_cut,span
                        
        COMMON apa_fov, gla_apa,glo_apa,azm_apa,ele_apa

;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	START_PSA
;
; PURPOSE:
;
;	Start up the Pulsating Aurora Data plotting environment on IDL
;

	PRO start_psa

          COMMON prf_psa
          COMMON cor_psa
          COMMON raw_psa
   
; arai san
          plot_aa_routine,'kemono_arai'

; Set environment things
	DEVICE,DECOMPOSED=0
        !P.NOERASE=1
        !PROMPT='PsA > '
        !MSG_PREFIX='IDL: 

; Set min and max values
	min_psa=2000
	max_psa=2500

; Set time resolution of img and satellite data
        tr_sat=8
        tr_img=100

; Set window shift interval of correlate coefficient
        cal_span=10

; Set orbit span
        YYYY_orb='2016'
        MM_orb='10'
        DD_orb='04'
        hour_orb='22'
        min_orb='00'
        dura_orb=3
; Format
	fmt_psa=0

; Data path
	dir_psa='/raid/data/Kiban-S/'

; Load the default colour table
	LOADCT,0

	END

;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	NO_TICKS
;
; FUNCTION:
;
;	Never mind the bollocks.
;

        FUNCTION no_ticks,      axis,index,time

        RETURN,''

        END

;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	AXIS_TIME
;
; FUNCTION:
;
;	Gives stupid x-axis	
;

        FUNCTION axis_time,axis,index,time

        hour=FIX(time) MOD 24
        mins=FIX(((time MOD 24)-hour)*60+0.01)
        secs=FIX((((time MOD 24)-hour)*3600+0.01) MOD 60)

        label=STRING(hour,mins,FORMAT='(I2.2,''!U'',I2.2,''!N'')')

        RETURN,label

        END

;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	FIND_NEXT_PREV_DAY
;
; FUNCTION:
;
;	Find next/prev date for given date
;

        FUNCTION find_next_prev_day,old_day,next=next,prev=prev

        IF KEYWORD_SET(next) THEN offset=+1
        IF KEYWORD_SET(prev) THEN offset=-1

        yy=FIX(STRMID(old_day,0,4))
        mm=FIX(STRMID(old_day,4,2))
        dd=FIX(STRMID(old_day,6,2))
        jday_tommorrow=JULDAY(mm,dd,yy)+offset
        CALDAT,jday_tommorrow,mf,df,yf
        new_day=STRING(yf,FORMAT='(I4.4)')+STRING(mf,FORMAT='(I2.2)')+STRING(df,FORMAT='(I2.2)')

        RETURN,new_day

        END

;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       RANGE_TO_COORDS
;
; NOTES:
;
;       Given a start "lat" and "lon" find the position of the point "range" in km away
;       at an azimuth of "azimuth" from north.
;

        FUNCTION range_to_coords,lat,lon,azimuth,range

        IF azimuth GT 180 THEN az=azimuth-360 ELSE az=azimuth
        Re=6370.0
        coLat=90-lat

        coLat_point=ACOS(COS(range/Re)*COS(coLat*!DTOR)+SIN(range/Re)* $
                SIN(coLat*!DTOR)*COS(az*!DTOR))*!RADEG
        lon_point=ACOS(MIN([(COS(range/Re)-COS(coLat_point*!DTOR)* $
                COS(coLat*!DTOR))/(SIN(coLat_point*!DTOR)*SIN(coLat*!DTOR)),1.0]))*!RADEG
        lat_point=90-coLat_point
        IF az GT 0 THEN lon_point=lon+lon_point ELSE lon_point=lon-lon_point

        RETURN,[lat_point,lon_point]

        END

;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       ASI_FOV_TO_GEOG
;
; NOTE:
;       Given an ASI position (geog coords) and azimuth and the zenith of observations
;       and the assumed emission altitude return the geographic lat and lon.
;

        FUNCTION asi_fov_to_geog,lat,lon,azimuth,zenith,altitude

; Earth radius
        Re=6370.0

; Assume that msp look direction maps to within +/-25 degrees of the station
; and determine heights that projection maps to within this range of lats
        delta_lat=25*FINDGEN(1000)/1000
        h=Re*(TAN((90-ABS(zenith)+delta_lat)*!pi/180)*SIN(delta_lat*!pi/180)+COS(delta_lat*!pi/180)-1)

; Find h that corresponds best to desired altitude
        min_dev=MIN(ABS(altitude-h),min_pos)

; Convert position to a distance
        range=Re*delta_lat(min_pos)*!pi/180

; Determine latitude and longitude of distance from msp location at given azimuth
        IF zenith LT 0 THEN az_point=azimuth+180 ELSE az_point=azimuth

        RETURN,range_to_coords(lat,lon,az_point,range)

        END

;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	PSYM8
;
; PURPOSE:
;
;	Setting original symbol

        PRO psym8,fill=fill,dot=dot

        cir=findgen(31)*(!pi*2/30)
        IF KEYWORD_SET(fill) THEN $
                USERSYM,cos(cir),sin(cir),/fill ELSE $
                USERSYM,cos(cir),sin(cir)
        END

;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	IDENTIFY_PSA_CAM
;

	PRO identify_psa_cam,camid,date,loc_out,flt_out,fps_out

	date_long=LONG(date)

	IF camid EQ 1 THEN BEGIN
		loc_out='Tromso'
		flt_out='777.4 nm'
		IF date_long GE 20170126 THEN flt_out='RG665'
		fps_out=10
		IF date_long GE 20170126 THEN fps_out=100
	ENDIF

	IF camid EQ 2 THEN BEGIN
		loc_out='Sodankyla'
		flt_out='RG665'
		fps_out=100
	ENDIF

	IF camid EQ 3 THEN BEGIN
		loc_out='Tromso'
		flt_out='844.6 nm'
		fps_out=10
	ENDIF

	IF camid EQ 4 THEN BEGIN
		loc_out='Tromso'
		flt_out='N21PG Wide'
		fps_out=10
	ENDIF

	IF camid EQ 5 THEN BEGIN
		loc_out='Tromso'
		flt_out='427.8 nm'
		fps_out=10
	ENDIF

	IF camid EQ 6 THEN BEGIN
		loc_out='Kevo'
		flt_out='RG665'
		fps_out=100
	ENDIF

	IF camid EQ 7 THEN BEGIN
		loc_out='Gakona'
		flt_out='RG665'
		fps_out=100
	ENDIF

	END

;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	FILE_PSA_RAW
;
; EXAMPLE:
;
;	To read the raw count data:
;
;	PsA > file_psa_raw,1,'20161030',2355,10
;
;	To read the absolute optical intensity data in Rayleigh:
;
;	PsA > file_psa_raw,1,'20161030',2355,10,/abs
;

	PRO file_psa_raw,camid,date,hhmm,dura,pxl=pxl,abs=abs,vbs=vbs

	COMMON prf_psa
	COMMON raw_psa

; set time for spedas

        YYYY=STRMID(date,0,4)
        MM=STRMID(date,4,2)
        DD=STRMID(date,6,2)
        hour=STRING(FIX(hhmm)/100)
	min=STRING(FIX(hhmm-100*hour))

        ;timespan,YYYY+'-'+MM+'-'+DD+'/'+hour+':'+min+':00',dura,/minute
; Check the camera identity
	identify_psa_cam,camid,date,loc_out,flt_out,fps_out
	loc_raw=loc_out
	flt_raw=flt_out
	fps_raw=fps_out
	PRINT,''
	PRINT,'Read the following data:'
	PRINT,' Date:      '+date
	PRINT,' Location:  '+loc_raw
	PRINT,' Filter:    '+flt_raw
	PRINT,' Sampling:  '+STRING(fps_raw,FORMAT='(I3.3)')+' fps'
	IF KEYWORD_SET(abs) THEN PRINT,' Data type: Absolute optical intensity from .raw file' $
		ELSE PRINT,' Data type: Raw count from .raw file'

; Invalid request
	IF fps_raw EQ 100 AND KEYWORD_SET(abs) THEN BEGIN
		PRINT,''
		PRINT,'ERROR: Absolute intensity data are not available for 100 Hz cameras'
		PRINT,''
		RETURN
	ENDIF

; Identify camera id
	cid_raw=camid

; Get the calibration table if needed
	IF KEYWORD_SET(abs) THEN BEGIN
		IF camid EQ 1 THEN cal_fname=dir_psa+'ctbl/C'+STRING(camid,FORMAT='(I1.1)')+'_7774_cal.sav'
		IF camid EQ 3 THEN cal_fname=dir_psa+'ctbl/C'+STRING(camid,FORMAT='(I1.1)')+'_8446_cal.sav'
		IF camid EQ 4 THEN cal_fname=dir_psa+'ctbl/C'+STRING(camid,FORMAT='(I1.1)')+'_670W_cal.sav'
		IF camid EQ 5 THEN cal_fname=dir_psa+'ctbl/C'+STRING(camid,FORMAT='(I1.1)')+'_4278_cal.sav'
		PRINT,''
		PRINT,'Read Calibration Table:'
		PRINT,' '+cal_fname
		RESTORE,cal_fname
		IF LONG(date) GE 20170126 THEN BEGIN
			PRINT,''
			PRINT,' NOTE:'
			PRINT,' Due to the change of the formation of the cameras in Tromso, the file you read'
			PRINT,' is tantative. The table will be modified after the calibration in July, 2017'
		ENDIF
	ENDIF

; Get the FOV information
	file_psa_fov,camid,date

; Keyword setting
	IF NOT KEYWORD_SET(pxl) THEN pxl=256

; Make arrays
	nmax=60L*fps_raw*dura
	img_raw=FLTARR(nmax,pxl,pxl)
	day_raw=STRARR(nmax)
	yr_raw=INTARR(nmax)
	mo_raw=INTARR(nmax)
	dy_raw=INTARR(nmax)
	hh_raw=INTARR(nmax)
	mm_raw=INTARR(nmax)
	ss_raw=INTARR(nmax)
	ms_raw=INTARR(nmax)
	ts_raw=FLTARR(nmax)
	ti_raw=STRARR(nmax)

; Start time
	hh=FIX(hhmm)/100
	mm=FIX(hhmm-100*hh)

; Flag for correct or wrong data format
	flg=0

; Loop for files
	PRINT,''
	PRINT,'Read raw files:'
	n=0L
	FOR i=0,dura-1 DO BEGIN

; Identify time (hhmm) of the raw file
		mm_next=mm+i
		xm=mm_next/60
		ym=mm_next MOD 60
		hh_next=hh+xm
		mm_next=ym
		date_next=date
		IF hh_next GE  24 THEN BEGIN
			hh_next=hh_next MOD 24
			date_next=find_next_prev_day(date,/next)
		ENDIF
		hhmm_next=hh_next*100+mm_next

; Generate path for the data
		yr_read=STRMID(date_next,0,4)
		mo_read=STRMID(date_next,4,2)
		dy_read=STRMID(date_next,6,2)
		path=dir_psa+'cam'+STRING(cid_raw,FORMAT='(I1.1)')+'/'+yr_read+'/'+mo_read+'/'+dy_read
		fname=path+'/C'+STRING(camid,FORMAT='(I1.1)')+'_'+date_next+'_'+STRING(hhmm_next,FORMAT='(I4.4)')+'.raw'
		;PRINT,''
		;PRINT,'Now reading '+fname
		PRINT,' '+fname

; Open the raw file
		OPENR,id,fname,/GET_LUN

; Loop for images
		WHILE NOT EOF(id) DO BEGIN

			IF n GE nmax THEN GOTO,skip0	

; Call ebiread_ym and read one image
			ebireaded_ym,id,itag,x,y,imgname,dat
	
			IF itag eq 2000 then BEGIN

; Preview
				IF n EQ 0 AND STRMID(imgname,5,1) EQ 'd' THEN flg=1

; Extract time info
				IF flg EQ 1 THEN BEGIN
					yr_raw(n)=FIX(STRMID(imgname,6,4))
					mo_raw(n)=FIX(STRMID(imgname,10,2))
					dy_raw(n)=FIX(STRMID(imgname,13,2))
					hh_raw(n)=FIX(STRMID(imgname,15,2))
					mm_raw(n)=FIX(STRMID(imgname,17,2))
					ss_raw(n)=FIX(STRMID(imgname,20,3))
					ms_raw(n)=FIX(STRMID(imgname,24,3))
				ENDIF ELSE BEGIN
					yr_raw(n)=FIX(STRMID(imgname,3,4))
					mo_raw(n)=FIX(STRMID(imgname,7,2))
					dy_raw(n)=FIX(STRMID(imgname,9,2))
					hh_raw(n)=FIX(STRMID(imgname,12,2))
					mm_raw(n)=FIX(STRMID(imgname,14,2))
					ss_raw(n)=FIX(STRMID(imgname,16,3))
					ms_raw(n)=FIX(STRMID(imgname,19,3))
				ENDELSE
				day_raw(n)=STRING(yr_raw(n),FORMAT='(I4.4)')+ $
                                           STRING(mo_raw(n),FORMAT='(I2.2)')+STRING(dy_raw(n),FORMAT='(I2.2)')
				ts_raw(n)=3600.0*hh_raw(n)+60.0*mm_raw(n)+ss_raw(n)+ms_raw(n)/1000.0
				IF date NE date_next THEN ts_raw(n)=ts_raw(n)+86400.0
				hh_str=STRING(hh_raw(n),FORMAT='(I2.2)')
				mm_str=STRING(mm_raw(n),FORMAT='(I2.2)')
				ss_str=STRING(ss_raw(n),FORMAT='(I2.2)')
				ms_str=STRING(ms_raw(n),FORMAT='(I3.3)')
        			ti_raw(n)=hh_str+mm_str+' '+ss_str+'s ('+ms_str+' ms)'

				IF KEYWORD_SET(vbs) THEN $
    					PRINT,STRING(n,FORMAT='(I6.6)')+' th image: ', $
						yr_raw(n),mo_raw(n),dy_raw(n),hh_raw(n),mm_raw(n),ss_raw(n),ms_raw(n)

; Put the data into array
				IF KEYWORD_SET(abs) THEN BEGIN
					img_raw(n,*,*)=(ROTATE(dat,7)-dk_fit)/at_fit
				ENDIF ELSE BEGIN
					img_raw(n,*,*)=ROTATE(dat,7)
				ENDELSE

				;IF n GT nmax-1 THEN GOTO,skip0
    				n=n+1
    			ENDIF
		ENDWHILE
		skip0:
		num_raw=n

; Close the file
		FREE_LUN,id

	ENDFOR

; Total number of files
	PRINT,''
	PRINT,' >>> The total number of stored images:'+STRCOMPRESS(STRING(num_raw,FORMAT='(I6)'))
	PRINT,''

; Setting type and corresponding limit
	IF KEYWORD_SET(abs) THEN BEGIN
		typ_raw=1
		IF min_psa GT 1500 THEN set_scale_psa,0,1000,/quiet
	ENDIF ELSE BEGIN
		typ_raw=0
		IF min_psa LT 1500 THEN set_scale_psa,2000,2500,/quiet
	ENDELSE

; Setting current image
	go_psa_raw,0
      	END
;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       FILE_PSA_APA

        PRO file_psa_apa,path=path

          COMMON apa_psa
          COMMON apa_fov
	
	
	IF NOT KEYWORD_SET(path) THEN BEGIN
	     path=''
	     READ,path,PROMPT='Enter file path: '
          ENDIF

        span=300
        shift=10
        time_stamp,/off

        found_file=FILE_SEARCH(path+'*.txt',COUNT=num_found_file)
        num_apa=num_found_file

; Values of fish eye lens

; Initialize arrays
        max_apa=4000
        min_apa=500
        max_cor=1.0
        min_cor=0.0
	siz_apa_x=383
	siz_apa_y=287
        img_apa=FLTARR(num_apa,siz_apa_x,siz_apa_y)
        keo_apa=FLTARR(num_apa,siz_apa_y)
        keo_col=FLTARR(num_apa,siz_apa_y)
        day_apa=STRARR(num_apa)
        yr_apa=STRARR(num_apa)
        mo_apa=STRARR(num_apa)
        dy_apa=STRARR(num_apa)
        hms_apa=STRARR(num_apa)
	hh_apa=STRARR(num_apa)
	mm_apa=STRARR(num_apa)
	ss_apa=STRARR(num_apa)
        tsc_apa=LONARR(num_apa)
        time_str_apa=STRARR(num_apa)

; Read the data
	FOR i=0,num_apa-1 DO BEGIN
           day_apa(i)=STRMID(found_file(i),STRLEN(path),8)
           yr_apa(i)=STRMID(day_apa(i),0,4)
           mo_apa(i)=STRMID(day_apa(i),4,2)
           dy_apa(i)=STRMID(day_apa(i),6,2)
           hh_apa(i)=STRMID(found_file(i),STRLEN(path)+8,2)
           mm_apa(i)=STRMID(found_file(i),STRLEN(path)+10,2)
           ss_apa(i)=STRMID(found_file(i),STRLEN(path)+12,2)
           hms_apa(i)=hh_apa(i)+mm_apa(i)+ss_apa(i)
           tsc_apa(i)=3600L*FIX(hh_apa(i))+60L*FIX(mm_apa(i))+FIX(ss_apa(i))
           time_str_apa(i)=yr_apa(i)+'-'+mo_apa(i)+'-'+dy_apa(i)+'/'+hh_apa(i)+':'+mm_apa(i)+':'+ss_apa(i)
           
           tmp=READ_ASCII(found_file(i))
           img_tmp=tmp.FIELD001
           ;img_tmp=MEAN(tmp_read,dimension=1)
           img_apa(i,*,*)=img_tmp(*,*)
           keo_apa(i,*)=REFORM(img_apa(i,siz_apa_x/2,*))
           
           PRINT,'scan'+STRING(i,FORMAT='(I6)')+' : '+day_apa(i)+'-'+hms_apa(i)
           
        ENDFOR
        go_apa,0
        file_apa_fov

   END

;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	FILE_APA_FOV
;

	PRO file_apa_fov,alt=alt

	COMMON apa_psa
	COMMON apa_fov

	IF NOT KEYWORD_SET(alt) THEN alt=100

; Message
	PRINT,''
	PRINT,'Read the FOV files for '+STRING(alt,FORMAT='(I3)')+' km altitude:'
        
; Selection of FOV files
        fname_base='star_cal'
	
; Make arrays
	azm_apa=FLTARR(siz_apa_x,siz_apa_y)
        ele_apa=FLTARR(siz_apa_x,siz_apa_y)
	gla_apa=FLTARR(siz_apa_x,siz_apa_y)
	glo_apa=FLTARR(siz_apa_x,siz_apa_y)

; Read azimuth and elevation angles
	found_file=FILE_SEARCH('/home/suguru/spedas/image_data/russia/apa_fov/'+fname_base+'_ang'+'*.txt',COUNT=num_found_file)
	PRINT,' ',found_file(0)
	PRINT,' ',found_file(1)
	OPENR,azm_lun,found_file(0),/GET_LUN
        OPENR,ele_lun,found_file(1),/GET_LUN
        
        READF,ele_lun,ele_apa
        READF,azm_lun,azm_apa
        
        FREE_LUN,azm_lun
        FREE_LUN,ele_lun

; Read glat and glon info
	found_file=FILE_SEARCH('/home/suguru/spedas/image_data/russia/apa_fov/'+fname_base+'_'+STRING(alt,FORMAT='(I3.3)')+'*.txt', $
		COUNT=num_found_file)
	PRINT,' ',found_file(0)
	PRINT,' ',found_file(1)
	OPENR,gla_lun,found_file(0),/GET_LUN
        OPENR,glo_lun,found_file(1),/GET_LUN
        
        READF,gla_lun,gla_apa
        READF,glo_lun,glo_apa

        FREE_LUN,gla_lun
        FREE_LUN,glo_lun
        
        gla_apa=REVERSE(gla_apa,1)
        glo_apa=REVERSE(glo_apa,1)
        gla_apa=REVERSE(gla_apa,2)
        glo_apa=REVERSE(glo_apa,2)
	END
;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       GO_APA

        PRO go_apa,image_to_go

	COMMON apa_psa
	
	now_apa=0L
	now_apa=image_to_go

	current_hms=STRING(hms_apa(now_apa),FORMAT='(I6.6)')
	PRINT,'Time of current image: '+current_hms+' UT'

	END
           
;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	SAVE_PSA_RAW_ABS
;
; PURPOSE:
;
;	Save the absolute optical intensity data as individual 16-bit TIFF files
;
; EXAMPLE:
;
;	PsA > save_psa_raw_abs,1,'20161030',0125,10
;
;	TIFF images are automatically saved in a directory "C1_20161030'
;

	PRO save_psa_raw_abs,camid,date,hhmm,dura

	COMMON prf_psa
	COMMON raw_psa

	file_psa_raw,camid,date,hhmm,dura,/abs

	out_dir='C'+STRING(cid_raw,FORMAT='(I1.1)')+'_'+day_raw(0)
	SPAWN,'mkdir -p '+out_dir

	FOR i=0,num_raw-1 DO BEGIN

		sav_img=LONG(REFORM(img_raw(i,*,*)))
		zero_pix=WHERE(sav_img LT 0,num_zero_pix)
		sav_img(zero_pix)=0
		fname='C'+STRING(cid_raw,FORMAT='(I1.1)')+'_'+day_raw(i)+'_'+ $
			STRING(hh_raw(i),FORMAT='(I2.2)')+STRING(mm_raw(i),FORMAT='(I2.2)')+ $
			STRING(ss_raw(i),FORMAT='(I2.2)')+'_'+STRING(ms_raw(i),FORMAT='(I3.3)')+'.tif'
		WRITE_TIFF,out_dir+'/'+fname,sav_img,/SHORT
	ENDFOR

	END

;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	FILE_PSA_TIF
;
; EXAMPLE:
;
;	To read the QL tif data:
;
;	IDL > file_psa_tif,1,'20161030'
;
;	To read the absolute optical intensity data in Rayleigh:
;
;	IDL > file_psa_tif,1,'20161030',01,2,/abs
;

	PRO file_psa_tif,camid,date,hh,dura,pxl=pxl,abs=abs,vbs=vbs

	COMMON prf_psa
	COMMON tif_psa

; Check the camera identity
        identify_psa_cam,camid,date,loc_out,flt_out,fps_out
        loc_tif=loc_out
        flt_tif=flt_out
        fps_tif=fps_out
        PRINT,''
        PRINT,'Read the following data:'
        PRINT,' Date:      '+date
        PRINT,' Location:  '+loc_tif
        PRINT,' Filter:    '+flt_tif
        PRINT,' Sampling:  Every 10 sec (original sampling: '+STRING(fps_tif,FORMAT='(I3)')+' fps)'
        IF KEYWORD_SET(abs) THEN PRINT,' Data type: Absolute optical intensity from .tif file' $
                ELSE PRINT,' Data type: Raw count from .tif file'

; Invalid request
        IF fps_tif EQ 100 AND KEYWORD_SET(abs) THEN BEGIN
                PRINT,''
                PRINT,'ERROR: Absolute intensity data are not available for 100 Hz cameras'
                PRINT,''
                RETURN
        ENDIF

; Identify camera id
	cid_tif=camid

; Get the calibration table if needed
        IF KEYWORD_SET(abs) THEN BEGIN
                IF camid EQ 1 THEN cal_fname=dir_psa+'ctbl/C'+STRING(camid,FORMAT='(I1.1)')+'_7774_cal.sav'
                IF camid EQ 3 THEN cal_fname=dir_psa+'ctbl/C'+STRING(camid,FORMAT='(I1.1)')+'_8446_cal.sav'
                IF camid EQ 4 THEN cal_fname=dir_psa+'ctbl/C'+STRING(camid,FORMAT='(I1.1)')+'_670W_cal.sav'
                IF camid EQ 5 THEN cal_fname=dir_psa+'ctbl/C'+STRING(camid,FORMAT='(I1.1)')+'_4278_cal.sav'
                PRINT,''
                PRINT,'Read Calibration Table:'
                PRINT,' '+cal_fname
                RESTORE,cal_fname
                IF LONG(date) GE 20170126 THEN BEGIN
                        PRINT,''
                        PRINT,' NOTE:'
                        PRINT,' Due to the change of the formation of the cameras in Tromso, the file you read'
                        PRINT,' is tantative. The table will be modified after the calibration in July, 2017'
                ENDIF
        ENDIF

; File field-of-view
	file_psa_fov,camid,date

; Keyword setting
	IF NOT KEYWORD_SET(pxl) THEN pxl=256

; Make arrays
	nmax=360*dura
	img_tif=FLTARR(nmax,pxl,pxl)
	day_tif=STRARR(nmax)
	yr_tif=INTARR(nmax)
	mo_tif=INTARR(nmax)
	dy_tif=INTARR(nmax)
	hh_tif=INTARR(nmax)
	mm_tif=INTARR(nmax)
	ss_tif=INTARR(nmax)
	ms_tif=INTARR(nmax)
	ts_tif=FLTARR(nmax)
	ti_tif=STRARR(nmax)

; List up files
        PRINT,''
        PRINT,'Read tif files:'
	n=0L
	FOR k=0,dura-1 DO BEGIN

		hh_next=hh+k
		date_next=date
		IF hh_next/24 EQ 1 THEN BEGIN
			hh_next=hh_next MOD 24
			date_next=find_next_prev_day(date,/next)
		ENDIF

; Break date
		yr_read=STRMID(date_next,0,4)
		mo_read=STRMID(date_next,4,2)
		dy_read=STRMID(date_next,6,2)

; Path for the data
		path=dir_psa+'cam'+STRING(cid_tif,FORMAT='(I1.1)')+'/'+yr_read+'/'+mo_read+'/'+dy_read+'_QL'
                fname=path+'/C'+STRING(camid,FORMAT='(I1.1)')+'_'+date_next+'_'+STRING(hh_next,FORMAT='(I2.2)')+'*.tif'
		found_files=FILE_SEARCH(fname,COUNT=num_found_files)

; Generate filenames
		FOR i=0,num_found_files-1 DO BEGIN
		
; Extract time info
			IF i EQ 0 THEN BEGIN
				PRINT,' ',found_files(i)
				PRINT,'           .... '
			ENDIF
			IF i EQ num_found_files-1 THEN PRINT,' ',found_files(i)
			yr_tif(n)=STRMID(found_files(i),STRLEN(path)+4,4)
			mo_tif(n)=STRMID(found_files(i),STRLEN(path)+8,2)
			dy_tif(n)=STRMID(found_files(i),STRLEN(path)+10,2)
			hh_tif(n)=STRMID(found_files(i),STRLEN(path)+13,2)
			mm_tif(n)=STRMID(found_files(i),STRLEN(path)+15,2)
			ss_tif(n)=STRMID(found_files(i),STRLEN(path)+17,2)
			ms_tif(n)=STRMID(found_files(i),STRLEN(path)+20,3)
			ts_tif(n)=3600.0*hh_tif(n)+60.0*mm_tif(n)+ss_tif(n)+ms_tif(n)/1000.0
			IF date NE date_next THEN ts_tif(n)=ts_tif(n)+86400.0
			hh_str=STRING(hh_tif(n),FORMAT='(I2.2)')
			mm_str=STRING(mm_tif(n),FORMAT='(I2.2)')
			ss_str=STRING(ss_tif(n),FORMAT='(I2.2)')
			ms_str=STRING(ms_tif(n),FORMAT='(I3.3)')
       			ti_tif(n)=hh_str+mm_str+' '+ss_str+'s ('+ms_str+' ms)'
			day_tif(n)=STRING(yr_tif(n),FORMAT='(I4.4)')+ $
				STRING(mo_tif(n),FORMAT='(I2.2)')+STRING(dy_tif(n),FORMAT='(I2.2)')

			IF KEYWORD_SET(vbs) THEN BEGIN
    				PRINT,STRING(n,FORMAT='(I4.4)')+' th image: ', $
					yr_tif(n),mo_tif(n),dy_tif(n),hh_tif(n),mm_tif(n),ss_tif(n),ms_tif(n)
			ENDIF

; Put the data into array
                	img=READ_TIFF(found_files(i))
			IF KEYWORD_SET(abs) THEN BEGIN
				img_tif(n,*,*)=(ROTATE(img,7)-dk_fit)/at_fit
			ENDIF ELSE BEGIN
				img_tif(n,*,*)=ROTATE(img,7)
			ENDELSE

			n=n+1
		ENDFOR
	ENDFOR
	num_tif=n
	PRINT,''
	PRINT,' >>> The total number of stored images:',STRCOMPRESS(STRING(num_tif,FORMAT='(I4)'))
	PRINT,''

; Setting type and corresponding limit
	IF KEYWORD_SET(abs) THEN BEGIN
		typ_tif=1
		IF min_psa GT 1500 THEN set_scale_psa,0,1000
	ENDIF ELSE BEGIN
		typ_tif=0
		IF min_psa LT 1500 THEN set_scale_psa,2000,2500
	ENDELSE

; Setting current image
	go_psa_tif,0

	END

;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	FILE_PSA_FOV
;

	PRO file_psa_fov,camid,date,pxl=pxl,alt=alt

	COMMON prf_psa
	COMMON fov_psa

; Check keywords
	IF NOT KEYWORD_SET(pxl) THEN pxl=256
	IF NOT KEYWORD_SET(alt) THEN alt=110

; Message
	PRINT,''
	PRINT,'Read the FOV files for '+STRING(alt,FORMAT='(I3)')+' km altitude:'

; Selection of FOV files
	IF camid EQ 1 AND LONG(date) LT 20170126 THEN fname_base='C1_20161021_1901'
	IF camid EQ 1 AND LONG(date) GE 20170126 THEN fname_base='C1_20170208_0455'
	IF camid EQ 2 AND LONG(date) LT 20170130 THEN fname_base='C2_20161002_0025'
	IF camid EQ 2 AND LONG(date) GE 20170130 THEN fname_base='C2_20170201_0455'
	IF camid EQ 3 AND LONG(date) LT 20170126 THEN fname_base='C3_20161031_0440'
	IF camid EQ 3 AND LONG(date) GE 20170126 THEN fname_base='C3_20161031_0440'   ; This should be fixed soon
	IF camid EQ 4 AND LONG(date) LT 20170101 THEN fname_base='C4_20160927_2200'
	IF camid EQ 4 AND LONG(date) GE 20170101 THEN fname_base='C4_20170208_0455'
	IF camid EQ 5 AND LONG(date) LT 20170126 THEN fname_base='C5_20161021_1901'
	IF camid EQ 5 AND LONG(date) GE 20170126 THEN fname_base='C5_20170208_0455'
	IF camid EQ 6 THEN fname_base='C6_20170127_1714'
	IF camid EQ 7 THEN fname_base='C7_20170304-0527'
	
; Make arrays
	azm_psa=FLTARR(pxl,pxl)
	ele_psa=FLTARR(pxl,pxl)
	gla_psa=FLTARR(pxl,pxl)
	glo_psa=FLTARR(pxl,pxl)

; Read azimuth and elevation angles
	found_file=FILE_SEARCH(dir_psa+'fovd/azm_ele/'+fname_base+'*.txt',COUNT=num_found_file)
	PRINT,' ',found_file(0)
	PRINT,' ',found_file(1)
	OPENR,azm_lun,found_file(0),/GET_LUN
	OPENR,ele_lun,found_file(1),/GET_LUN
        FOR i=0,pxl-1 DO BEGIN
                tmp_read=FLTARR(pxl)
                READF,ele_lun,tmp_read
                ele_psa(*,pxl-1-i)=tmp_read
                READF,azm_lun,tmp_read
                azm_psa(*,pxl-1-i)=tmp_read
        ENDFOR
        FREE_LUN,azm_lun
        FREE_LUN,ele_lun

; Read glat and glon info
	found_file=FILE_SEARCH(dir_psa+'fovd/gla_glo/'+fname_base+'_'+STRING(alt,FORMAT='(I3.3)')+'*.txt', $
		COUNT=num_found_file)
	PRINT,' ',found_file(0)
	PRINT,' ',found_file(1)
	OPENR,gla_lun,found_file(0),/GET_LUN
	OPENR,glo_lun,found_file(1),/GET_LUN
        FOR i=0,pxl-1 DO BEGIN
                tmp_read=FLTARR(pxl)
                READF,gla_lun,tmp_read
                gla_psa(*,pxl-1-i)=tmp_read
                READF,glo_lun,tmp_read
                glo_psa(*,pxl-1-i)=tmp_read
        ENDFOR
        FREE_LUN,gla_lun
        FREE_LUN,glo_lun

	END

;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       GO_PSA_RAW
;
; EXAMPLE:
;
;       Go > go_psa_raw,10
;

        PRO go_psa_raw,img_to_go

        COMMON raw_psa

; Set new current image
        now_raw=0L
        now_raw=img_to_go

; Verbose
	hh_str=STRING(hh_raw(now_raw),FORMAT='(I2.2)')
	mm_str=STRING(mm_raw(now_raw),FORMAT='(I2.2)')
	ss_str=STRING(ss_raw(now_raw),FORMAT='(I2.2)')
	ms_str=STRING(ms_raw(now_raw),FORMAT='(I3.3)')
        current_hms=hh_str+mm_str+' '+ss_str+'s ('+ms_str+' ms)'
	current_num=STRING(now_raw,FORMAT='(I6.6)')
        PRINT,'Time of current RAW image: '+current_hms+' UT - '+current_num+' th image'

        END

;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       GO_PSA_TIF
;
; EXAMPLE:
;
;       Go > go_psa_tif,10
;

        PRO go_psa_tif,img_to_go

        COMMON tif_psa

; Set new current image
        now_tif=0L
        now_tif=img_to_go

; Verbose
	hh_str=STRING(hh_tif(now_tif),FORMAT='(I2.2)')
	mm_str=STRING(mm_tif(now_tif),FORMAT='(I2.2)')
	ss_str=STRING(ss_tif(now_tif),FORMAT='(I2.2)')
	ms_str=STRING(ms_tif(now_tif),FORMAT='(I3.3)')
        current_hms=hh_str+mm_str+' '+ss_str+'s ('+ms_str+' ms)'
	current_num=STRING(now_tif,FORMAT='(I6.6)')
        PRINT,'Time of current TIF image: '+current_hms+' UT - '+current_num+' th image'

        END

;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       GO_TIME_PSA
;
; EXAMPLE:
;
;       PsA > go_time_psa,025800
;

        PRO go_time_psa,image_time

        COMMON raw_psa
        COMMON tif_psa

; Identify total sec
	hh=FIX(STRMID(STRING(image_time,FORMAT='(I6.6)'),0,2))
	mm=FIX(STRMID(STRING(image_time,FORMAT='(I6.6)'),2,2))
	ss=FIX(STRMID(STRING(image_time,FORMAT='(I6.6)'),4,2))
	ts=3600L*hh+60L*mm+ss

; Find RAW image
	IF N_ELEMENTS(num_raw) NE 0 THEN BEGIN

        	found_image=WHERE(LONG(ts_raw) EQ ts,num_found_image)

; Jump image
        	IF num_found_image EQ 0 THEN BEGIN
                	PRINT,'Corresponding RAW image was not found!'
        	ENDIF ELSE BEGIN
                	go_psa_raw,found_image(0)
        	ENDELSE
	ENDIF ELSE BEGIN
                PRINT,'RAW images have not been filed!'
	ENDELSE

; Find TIF image
	IF N_ELEMENTS(num_tif) NE 0 THEN BEGIN

	        image_found=-1
        	FOR i=1,num_tif-1 DO BEGIN
                	IF ts_tif(i-1) LT ts AND ts_tif(i) GE ts THEN BEGIN

                        	image_found=i
                	ENDIF
        	ENDFOR
	
; Jump image
        	IF image_found EQ -1 THEN BEGIN
                	PRINT,'Corresponding RAW image was not found!'
        	ENDIF ELSE BEGIN
                	go_psa_tif,image_found
        	ENDELSE
	ENDIF ELSE BEGIN
                PRINT,'TIF images have not been filed!'
	ENDELSE

        END

;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	TIME_PSA
;
; EXAMPLE:
;
;       IDL > time_psa,030000,031000
;
; NOTE:
;
;       Time should be specified in 6-digits (hhmmss).
;

	PRO time_psa,stime,etime

        COMMON prf_psa

; Break down of stime
	hh_stime=stime/10000L
	mm_stime=(stime-hh_stime*10000L)/100L
	ss_stime=stime-hh_stime*10000L-mm_stime*100L

; Break down of etime
	hh_etime=etime/10000L
	mm_etime=(etime-hh_etime*10000L)/100L
	ss_etime=etime-hh_etime*10000L-mm_etime*100L

	stm_psa=3600.0*hh_stime+60.0*mm_stime+ss_stime
	etm_psa=3600.0*hh_etime+60.0*mm_etime+ss_etime

; Verbose
        PRINT,'Start time: '+STRMID(STRING(stime,FORMAT='(I6.6)'),0,4)+' '+ $
                STRMID(STRING(stime,FORMAT='(I6.6)'),4,2)+'s UT'
        PRINT,'End   time: '+STRMID(STRING(etime,FORMAT='(I6.6)'),0,4)+' '+ $
                STRMID(STRING(etime,FORMAT='(I6.6)'),4,2)+'s UT'

	END

;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	SET_SCALE_PSA
;
; EXAMPLE:
;
;	PsA > set_scale_psa,2200,2300
;

	PRO set_scale_psa,min_val,max_val,quiet=quiet

	COMMON prf_psa

	min_psa=min_val
	max_psa=max_val

	IF NOT KEYWORD_SET(quiet) THEN PRINT,'Scale limits: ',STRTRIM(FIX(min_psa),2),' to ',STRTRIM(FIX(max_psa),2)

	END

;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       PLOT_PSA_RAW_KEO
;
; EXAMPLE:
;
;	Read 1 min raw images from 2350 to 2351 UT on Oct 4, 2017 and plot the N-S and E-W keograms
;
;       PsA > file_psa_raw,2,'20161004',2350,1
;	PsA > time_psa,235000,235100
;       PsA > plot_psa_raw_keo
;
;	If you want to smooth the data with sliding window of 10 points:
;
;       PsA > file_psa_raw,2,'20161004',2330,1
;	PsA > time_psa,233000,233100
;       PsA > plot_psa_raw_keo,smo=10
;

        PRO plot_psa_raw_keo,x=x,y=y,smo=smo

	IF NOT KEYWORD_SET(x) THEN x=128
	IF NOT KEYWORD_SET(y) THEN y=128
        IF NOT KEYWORD_SET(smo) THEN smo=0

	window_psa,/landscape
	ERASE
	dpanel_psa,1,2,0,0
	plot_psa_raw_keo_panel,x=x,smo=smo
 	dpanel_psa,1,2,0,1
	plot_psa_raw_keo_panel,y=y,smo=smo

	END

;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       PLOT_PSA_RAW_KEO_PANEL
;
; NOTE:
;
;	IF you use a keyword x like x=128, N-S keogram is generated. If you use y=128, you have E-W keogram
;
; EXAMPLE:
;
;       PsA > window_psa,/landscape
;       PsA > dpanel_psa,1,2,0,0
;	PsA > plot_psa_raw_keo_panel,x=128
;       PsA > dpanel_psa,1,2,0,1
;	PsA > plot_psa_raw_keo_panel,y=128
;
; NOTE:
;
;	If you do not specify x and y, plain N-S keogram is plotted
;

        PRO plot_psa_raw_keo_panel,x=x,y=y,position=position,no_x=no_x,no_y=no_y,no_bar=no_bar,smo=smo

	COMMON prf_psa
        COMMON raw_psa

; Max value setting
        min_val=min_psa
        max_val=max_psa

; Plot area
	IF KEYWORD_SET(position) THEN !P.POSITION=position

; Determine time spacing and UT axis
        ts=stm_psa/3600.0 & te=etm_psa/3600.0
        time_axis_intervals,ts,te,minor_ticks,where_tick,no_of_ticks

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
        ENDIF ELSE BEGIN
		ytit='South to North Keogram!C!DCut at X = 128'
		IF KEYWORD_SET(x) THEN BEGIN
			ytit='South to North Keogram!C!DCut at X = '+STRING(x,FORMAT='(I3)')
			ytickname=['S',' ','Z',' ','N']
		ENDIF
		IF KEYWORD_SET(y) THEN BEGIN
			ytit='East to West Keogram!C!DCut at Y = '+STRING(y,FORMAT='(I3)')
			ytickname=['E',' ','Z',' ','W']
		ENDIF
	ENDELSE

; Plot frame
        yticks=4 & yminor=4
        pos=!P.POSITION
        PLOT,[0],[0],XRANGE=[ts,te],/XSTYLE,YRANGE=[0,256],/YSTYLE,XTITLE=xtit,YTITLE=ytit, $
                XTICKFORMAT=xtickformat,YTICKNAME=ytickname,YTICKS=yticks,YMINOR=yminor,XTICKS=no_of_ticks, $
                XTICKV=where_tick,XMINOR=minor_ticks,YTICKFORMAT=ytickformat
        
; Selecting data for contouring
	plot_inx=WHERE(ts_raw GE stm_psa AND ts_raw LE etm_psa,num_plot_inx)
        keo_plot=REFORM(img_raw(plot_inx,128,*))
	IF KEYWORD_SET(x) THEN keo_plot=REFORM(img_raw(plot_inx,x,*))
	IF KEYWORD_SET(y) THEN keo_plot=REFORM(img_raw(plot_inx,*,y))
	IF KEYWORD_SET(smo) THEN FOR i=0,255 DO keo_plot(*,i)=SMOOTH(keo_plot(*,i),smo)
        mat_x=ts_raw(plot_inx)/3600.0
        mat_y=FINDGEN(256)+0.5
        min_vals=WHERE(keo_plot LT min_val,num_min_vals)
        IF num_min_vals NE 0 THEN keo_plot(min_vals)=min_val
        max_vals=WHERE(keo_plot GT max_val,num_max_vals)
        IF num_max_vals NE 0 THEN keo_plot(max_vals)=max_val

; Plot the data in contour form
        contour_lines=30
        levels_set=(max_val-min_val)*FINDGEN(contour_lines)/FLOAT(contour_lines)+min_val
        colour_set=!D.TABLE_SIZE*INDGEN(contour_lines)/contour_lines
        CONTOUR,keo_plot,mat_x,mat_y,/OVERPLOT,/FILL,LEVELS=levels_set,C_COLORS=colour_set

; Overplot frame
	OPLOT,[0,72],[128,128],LINE=1
        PLOT,[0],[0],XRANGE=[ts,te],/XSTYLE,YRANGE=[0,256],/YSTYLE,XTITLE=xtit,YTITLE=ytit, $
                XTICKFORMAT=xtickformat,YTICKNAME=ytickname,YTICKS=yticks,YMINOR=yminor,XTICKS=no_of_ticks, $
                XTICKV=where_tick,XMINOR=minor_ticks,YTICKFORMAT=ytickformat

; Plot some messages
	pos=!P.POSITION
	XYOUTS,pos(0),pos(3)+0.01,day_raw,/NORMAL
	mesg='Cam'+STRING(cid_raw,FORMAT='(I1.1)')+' in '+loc_raw+' ('+flt_raw+')'
	XYOUTS,pos(2),pos(3)+0.01,mesg,/NORMAL,ALI=1

; Plot colour bar
        IF NOT KEYWORD_SET(no_bar) THEN BEGIN
		position=[!P.POSITION(2)+0.02,!P.POSITION(1),!P.POSITION(2)+0.03,!P.POSITION(3)]
		IF typ_raw EQ 1 THEN legend='Absolute Optical Intensity (R)' ELSE legend='Raw Count'
		plot_psa_colour_bar,position=position,legend=legend
	ENDIF

        END

;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       PLOT_PSA_TIF_KEO
;
; EXAMPLE:
;
;       Read 1 hour tif images from 2300 to 2400 UT on Oct 4, 2017 and plot the N-S and E-W keograms
;
;       PsA > file_psa_raw,2,'20161004',23,1
;       PsA > time_psa,230000,240000
;       PsA > plot_psa_tif_keo
;
;       If you want to smooth the data with sliding window of 10 points:
;
;       PsA > plot_psa_tif_keo,smo=10
;

        PRO plot_psa_tif_keo,x=x,y=y,smo=smo

; Keyword handling
        IF NOT KEYWORD_SET(x) THEN x=128
        IF NOT KEYWORD_SET(y) THEN y=128
        IF NOT KEYWORD_SET(smo) THEN smo=0

; Make window
        ;window_psa,/landscape
        ;ERASE

; Plot N-S keogram
        dpanel_psa,1,2,0,0
        plot_psa_tif_keo_panel,x=x,smo=smo

; Plot E-W keogram
        dpanel_psa,1,2,0,1
        plot_psa_tif_keo_panel,y=y,smo=smo

        END

;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       PLOT_PSA_TIF_KEO_PANEL
;
; EXAMPLE:
;
;       PsA > window_psa,/landscape
;       PsA > dpanel_psa,1,2,0,0
;	PsA > plot_psa_tif_keo_panel,x=128
;       PsA > dpanel_psa,1,2,0,1
;	PsA > plot_psa_tif_keo_panel,y=128
;

        PRO plot_psa_tif_keo_panel,x=x,y=y,position=position,no_x=no_x,no_y=no_y,no_bar=no_bar,smo=smo

	COMMON prf_psa
        COMMON tif_psa

; Max value setting
        min_val=min_psa
        max_val=max_psa

; Plot area
	IF KEYWORD_SET(position) THEN !P.POSITION=position

; Determine time spacing and UT axis
        ts=stm_psa/3600.0 & te=etm_psa/3600.0
        time_axis_intervals,ts,te,minor_ticks,where_tick,no_of_ticks

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
        ENDIF ELSE BEGIN
		ytit='South to North Keogram!C!DCut at X = 128'
		IF KEYWORD_SET(x) THEN BEGIN
			ytit='South to North Keogram!C!DCut at X = '+STRING(x,FORMAT='(I3)')
			ytickname=['S',' ','Z',' ','N']
		ENDIF
		IF KEYWORD_SET(y) THEN BEGIN
			ytit='East to West Keogram!C!DCut at Y = '+STRING(y,FORMAT='(I3)')
			ytickname=['E',' ','Z',' ','W']
		ENDIF
	ENDELSE

; Plot frame
        yticks=4 & yminor=4
        pos=!P.POSITION
        PLOT,[0],[0],XRANGE=[ts,te],/XSTYLE,YRANGE=[0,256],/YSTYLE,XTITLE=xtit,YTITLE=ytit, $
                XTICKFORMAT=xtickformat,YTICKNAME=ytickname,YTICKS=yticks,YMINOR=yminor,XTICKS=no_of_ticks, $
                XTICKV=where_tick,XMINOR=minor_ticks,YTICKFORMAT=ytickformat
        
; Selecting data for contouring
	plot_inx=WHERE(ts_tif GE stm_psa AND ts_tif LE etm_psa,num_plot_inx)
        keo_plot=REFORM(img_tif(plot_inx,128,*))
	IF KEYWORD_SET(x) THEN keo_plot=REFORM(img_tif(plot_inx,x,*))
	IF KEYWORD_SET(y) THEN keo_plot=REFORM(img_tif(plot_inx,*,y))
        IF KEYWORD_SET(smo) THEN FOR i=0,255 DO keo_plot(*,i)=SMOOTH(keo_plot(*,i),smo)
        mat_x=ts_tif(plot_inx)/3600.0
        mat_y=FINDGEN(256)+0.5
        min_vals=WHERE(keo_plot LT min_val,num_min_vals)
        IF num_min_vals NE 0 THEN keo_plot(min_vals)=min_val
        max_vals=WHERE(keo_plot GT max_val,num_max_vals)
        IF num_max_vals NE 0 THEN keo_plot(max_vals)=max_val

; Plot the data in contour form
        contour_lines=30
        levels_set=(max_val-min_val)*FINDGEN(contour_lines)/FLOAT(contour_lines)+min_val
        colour_set=!D.TABLE_SIZE*INDGEN(contour_lines)/contour_lines
        CONTOUR,keo_plot,mat_x,mat_y,/OVERPLOT,/FILL,LEVELS=levels_set,C_COLORS=colour_set

; Overplot frame
	OPLOT,[0,72],[128,128],LINE=1
        PLOT,[0],[0],XRANGE=[ts,te],/XSTYLE,YRANGE=[0,256],/YSTYLE,XTITLE=xtit,YTITLE=ytit, $
                XTICKFORMAT=xtickformat,YTICKNAME=ytickname,YTICKS=yticks,YMINOR=yminor,XTICKS=no_of_ticks, $
                XTICKV=where_tick,XMINOR=minor_ticks,YTICKFORMAT=ytickformat

; Plot some messages
	pos=!P.POSITION
	XYOUTS,pos(0),pos(3)+0.01,day_tif,/NORMAL
        mesg='Cam'+STRING(cid_tif,FORMAT='(I1.1)')+' in '+loc_tif+' ('+flt_tif+')'
        XYOUTS,pos(2),pos(3)+0.01,mesg,/NORMAL,ALI=1

; Plot colour bar
        IF NOT KEYWORD_SET(no_bar) THEN BEGIN
		position=[!P.POSITION(2)+0.02,!P.POSITION(1),!P.POSITION(2)+0.03,!P.POSITION(3)]
		IF typ_tif EQ 1 THEN legend='Absolute Optical Intensity (R)' ELSE legend='Raw Count'
		plot_psa_colour_bar,position=position,legend=legend
	ENDIF

        END

;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       PLOT_PSA_RAW_LINE
;
; EXAMPLE:
;
;       Read 1 min raw images from 2350 to 2351 UT on Oct 4, 2017 and plot the optical intensity at the central pixel:
;
;       PsA > file_psa_raw,2,'20161004',2350,1
;       PsA > time_psa,235000,235100
;       PsA > plot_psa_raw_line
;
;       If you want to smooth the data with sliding window of 10 points:
;
;       PsA > file_psa_raw,2,'20161004',2330,1
;       PsA > time_psa,233000,233100
;       PsA > plot_psa_raw_line,smo=10
;

        PRO plot_psa_raw_line,x=x,y=y,smo=smo

        IF NOT KEYWORD_SET(x) THEN x=128
        IF NOT KEYWORD_SET(y) THEN y=128
        IF NOT KEYWORD_SET(smo) THEN smo=0

        window_psa,/landscape
	ERASE
        dpanel_psa,1,1,0,0
        plot_psa_raw_line_panel,x=x,y=y,smo=smo

        END

;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       PLOT_PSA_RAW_LINE_PANEL
;
; EXAMPLE:
;
;       PsA > window_psa,/landscape
;       PsA > dpanel_psa,1,1,0,0
;	PsA > plot_psa_raw_line_panel,x=108,y=200
;

        PRO plot_psa_raw_line_panel,x=x,y=y,position=position,no_x=no_x,no_y=no_y,smo=smo,col=col

	COMMON prf_psa
        COMMON raw_psa

; Max value setting
        min_val=min_psa
        max_val=max_psa

; Plot area
	IF KEYWORD_SET(position) THEN !P.POSITION=position

; Which point to be plotted
	IF NOT KEYWORD_SET(x) THEN x=128
	IF NOT KEYWORD_SET(y) THEN y=128

; Determine time spacing and UT axis
        ts=stm_psa/3600.0 & te=etm_psa/3600.0
        time_axis_intervals,ts,te,minor_ticks,where_tick,no_of_ticks

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
        ENDIF ELSE BEGIN
		IF typ_raw EQ 0 THEN ytit='Raw Count'
		IF typ_raw EQ 1 THEN ytit='Absolute Optical Intensity (R)'
	ENDELSE

; Plot frame
        yticks=4 & yminor=4
        pos=!P.POSITION
        PLOT,[0],[0],XRANGE=[ts,te],/XSTYLE,YRANGE=[min_val,max_val],/YSTYLE,XTITLE=xtit,YTITLE=ytit, $
                XTICKFORMAT=xtickformat,YTICKNAME=ytickname,YTICKS=yticks,YMINOR=yminor,XTICKS=no_of_ticks, $
                XTICKV=where_tick,XMINOR=minor_ticks,YTICKFORMAT=ytickformat
        
; Selecting data for contouring
	plot_inx=WHERE(ts_raw GE stm_psa AND ts_raw LE etm_psa,num_plot_inx)
        pnt_plot=REFORM(img_raw(plot_inx,x,y))
	IF KEYWORD_SET(smo) THEN pnt_plot=SMOOTH(pnt_plot,smo)
        tim_plot=ts_raw(plot_inx)/3600.0
	OPLOT,tim_plot,pnt_plot,COL=col
	SAVE,pnt_plot,FILENAME='Cam'+STRING(cid_raw,FORMAT='(I1.1)')+'.sav'

; Overplot frame
        PLOT,[0],[0],XRANGE=[ts,te],/XSTYLE,YRANGE=[min_val,max_val],/YSTYLE,XTITLE=xtit,YTITLE=ytit, $
                XTICKFORMAT=xtickformat,YTICKNAME=ytickname,YTICKS=yticks,YMINOR=yminor,XTICKS=no_of_ticks, $
                XTICKV=where_tick,XMINOR=minor_ticks,YTICKFORMAT=ytickformat

; Plot some messages
	pos=!P.POSITION
	XYOUTS,pos(0),pos(3)+0.01,day_raw,/NORMAL
	mesg='Cam'+STRING(cid_raw,FORMAT='(I1.1)')+' in '+loc_raw+' ('+flt_raw+')'
	XYOUTS,pos(2),pos(3)+0.01,mesg,/NORMAL,ALI=1
	
; Show point to be plotted
	mesg='X = '+STRING(x,FORMAT='(I3)')+', Y = '+STRING(y,FORMAT='(I3)')
	XYOUTS,0.5*(pos(0)+pos(2)),pos(1)+0.9*(pos(3)-pos(1)),mesg,/NORMAL,ALI=0.5

        END

;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       PLOT_PSA_TIF_LINE
;
; EXAMPLE:
;
;       Read 1 hour tif images from 2300 to 2400 UT on Oct 4, 2017 and plot the optical intensity at the central pixel:
;
;       PsA > file_psa_tif,2,'20161004',23,1
;       PsA > time_psa,230000,240000
;       PsA > plot_psa_tif_line
;
;       If you want to smooth the data with sliding window of 10 points:
;
;       PsA > plot_psa_tif_line,smo=10
;

        PRO plot_psa_tif_line,x=x,y=y,smo=smo

        IF NOT KEYWORD_SET(x) THEN x=128
        IF NOT KEYWORD_SET(y) THEN y=128
        IF NOT KEYWORD_SET(smo) THEN smo=0

; Make window
        window_psa,/landscape
        ERASE

; Plot the data
        dpanel_psa,1,1,0,0
        plot_psa_tif_line_panel,x=x,y=y,smo=smo

        END

;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       PLOT_PSA_TIF_LINE_PANEL
;
; EXAMPLE:
;
;       PsA > window_psa,/landscape
;       PsA > dpanel_psa,1,1,0,0
;	PsA > plot_psa_tif_line_panel,x=108,y=200
;

        PRO plot_psa_tif_line_panel,x=x,y=y,position=position,no_x=no_x,no_y=no_y,smo=smo,col=col

	COMMON prf_psa
        COMMON tif_psa

; Max value setting
        min_val=min_psa
        max_val=max_psa

; Plot area
	IF KEYWORD_SET(position) THEN !P.POSITION=position

; Which point to be plotted
	IF NOT KEYWORD_SET(x) THEN x=128
	IF NOT KEYWORD_SET(y) THEN y=128

; Determine time spacing and UT axis
        ts=stm_psa/3600.0 & te=etm_psa/3600.0
        time_axis_intervals,ts,te,minor_ticks,where_tick,no_of_ticks

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
        ENDIF ELSE BEGIN
		IF typ_tif EQ 0 THEN ytit='Raw Count'
		IF typ_tif EQ 1 THEN ytit='Absolute Optical Intensity (R)'
	ENDELSE

; Plot frame
        yticks=4 & yminor=4
        pos=!P.POSITION
        PLOT,[0],[0],XRANGE=[ts,te],/XSTYLE,YRANGE=[min_val,max_val],/YSTYLE,XTITLE=xtit,YTITLE=ytit, $
                XTICKFORMAT=xtickformat,YTICKNAME=ytickname,YTICKS=yticks,YMINOR=yminor,XTICKS=no_of_ticks, $
                XTICKV=where_tick,XMINOR=minor_ticks,YTICKFORMAT=ytickformat
        
; Selecting data for contouring
	plot_inx=WHERE(ts_tif GE stm_psa AND ts_tif LE etm_psa,num_plot_inx)
        pnt_plot=REFORM(img_tif(plot_inx,x,y))
	IF KEYWORD_SET(smo) THEN pnt_plot=SMOOTH(pnt_plot,smo)
        tim_plot=ts_tif(plot_inx)/3600.0
	OPLOT,tim_plot,pnt_plot,COL=col

; Overplot frame
        PLOT,[0],[0],XRANGE=[ts,te],/XSTYLE,YRANGE=[min_val,max_val],/YSTYLE,XTITLE=xtit,YTITLE=ytit, $
                XTICKFORMAT=xtickformat,YTICKNAME=ytickname,YTICKS=yticks,YMINOR=yminor,XTICKS=no_of_ticks, $
                XTICKV=where_tick,XMINOR=minor_ticks,YTICKFORMAT=ytickformat

; Plot some messages
	pos=!P.POSITION
	XYOUTS,pos(0),pos(3)+0.01,day_tif,/NORMAL
        mesg='Cam'+STRING(cid_raw,FORMAT='(I1.1)')+' in '+loc_raw+' ('+flt_raw+')'
        XYOUTS,pos(2),pos(3)+0.01,mesg,/NORMAL,ALI=1
	
; Show point to be plotted
	mesg='X = '+STRING(x,FORMAT='(I3)')+', Y = '+STRING(y,FORMAT='(I3)')
	XYOUTS,0.5*(pos(0)+pos(2)),pos(1)+0.9*(pos(3)-pos(1)),mesg,/NORMAL,ALI=0.5

        END

;--------------------------------------------------------------------------------------------------------------------------
; NAME
;
;	PLOT_PSA_RAW
;
; EXAMPLE:
;
;       PsA > window_psa,/landscape
;	PsA > file_psa,1,'20161004',2350,1
;	PsA > go_time_psa,235010
;       PsA > plot_psa_raw,2,3
;

	PRO plot_psa_raw,xmaps,ymaps,smo=smo,skip=skip

	COMMON raw_psa

; Panel positioning
        IF N_PARAMS() NE 2 THEN BEGIN
                xmaps=1 & ymaps=1
        ENDIF

; Skip setting
	IF NOT KEYWORD_SET(skip) THEN skip=1

; Make window
	window_psa,/landscape
	ERASE

; Define panels and plot the data
        FOR i=0,xmaps*ymaps-1 DO BEGIN

; End of file identification
                IF now_raw GT num_raw-1 THEN BEGIN
                        PRINT,'End of images reached'
                        RETURN
                ENDIF

; Panel management
                y_div=i/xmaps & x_div=i-xmaps*(i/xmaps)

; Plot the panel
		dpanel_psa,xmaps,ymaps,x_div,y_div,/square
		IF x_div EQ xmaps-1 THEN no_bar=0 ELSE no_bar=1
                plot_psa_raw_panel,smo=smo,no_bar=no_bar

; Progress the scan
                now_raw=now_raw+skip
        ENDFOR

	END

;--------------------------------------------------------------------------------------------------------------------------
; NAME
;
;	PLOT_PSA_RAW_PANEL
;
; EXAMPLE:
;
;       PsA > window_psa,/landscape
;	PsA > file_psa,1,'20161004',2350,1
;	PsA > go_time_psa,235010
;       PsA > dpanel_psa,2,1,0,0,/square
;	PsA > plot_psa_raw_panel
;       PsA > dpanel_psa,2,1,1,0,/square
;	PsA > plot_psa_raw_panel,smo=10
;

	PRO plot_psa_raw_panel,position=position,smo=smo,no_bar=no_bar

	COMMON prf_psa
        COMMON raw_psa

; Plot area
	IF KEYWORD_SET(position) THEN !P.POSITION=position
	pos=!P.POSITION

; Plot the data
	;tmpimg1=SMOOTH(REFORM(img_raw(now_raw,*,*)),smo)
	IF KEYWORD_SET(smo) THEN BEGIN
		tmp_sum=FLTARR(256,256)
		FOR i=0,smo-1 DO BEGIN
			tmp_sum=tmp_sum+REFORM(img_raw(now_raw+i,*,*))
		ENDFOR
		tmpimg1=tmp_sum/smo
	ENDIF ELSE BEGIN
		tmpimg1=REFORM(img_raw(now_raw,*,*))
	ENDELSE
        tmpimg2=CONGRID(tmpimg1,FIX((pos(2)-pos(0))*!D.X_Size),FIX((pos(3)-pos(1))*!D.Y_Size))
       	tmpimg3=BYTSCL(tmpimg2,MAX=max_psa,MIN=min_psa,TOP=!D.TABLE_SIZE)
       	TV,tmpimg3,FIX(pos(0)*!D.X_Size),FIX(pos(1)*!D.Y_Size), $
               	XSIZE=FIX((pos(2)-pos(0))*!D.X_Size),YSIZE=FIX((pos(3)-pos(1))*!D.Y_Size)
	PLOT,[0],[0],XRANGE=[0,256],YRANGE=[0,256],/XSTYLE,/YSTYLE,XTICKS=4,YTICKS=4, $
		XTICKFORMAT='no_ticks',YTICKFORMAT='no_ticks'
	OPLOT,[128,128],[0,256],LINE=1
	OPLOT,[0,256],[128,128],LINE=1

; Mask no data area
        FOR i=0,180 DO BEGIN
                x1=128+128*COS(2*(i-1)*!PI/180.0)
                x2=128+128*COS(2*(i+1)*!PI/180.0)
                x3=128+198*COS(2*(i+1)*!PI/180.0)
                x4=128+198*COS(2*(i-1)*!PI/180.0)
                y1=128+128*SIN(2*(i-1)*!PI/180.0)
                y2=128+128*SIN(2*(i+1)*!PI/180.0)
                y3=128+198*SIN(2*(i+1)*!PI/180.0)
                y4=128+198*SIN(2*(i-1)*!PI/180.0)
                POLYFILL,[x1,x2,x3,x4],[y1,y2,y3,y4],COL=100,NOCLIP=0
        ENDFOR

; Plot frame once again
	PLOT,[0],[0],XRANGE=[0,256],YRANGE=[0,256],/XSTYLE,/YSTYLE,XTICKS=4,YTICKS=4, $
		XTICKFORMAT='no_ticks',YTICKFORMAT='no_ticks'

; Plot messages
	pos=!P.POSITION
	XYOUTS,pos(0),pos(3)+0.01,day_raw(now_raw),/NORMAL
	XYOUTS,pos(2),pos(3)+0.01,ti_raw(now_raw),/NORMAL,ALI=1
	mesg='Cam'+STRING(cid_raw,FORMAT='(I1.1)')+' in '+loc_raw+' ('+flt_raw+')'
	;XYOUTS,pos(0)-0.01,pos(3),'Cam'+STRING(cid_raw,FORMAT='(I1.1)'),/NORMAL,ALI=1,ORI=90
	XYOUTS,pos(0)-0.01,pos(3),mesg,/NORMAL,ALI=1,ORI=90

; Plot colour bar
        IF NOT KEYWORD_SET(no_bar) THEN BEGIN
		position=[!P.POSITION(2)+0.01,!P.POSITION(1),!P.POSITION(2)+0.02,!P.POSITION(3)]
		IF typ_raw EQ 1 THEN legend='Absolute Optical Intensity (R)' ELSE legend='Raw Count'
		plot_psa_colour_bar,position=position,legend=legend
	ENDIF

	END

;--------------------------------------------------------------------------------------------------------------------------
; NAME
;
;	PLOT_PSA_TIF
;
; EXAMPLE:
;
;       PsA > window_psa,/landscape
;	PsA > file_psa_tif,1,'20161004',23,1
;	PsA > go_time_psa,233000
;       PsA > plot_psa_tif,2,3
;

	PRO plot_psa_tif,xmaps,ymaps,smo=smo,skip=skip

	COMMON tif_psa

; Panel positioning
        IF N_PARAMS() NE 2 THEN BEGIN
                xmaps=1 & ymaps=1
        ENDIF

; Skip setting
	IF NOT KEYWORD_SET(skip) THEN skip=1

; Make window
	;window_psa,/landscape
	;ERASE

; Define panels and plot the data
        FOR i=0,xmaps*ymaps-1 DO BEGIN

; End of file identification
                IF now_tif GT num_tif-1 THEN BEGIN
                        PRINT,'End of images reached'
                        RETURN
                ENDIF

; Panel management
                y_div=i/xmaps & x_div=i-xmaps*(i/xmaps)

; Plot the panel
		dpanel_psa,xmaps,ymaps,x_div,y_div,/square
		IF x_div EQ xmaps-1 THEN no_bar=0 ELSE no_bar=1
                plot_psa_tif_panel,smo=smo,no_bar=no_bar

; Progress the scan
                now_tif=now_tif+skip
        ENDFOR

	END

;--------------------------------------------------------------------------------------------------------------------------
; NAME
;
;	PLOT_PSA_TIF_PANEL
;

	PRO plot_psa_tif_panel,position=position,smo=smo,no_bar=no_bar

	COMMON prf_psa
        COMMON tif_psa

; Plot area
	IF KEYWORD_SET(position) THEN !P.POSITION=position
	pos=!P.POSITION

; Plot the data
	;tmpimg1=SMOOTH(REFORM(img_tif(now_tif,*,*)),smo)
	IF KEYWORD_SET(smo) THEN BEGIN
		tmp_sum=FLTARR(256,256)
		FOR i=0,smo-1 DO BEGIN
			tmp_sum=tmp_sum+REFORM(img_tif(now_tif+i,*,*))
		ENDFOR
		tmpimg1=tmp_sum/smo
	ENDIF ELSE BEGIN
		tmpimg1=REFORM(img_tif(now_tif,*,*))
	ENDELSE
        tmpimg2=CONGRID(tmpimg1,FIX((pos(2)-pos(0))*!D.X_Size),FIX((pos(3)-pos(1))*!D.Y_Size))
       	tmpimg3=BYTSCL(tmpimg2,MAX=max_psa,MIN=min_psa,TOP=!D.TABLE_SIZE)
       	TV,tmpimg3,FIX(pos(0)*!D.X_Size),FIX(pos(1)*!D.Y_Size), $
               	XSIZE=FIX((pos(2)-pos(0))*!D.X_Size),YSIZE=FIX((pos(3)-pos(1))*!D.Y_Size)
	PLOT,[0],[0],XRANGE=[0,256],YRANGE=[0,256],/XSTYLE,/YSTYLE,XTICKS=4,YTICKS=4, $
		XTICKFORMAT='no_ticks',YTICKFORMAT='no_ticks'
	OPLOT,[128,128],[0,256],LINE=1
	OPLOT,[0,256],[128,128],LINE=1

; Mask no data area
        FOR i=0,180 DO BEGIN
                x1=128+128*COS(2*(i-1)*!PI/180.0)
                x2=128+128*COS(2*(i+1)*!PI/180.0)
                x3=128+198*COS(2*(i+1)*!PI/180.0)
                x4=128+198*COS(2*(i-1)*!PI/180.0)
                y1=128+128*SIN(2*(i-1)*!PI/180.0)
                y2=128+128*SIN(2*(i+1)*!PI/180.0)
                y3=128+198*SIN(2*(i+1)*!PI/180.0)
                y4=128+198*SIN(2*(i-1)*!PI/180.0)
                POLYFILL,[x1,x2,x3,x4],[y1,y2,y3,y4],COL=100,NOCLIP=0
        ENDFOR

; Plot frame once again
	PLOT,[0],[0],XRANGE=[0,256],YRANGE=[0,256],/XSTYLE,/YSTYLE,XTICKS=4,YTICKS=4, $
		XTICKFORMAT='no_ticks',YTICKFORMAT='no_ticks'

; Plot messages
	pos=!P.POSITION
	XYOUTS,pos(0),pos(3)+0.01,day_tif(now_tif),/NORMAL
	XYOUTS,pos(2),pos(3)+0.01,ti_tif(now_tif),/NORMAL,ALI=1
        mesg='Cam'+STRING(cid_tif,FORMAT='(I1.1)')+' in '+loc_tif+' ('+flt_tif+')'
        XYOUTS,pos(0)-0.01,pos(3),mesg,/NORMAL,ALI=1,ORI=90

; Plot colour bar
        IF NOT KEYWORD_SET(no_bar) THEN BEGIN
		position=[!P.POSITION(2)+0.01,!P.POSITION(1),!P.POSITION(2)+0.02,!P.POSITION(3)]
		IF typ_tif EQ 1 THEN legend='Absolute Optical Intensity (R)' ELSE legend='Raw Count'
		plot_psa_colour_bar,position=position,legend=legend
	ENDIF

	END

;--------------------------------------------------------------------------------------------------------------------------
; NAME
;
;	PLOT_PSA_RAW_MAP
;
; EXAMPLE:
;
;       PsA > window_psa,/landscape
;	PsA > file_psa,1,'20161004',2350,1
;	PsA > go_time_psa,235010
;       PsA > plot_psa_raw_map,2,3
;

	PRO plot_psa_raw_map,xmaps,ymaps,smo=smo,skip=skip

	COMMON raw_psa

; Panel positioning
        IF N_PARAMS() NE 2 THEN BEGIN
                xmaps=1 & ymaps=1
        ENDIF

; Skip setting
	IF NOT KEYWORD_SET(skip) THEN skip=1

; Make window
	window_psa,/landscape
	ERASE

; Define panels and plot the data
        FOR i=0,xmaps*ymaps-1 DO BEGIN

; End of file identification
                IF now_raw GT num_raw-1 THEN BEGIN
                        PRINT,'End of images reached'
                        RETURN
                ENDIF

; Panel management
                y_div=i/xmaps & x_div=i-xmaps*(i/xmaps)

; Plot the panel
		dpanel_psa,xmaps,ymaps,x_div,y_div,/square
		IF x_div EQ xmaps-1 THEN no_bar=0 ELSE no_bar=1
                plot_psa_raw_map_panel,smo=smo,no_bar=no_bar

; Progress the scan
                now_raw=now_raw+skip
        ENDFOR

	END

;--------------------------------------------------------------------------------------------------------------------------
; NAME
;
;	PLOT_PSA_RAW_MAP_PANEL
;

	PRO plot_psa_raw_map_panel,position=position,smo=smo,no_bar=no_bar,ele_lim=ele_lim,xrange=xrange,yrange=yrange

	COMMON prf_psa
	COMMON fov_psa
        COMMON raw_psa

; Plot area
	IF KEYWORD_SET(position) THEN !P.POSITION=position
	pos=!P.POSITION

; Set up the elevation angle limitation
	IF NOT KEYWORD_SET(ele_lim) THEN ele_lim=10.0

; Set up the plotting area
	IF NOT KEYWORD_SET(xrange) THEN BEGIN
		xrange=[4,34]
		IF cid_raw EQ 2 THEN xrange=[12,41] ; offset is 29
	ENDIF
	IF NOT KEYWORD_SET(yrange) THEN BEGIN
		yrange=[64,75]
		IF cid_raw EQ 2 THEN yrange=[62,73] ; offset is 11
	ENDIF

; Plot frame
        PLOT,[0],[0],XRANGE=xrange,YRANGE=yrange,/XSTYLE,/YSTYLE,XTITLE='Geo. Lon.',YTITLE='Geo. Lat'

; Plot the data
	IF KEYWORD_SET(smo) THEN BEGIN
		tmp_sum=FLTARR(256,256)
		FOR i=0,smo-1 DO tmp_sum=tmp_sum+REFORM(img_raw(now_raw+i,*,*))
		tmp_img=tmp_sum/smo
	ENDIF ELSE BEGIN
		tmp_img=REFORM(img_raw(now_raw,*,*))
	ENDELSE

; Prepare array for contouring
        tmp_bin=WHERE(gla_psa NE 99999.9 AND gla_psa NE 99999.9 AND ele_psa GE ele_lim,no_tmp_bin)
        x_cont=MAKE_ARRAY(no_tmp_bin,/FLOAT,VALUE=0.0)
        y_cont=MAKE_ARRAY(no_tmp_bin,/FLOAT,VALUE=0.0)
        z_cont=MAKE_ARRAY(no_tmp_bin,/FLOAT,VALUE=0.0)

; Put values into the arrays for contouring
        num_cnt=0L
        FOR x=0,255 DO BEGIN
                FOR y=0,255 DO BEGIN
                        IF gla_psa(x,y) NE 99999.9 AND glo_psa(x,y) NE 99999.9 AND ele_psa(x,y) GE ele_lim THEN BEGIN

                                z_cont(num_cnt)=FLOAT(tmp_img(x,y));*COS(zan_psa(x,y)*!DTOR)
                                IF z_cont(num_cnt) GT max_psa THEN z_cont(num_cnt)=max_psa
                                IF z_cont(num_cnt) LT min_psa THEN z_cont(num_cnt)=min_psa

                                x_cont(num_cnt)=glo_psa(x,y)
                                y_cont(num_cnt)=gla_psa(x,y)

                                num_cnt++
                        ENDIF
                ENDFOR
        ENDFOR

; Contouring
        contour_lines=30
        levels_set=(max_psa-min_psa)*FINDGEN(contour_lines)/FLOAT(contour_lines)+min_psa
        colour_set=!D.TABLE_SIZE*INDGEN(contour_lines)/contour_lines
        CONTOUR,z_cont,x_cont,y_cont,/OVERPLOT,/FILL,LEVELS=levels_set,C_COLORS=colour_set,/IRREGULAR

; Plot coastlines
	overlay_coast

; Plot messages
	pos=!P.POSITION
	XYOUTS,pos(0),pos(3)+0.01,day_raw(now_raw)+': Cam'+STRING(cid_raw,FORMAT='(I1.1)'),/NORMAL
	XYOUTS,pos(2),pos(3)+0.01,ti_raw(now_raw),/NORMAL,ALI=1

; Plot colour bar
        IF NOT KEYWORD_SET(no_bar) THEN BEGIN
		position=[!P.POSITION(2)+0.01,!P.POSITION(1),!P.POSITION(2)+0.02,!P.POSITION(3)]
		IF typ_raw EQ 1 THEN legend='Absolute Optical Intensity (R)' ELSE legend='Raw Count'
		plot_psa_colour_bar,position=position,legend=legend
	ENDIF

	END

;--------------------------------------------------------------------------------------------------------------------------
; NAME
;
;	PLOT_PSA_TIF_MAP
;
; EXAMPLE:
;
;       PsA > window_psa,/landscape
;	PsA > file_psa_tif,1,'20161004',23,1
;	PsA > go_time_psa,235010
;       PsA > plot_psa_tif_map,2,3
;

	PRO plot_psa_tif_map,xmaps,ymaps,smo=smo,skip=skip

	COMMON tif_psa

; Panel positioning
        IF N_PARAMS() NE 2 THEN BEGIN
                xmaps=1 & ymaps=1
        ENDIF

; Skip setting
	IF NOT KEYWORD_SET(skip) THEN skip=1

; Make window
	window_psa,/landscape
	ERASE

; Define panels and plot the data
        FOR i=0,xmaps*ymaps-1 DO BEGIN

; End of file identification
                IF now_tif GT num_tif-1 THEN BEGIN
                        PRINT,'End of images reached'
                        RETURN
                ENDIF

; Panel management
                y_div=i/xmaps & x_div=i-xmaps*(i/xmaps)

; Plot the panel
		dpanel_psa,xmaps,ymaps,x_div,y_div,/square
		IF x_div EQ xmaps-1 THEN no_bar=0 ELSE no_bar=1
                plot_psa_tif_map_panel,smo=smo,no_bar=no_bar

; Progress the scan
                now_tif=now_tif+skip
        ENDFOR

	END

;--------------------------------------------------------------------------------------------------------------------------
; NAME
;
;	PLOT_PSA_TIF_MAP_PANEL
;

	PRO plot_psa_tif_map_panel,position=position,smo=smo,no_bar=no_bar,ele_lim=ele_lim,xrange=xrange,yrange=yrange

	COMMON prf_psa
	COMMON fov_psa
        COMMON tif_psa

; Plot area
	IF KEYWORD_SET(position) THEN !P.POSITION=position
	pos=!P.POSITION

; Set up the elevation angle limitation
	IF NOT KEYWORD_SET(ele_lim) THEN ele_lim=10.0

; Set up the plotting area
	IF NOT KEYWORD_SET(xrange) THEN BEGIN
		xrange=[4,34]
		IF cid_tif EQ 2 THEN xrange=[12,41] ; offset is 29
	ENDIF
	IF NOT KEYWORD_SET(yrange) THEN BEGIN
		yrange=[64,75]
		IF cid_tif EQ 2 THEN yrange=[62,73] ; offset is 11
	ENDIF

; Plot frame
        PLOT,[0],[0],XRANGE=xrange,YRANGE=yrange,/XSTYLE,/YSTYLE,XTITLE='Geo. Lon.',YTITLE='Geo. Lat'

; Plot the data
	IF KEYWORD_SET(smo) THEN BEGIN
		tmp_sum=FLTARR(256,256)
		FOR i=0,smo-1 DO tmp_sum=tmp_sum+REFORM(img_tif(now_tif+i,*,*))
		tmp_img=tmp_sum/smo
	ENDIF ELSE BEGIN
		tmp_img=REFORM(img_tif(now_tif,*,*))
	ENDELSE

; Prepare array for contouring
        tmp_bin=WHERE(gla_psa NE 99999.9 AND gla_psa NE 99999.9 AND ele_psa GE ele_lim,no_tmp_bin)
        x_cont=MAKE_ARRAY(no_tmp_bin,/FLOAT,VALUE=0.0)
        y_cont=MAKE_ARRAY(no_tmp_bin,/FLOAT,VALUE=0.0)
        z_cont=MAKE_ARRAY(no_tmp_bin,/FLOAT,VALUE=0.0)

; Put values into the arrays for contouring
        num_cnt=0L
        FOR x=0,255 DO BEGIN
                FOR y=0,255 DO BEGIN
                        IF gla_psa(x,y) NE 99999.9 AND glo_psa(x,y) NE 99999.9 AND ele_psa(x,y) GE ele_lim THEN BEGIN

                                z_cont(num_cnt)=FLOAT(tmp_img(x,y));*COS(zan_psa(x,y)*!DTOR)
                                IF z_cont(num_cnt) GT max_psa THEN z_cont(num_cnt)=max_psa
                                IF z_cont(num_cnt) LT min_psa THEN z_cont(num_cnt)=min_psa

                                x_cont(num_cnt)=glo_psa(x,y)
                                y_cont(num_cnt)=gla_psa(x,y)

                                num_cnt++
                        ENDIF
                ENDFOR
        ENDFOR

; Contouring
        contour_lines=30
        levels_set=(max_psa-min_psa)*FINDGEN(contour_lines)/FLOAT(contour_lines)+min_psa
        colour_set=!D.TABLE_SIZE*INDGEN(contour_lines)/contour_lines
        CONTOUR,z_cont,x_cont,y_cont,/OVERPLOT,/FILL,LEVELS=levels_set,C_COLORS=colour_set,/IRREGULAR

; Plot coastlines
	overlay_coast

; Plot messages
	pos=!P.POSITION
	XYOUTS,pos(0),pos(3)+0.01,day_tif(now_tif)+': Cam'+STRING(cid_tif,FORMAT='(I1.1)'),/NORMAL
	XYOUTS,pos(2),pos(3)+0.01,ti_tif(now_tif),/NORMAL,ALI=1

; Plot colour bar
        IF NOT KEYWORD_SET(no_bar) THEN BEGIN
		position=[!P.POSITION(2)+0.01,!P.POSITION(1),!P.POSITION(2)+0.02,!P.POSITION(3)]
		IF typ_tif EQ 1 THEN legend='Absolute Optical Intensity (R)' ELSE legend='Raw Count'
		plot_psa_colour_bar,position=position,legend=legend
	ENDIF

	END

;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	TIME_AXIS_INTERVAL
;

        PRO time_axis_intervals,start,stop,minor_ticks,where_tick,no_of_ticks

; Calculate length of data set in seconds
        length_in_time=LONG((stop-start)*3600+2)

; create array
        temp_time=LINDGEN(length_in_time) + long(start*3600)

; tick intervals depend on length of data set
        CASE 1 OF
                (length_in_time LE 402L) : BEGIN
                        time_interval = 60
                        minor_ticks = 2
                        END
                (length_in_time GT 402L) AND (length_in_time LE 902L) : BEGIN
                        time_interval = 120
                        minor_ticks = 4
                        END
                (length_in_time GT 902L) AND (length_in_time LE 1802L) : BEGIN
                        time_interval = 300
                        minor_ticks = 5
                        END
                (length_in_time GT 1802L) AND (length_in_time LE 3602L) : BEGIN
                        time_interval = 600
                        minor_ticks = 10
                        END
                (length_in_time GT 3602L) AND (length_in_time LE 7202L) : BEGIN
                        time_interval = 1200
                        minor_ticks=4
                        END
                (length_in_time GT 7202L) AND (length_in_time LE 14402L) : BEGIN
                        time_interval = 1800
                        minor_ticks=3
                        END
                (length_in_time GT 14402L) AND (length_in_time LE 28802L) : BEGIN
                        time_interval = 1*3600
                        minor_ticks=6
                        END
                (length_in_time GT 28802L) AND (length_in_time LE 57602L) : BEGIN
                        time_interval = 2*3600
                        minor_ticks=4
                        END
                (length_in_time GT 57602L) AND (length_in_time LE 86402L) : BEGIN
                        time_interval = 4*3600
                        minor_ticks=4
                        END
                (length_in_time GT 86402L) AND (length_in_time LE 115202L) : BEGIN
                        time_interval = 6*3600
                        minor_ticks=6
                        END
                (length_in_time GT 115202L) AND (length_in_time LE 172802L) : BEGIN
                        time_interval = 8*3600L
                        minor_ticks=4
                        END
                (length_in_time GT 172802L) : BEGIN
                        time_interval = 12*3600L
                        minor_ticks=6
                        END
                ELSE: BEGIN
                        time_interval = 7200L
                        minor_ticks=2
                END

        ENDCASE

; find where ticks should go - I'm sure this could be made more efficient
        no_of_ticks=0
        where_tick=DBLARR(30)
        FOR k=0L,length_in_time-1 DO BEGIN
                IF ((temp_time(k) MOD time_interval) EQ 0) THEN BEGIN
                        IF (no_of_ticks GT 29) THEN STOP, 'ERROR: too many ticks'
                        where_tick(no_of_ticks)=temp_time(k)/DOUBLE(3600)
                        no_of_ticks=no_of_ticks+1
                ENDIF
        ENDFOR

        no_of_ticks=no_of_ticks-1
        where_tick=where_tick(0:no_of_ticks)

        END

;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	DPANEL_PSA
;

        PRO dpanel_psa,xmaps,ymaps,xmap,ymap,square=square,landscape=landscape

	COMMON prf_psa

; Default is one panel on screen
        IF N_PARAMS() NE 4 THEN BEGIN & xmaps=1 & ymaps=1 & xmap=0 & ymap=0 & ENDIF

; Check for bad x and y
        IF xmap GE xmaps OR ymap GE ymaps THEN BEGIN
                PRINT,'WARNING: panel positions out of range'
                xmap=0
                ymap=0
        ENDIF
        
        last_xmaps=xmaps & last_ymaps=ymaps & last_xmap=xmap & last_ymap=ymap

; Initialize plotting preferences
; x,ysize       -       proportion of screen to put panels in
; x,yorigin     -       where to start page from
; l,r,t,bmargin -       left, right, top and bottom margins around plot
;                       window as fractions of the panel
; If set_format,/sardines is in force then tmargin=bmargin=0 and move
; things about slightly
	xsize=0.83
	ysize=0.87
        xorigin=0.02
        yorigin=0.03
        lmargin=0.15
        rmargin=0.05
        tmargin=0.05
        bmargin=0.15
        IF fmt_psa EQ 1 THEN BEGIN
                ysize=0.80
                yorigin=0.08
        ENDIF

; Calculate size of each panel
        xframe=xsize/xmaps
        yframe=ysize/ymaps

; If /SQUARE option is set then constrain plotting window to be square -
; recalculate xframe and yframe accordingly, taking into account the
; device aspect ratio
        IF KEYWORD_SET(square) THEN BEGIN
                aspect_ratio=1.0*!D.Y_SIZE/!D.X_SIZE
                xpanel=xframe*(1-lmargin-rmargin)
                ypanel=yframe*(1-tmargin-bmargin)
                IF xmaps*ypanel*aspect_ratio/(1-lmargin-rmargin) GT xsize THEN $
                        ypanel=xpanel/aspect_ratio ELSE xpanel=ypanel*aspect_ratio      
                xframe=xpanel/(1-lmargin-rmargin)
                yframe=ypanel/(1-tmargin-bmargin)
        ENDIF

        x1=(xmap+lmargin)*xframe
        y1=(ymaps-ymap-1+bmargin)*yframe
        x2=(xmap+1-rmargin)*xframe
        y2=(ymaps-ymap-tmargin)*yframe

; If panels are forced square, then centre plotting area
        xcentre=(xsize-xframe*xmaps)*0.5
        ycentre=(ysize-yframe*ymaps)*0.5

        !P.POSITION=[x1+xcentre+xorigin,y1+ycentre+yorigin,x2+xcentre+xorigin,y2+ycentre+yorigin]

        END

;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	WINDOW_PSA
;
;	PsA > window_psa
;	PsA > window_psa,/landscape
;

        PRO window_psa,size,landscape=landscape

	COMMON prf_psa

        IF N_PARAMS() EQ 0 THEN size=900
        xsize=size
        ysize=size

        IF KEYWORD_SET(landscape) THEN xsize=xsize*SQRT(2) ELSE ysize=ysize*SQRT(2)
	IF KEYWORD_SET(landscape) THEN fmt_psa=1 ELSE fmt_psa=0

        IF !D.X_SIZE NE FIX(xsize) OR !D.Y_SIZE NE FIX(ysize) OR !D.WINDOW EQ -1 THEN $
		WINDOW,XSIZE=xsize,YSIZE=ysize,TITLE='Kiban-S output:',RETAIN=2
        ERASE

        END

;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       PLOT_PSA_COLOUR_BAR
;
; PURPOSE:
;
;       Plot colour bar for PsA plotting routines.
;

        PRO plot_psa_colour_bar,ymaps,ymap,position=position,leg_pos=leg_pos, $
                legend=legend,no_ticks=no_ticks,force_char=force_char,no_colours=no_colours

        COMMON prf_psa

; Handling the keywords
        IF NOT KEYWORD_SET(no_colours) THEN no_colours=10
        IF NOT KEYWORD_SET(no_ticks) THEN no_ticks=no_colours

; Allow several colour bars to be stacked
        IF N_PARAMS() NE 2 THEN BEGIN & ymaps=1 & ymap=0 & ENDIF

; Initialize colour bar position
        ysize=0.83
        xorigin=0.05
        yorigin=0.0
        IF fmt_psa EQ 1 THEN BEGIN
                ysize=0.73
                yorigin=0.1
        ENDIF
        xpos=0.85+xorigin
        xbox=MAX([0.02/ymaps,0.01])
        ypos=(ymaps-ymap-0.75)*ysize/ymaps+yorigin
        ybox_cols =0.6*ysize/(ymaps*no_colours)
        ybox_ticks=0.6*ysize/(ymaps*no_ticks)
        ybox_gnd  =0.6*ysize/(ymaps*10)
        IF KEYWORD_SET(force_char) THEN char=force_char ELSE char=!P.CHARSIZE
        
        IF KEYWORD_SET(position) THEN BEGIN
                xpos= position(0)
                ypos= position(1)
                xbox= position(2)-position(0)
                ybox_cols =(position(3)-position(1))/no_colours
                ybox_ticks=(position(3)-position(1))/no_ticks
                char=!P.CHARSIZE
        ENDIF

; Setting colours
	ncol=252
        cin=FIX(FINDGEN(no_colours)*(ncol-4)/(no_colours-1))+1
        lvl=min_psa+FINDGEN(no_colours)*(max_psa-min_psa)/no_colours

; Draw coloured boxes
        FOR level=0,no_colours-1 DO                                     $
                POLYFILL,[xpos,xpos+xbox,xpos+xbox,xpos],               $
                         [ypos+ybox_cols*level,ypos+ybox_cols*level,            $
                          ypos+ybox_cols*(level+1),ypos+ybox_cols*(level+1)],   $
                          COLOR=cin(level),/NORMAL

; Draw outline
        FOR level=0,no_ticks-1 DO                                       $
                PLOTS,xpos+xbox*[0,1,1,0,0],ypos+ybox_ticks*(level+[0,0,1,1,0]), $
                        COLOR=foreground,/NORMAL

; Plot levels
        IF FIX((max_psa-min_psa)/no_colours) NE FLOAT((max_psa-min_psa))/no_colours THEN BEGIN
                level_format='(F10.1)'
        ENDIF ELSE BEGIN
                level_format='(I)'
        ENDELSE
        lvl=min_psa+FINDGEN(no_ticks)*(max_psa-min_psa)/no_ticks
        FOR level=0,no_ticks-1 DO BEGIN
                numb=STRTRIM(FIX(ABS(lvl(level))),2)+'.'+STRTRIM(ABS(FIX(lvl(level)*10)) MOD 10,2)
                IF lvl(level) LT 0 THEN numb='-'+numb
                numb=STRTRIM(STRING(lvl(level),FORMAT=level_format),2)
                XYOUTS,xpos+1.4*xbox,ypos+ybox_ticks*level,                $
                        numb,COLOR=foreground,CHARSIZE=0.8*char,/NORMAL,ORI=ori_char,ALI=ali
        ENDFOR

; Plot title
        title='Raw Count'
        IF KEYWORD_SET(legend) THEN title=legend
        IF NOT KEYWORD_SET(leg_pos) THEN leg_pos=1.0
        XYOUTS,xpos+0.075*leg_pos,ypos+no_colours*ybox_cols*0.5,title,COLOR=foreground, $
                ORIENTATION=270,ALIGNMENT=0.5,CHARSIZE=char,/NORMAL

        END
;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       PLOT_COLOUR_BAR
;
; PURPOSE:
;
;       Plot colour bar for plotting routines.
;

        PRO plot_colour_bar,ymaps,ymap,position=position,leg_pos=leg_pos, $
                legend=legend,no_ticks=no_ticks,force_char=force_char,no_colours=no_colours,bar_title=bar_title,max_val=max_val,min_val=min_val

          COMMON prf_psa
          COMMON cor_psa

; Handling the keywords
        IF NOT KEYWORD_SET(no_colours) THEN no_colours=10
        IF NOT KEYWORD_SET(no_ticks) THEN no_ticks=no_colours

; Allow several colour bars to be stacked
        IF N_PARAMS() NE 2 THEN BEGIN & ymaps=1 & ymap=0 & ENDIF

; Initialize colour bar position
        ysize=0.82
        xorigin=0.05
        yorigin=0.0
        ;max_cor=1.0
        ;min_cor=0.0
        IF fmt_psa EQ 1 THEN BEGIN
                ysize=0.72
                yorigin=0.1
        ENDIF
        xpos=0.85+xorigin
        xbox=MAX([0.02/ymaps,0.01])
        ypos=(ymaps-ymap-0.75)*ysize/ymaps+yorigin
        ybox_cols =0.6*ysize/(ymaps*no_colours)
        ybox_ticks=0.6*ysize/(ymaps*no_ticks)
        ybox_gnd  =0.6*ysize/(ymaps*10)
        IF KEYWORD_SET(force_char) THEN char=force_char ELSE char=!P.CHARSIZE
        
        IF KEYWORD_SET(position) THEN BEGIN
                xpos= position(0)
                ypos= position(1)
                xbox= position(2)-position(0)
                ybox_cols =(position(3)-position(1))/no_colours
                ybox_ticks=(position(3)-position(1))/no_ticks
                char=!P.CHARSIZE
        ENDIF

; Setting colours
	ncol=252
        cin=FIX(FINDGEN(no_colours)*(ncol-4)/(no_colours-1))+1
        lvl=min_psa+FINDGEN(no_colours)*(max_val-min_val)/no_colours

; Draw coloured boxes
        FOR level=0,no_colours-1 DO                                     $
                POLYFILL,[xpos,xpos+xbox,xpos+xbox,xpos],               $
                         [ypos+ybox_cols*level,ypos+ybox_cols*level,            $
                          ypos+ybox_cols*(level+1),ypos+ybox_cols*(level+1)],   $
                          COLOR=cin(level),/NORMAL

; Draw outline
        FOR level=0,no_ticks-1 DO                                       $
                PLOTS,xpos+xbox*[0,1,1,0,0],ypos+ybox_ticks*(level+[0,0,1,1,0]), $
                        COLOR=foreground,/NORMAL

; Plot levels
        IF FIX((max_val-min_val)/no_colours) NE FLOAT((max_val-min_val))/no_colours THEN BEGIN
                level_format='(F10.1)'
        ENDIF ELSE BEGIN
                level_format='(I)'
        ENDELSE
        lvl=min_val+FINDGEN(no_ticks)*(max_val-min_val)/no_ticks
        FOR level=0,no_ticks-1 DO BEGIN
                numb=STRTRIM(FIX(ABS(lvl(level))),2)+'.'+STRTRIM(ABS(FIX(lvl(level)*10)) MOD 10,2)
                IF lvl(level) LT 0 THEN numb='-'+numb
                numb=STRTRIM(STRING(lvl(level),FORMAT=level_format),2)
                XYOUTS,xpos+1.4*xbox,ypos+ybox_ticks*level,                $
                        numb,COLOR=foreground,CHARSIZE=0.8*char,/NORMAL,ORI=ori_char,ALI=ali
        ENDFOR

; Plot title
        IF bar_title EQ 'cor' THEN title='Correlate Coefficient'
        IF bar_title EQ 'fft' THEN title='Peak Period'
        IF KEYWORD_SET(legend) THEN title=legend
        IF NOT KEYWORD_SET(leg_pos) THEN leg_pos=1.0
        XYOUTS,xpos+0.075*leg_pos,ypos+no_colours*ybox_cols*0.5,title,COLOR=foreground, $
                ORIENTATION=270,ALIGNMENT=0.5,CHARSIZE=char,/NORMAL

        END

;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       OVERLAY_COAST
;

        PRO overlay_coast,col=col

	COMMON prf_psa

        OPENR,map_unit,dir_psa+'fovd/world_data',/GET_LUN

        no_blocks=17
        block_len=INTARR(17)
        READF,map_unit,no_blocks,block_len

        coast=FLTARR(2,10000) & pts=0
        FOR read_block=0,no_blocks-1 DO BEGIN
                read_coast=FLTARR(2,block_len(read_block))
                READF,map_unit,read_coast
                coast(*,pts:pts+block_len(read_block)-1)=read_coast(*,*)
                pts=pts+block_len(read_block)
        ENDFOR
        CLOSE,map_unit
        FREE_LUN,map_unit
        
        plot_coast=FLTARR(2,5000) & plot_pts=0
        FOR i=0,pts-1 DO BEGIN
                IF (coast(0,i) NE 0 OR coast(1,i) NE 0) THEN BEGIN
                        IF coast(0,i) GT 0 THEN BEGIN
                                IF coast(1,i) LT 0.0 THEN coast(1,i)+=360.0
                                plot_coast(*,plot_pts)=[coast(1,i),coast(0,i)]
                                plot_pts=plot_pts+1
                        ENDIF
                ENDIF ELSE BEGIN
                        IF plot_pts GT 0 THEN BEGIN
                                OPLOT,plot_coast(0,[INDGEN(plot_pts),0]),plot_coast(1,[INDGEN(plot_pts),0]),COL=col,NOCLIP=0
                        ENDIF
                        plot_pts=0
                ENDELSE
        ENDFOR

        END

;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	EBIREADED_YM - Just read one image block in the data file
;

	PRO ebireaded_ym,id,itag,x,y,imgname,dat

  	iTag=0L
  	READU,id,iTag
  	CASE iTag of
    	1000: BEGIN ; header
      	iSize=0L
      	READU,id,iSize
      	tmp=bytarr(iSize)
      	READU,id,tmp
    	END
    	1001: BEGIN ; header
      	iSize=0L
      	READU,id,iSize
      	iDepth=0L
      	READU,id,iDepth
      	x=0L
      	READU,id,x
      	y=0L
      	READU,id,y
    	END
    	1002: BEGIN ; filename
      	iSize=0L
      	READU,id,iSize
      	cname=bytarr(iSize)
      	READU,id,cname
      	imgname=''
      	for i=0L, iSize-1 do BEGIN
        	if (cname[i] eq 0) then break
        	imgname=imgname+string(cname[i])
      	ENDfor
    	END
    	2000: BEGIN
      	iSize=0L
      	READU,id,iSize
      	dat=uintarr(x,y)
      	READU,id,dat ; read 16 bit data
    	END
    	9999: BEGIN
    	END
  	ENDCASE
	
	END

;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	PSYM8
;
; PURPOSE:
;
;	Setting original symbol
;
; EXAMPLE:
;
;	PsA > psym8,/FILL
;	PsA > OPLOT,[10],[10],PSYM=8
;

        PRO psym8,fill=fill,dot=dot

        cir=findgen(31)*(!pi*2/30)
        IF KEYWORD_SET(fill) THEN $
                USERSYM,cos(cir),sin(cir),/fill ELSE $
                USERSYM,cos(cir),sin(cir)
        END

;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	PLOT_PSA_TIF_MULTI_FOR_MOVIE
;

	PRO plot_psa_tif_multi_for_movie,date,shour,ehour,fin=fin,ask=ask,erg=erg

; Charsize setting
	!P.CHARSIZE=1.5

; Set device to Z-buffer and make imaginary window
       	SET_PLOT,'Z'
       	DEVICE,SET_RESOLUTION=[800,800],DECOMPOSED=0,SET_PIXEL_DEPTH=24
	dpanel_psa,1,1,0,0,/square
	set_scale_psa,2000,4000

; Make directory for storing png files
	out_dir=date+'_'+STRING(shour,FORMAT='(I2.2)')+'-'+STRING(ehour,FORMAT='(I2.2)')
	SPAWN,'mkdir -p '+out_dir

; Loop for images
	n=0
	FOR hh=shour,ehour-1 DO BEGIN
		FOR mm=0,59 DO BEGIN
			FOR ss=0,50,10 DO BEGIN

				hms=hh*10000L+mm*100L+ss
				ERASE
				plot_psa_tif_multi_panel,date,hms,ele_lim=20.0,fin=fin,ask=ask,erg=erg
				WRITE_PNG,'./'+out_dir+'/'+STRING(n,FORMAT='(I5.5)')+'.png',TVRD(true=1)
				n++

			ENDFOR
		ENDFOR
	ENDFOR

	END

;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	PLOT_PSA_TIF_MULTI_PANEL
;

	PRO plot_psa_tif_multi_panel,date,hhmmss,abs=abs,alt=alt,position=position,ele_lim=ele_lim,erg=erg,fin=fin,ask=ask

	COMMON prf_psa

        IF NOT KEYWORD_SET(alt) THEN alt=110

; If Finland
	IF KEYWORD_SET(fin) THEN BEGIN
		plot_psa_tif_multi_panel_fin,date,hhmmss,abs=abs,alt=alt,position=position,ele_lim=ele_lim,erg=erg
		RETURN
	ENDIF

; If Alaska
	IF KEYWORD_SET(ask) THEN BEGIN
		plot_psa_tif_multi_panel_ask,date,hhmmss,abs=abs,alt=alt,position=position,ele_lim=ele_lim,erg=erg
		RETURN
	ENDIF

; The number of cameras and pixel size
	num_cam=3
	pxl=256
	
; Break date
	yr_read=STRMID(date,0,4)
	mo_read=STRMID(date,4,2)
	dy_read=STRMID(date,6,2)
	hms_sel=STRMID(STRING(hhmmss,FORMAT='(I6.6)'),0,5)
	hh_read=STRMID(STRING(hhmmss,FORMAT='(I6.6)'),0,2)
	mm_read=STRMID(STRING(hhmmss,FORMAT='(I6.6)'),2,2)
	erg_tim=yr_read+'-'+mo_read+'-'+dy_read+'/'+hh_read+':'+mm_read+':00'

; Print message
	PRINT,date+': '+STRING(hhmmss,FORMAT='(I6.6)')+' UT'

; Make arrays for images
	avl_map=INTARR(num_cam)
	tif_map=FLTARR(num_cam,pxl,pxl)
        azm_map=FLTARR(num_cam,pxl,pxl)
        ele_map=FLTARR(num_cam,pxl,pxl)
        gla_map=FLTARR(num_cam,pxl,pxl)
        glo_map=FLTARR(num_cam,pxl,pxl)

; Path for the data
	;camids=[1,2,6]
	camids=[1,6,2]
	FOR n=0,num_cam-1 DO BEGIN
		camid=camids(n)

; Read the data
		path=dir_psa+'cam'+STRING(camid,FORMAT='(I1.1)')+'/'+yr_read+'/'+mo_read+'/'+dy_read+'_QL'
        	fname=path+'/C'+STRING(camid,FORMAT='(I1.1)')+'_'+date+'_'+hms_sel+'*.tif'
		found_files=FILE_SEARCH(fname,COUNT=num_found_files)
		IF num_found_files EQ 1 THEN BEGIN

		avl_map(n)=1
               	img=READ_TIFF(found_files(0))
		tif_map(n,*,*)=ROTATE(img,7)

; Read azimuth and elevation
	        found_file=FILE_SEARCH(dir_psa+'fovd/azm_ele/C'+STRING(camid,FORMAT='(I1.1)')+'*az.txt',COUNT=num_found_file)
        	OPENR,azm_lun,found_file(num_found_file-1),/GET_LUN
		;PRINT,'Azimuth FOV file:',found_file(num_found_file-1)
	        found_file=FILE_SEARCH(dir_psa+'fovd/azm_ele/C'+STRING(camid,FORMAT='(I1.1)')+'*el.txt',COUNT=num_found_file)
        	OPENR,ele_lun,found_file(num_found_file-1),/GET_LUN
		;PRINT,'Elevation FOV file:',found_file(num_found_file-1)
        	FOR i=0,pxl-1 DO BEGIN
                	tmp_read=FLTARR(pxl)
                	READF,ele_lun,tmp_read
                	ele_map(n,*,pxl-1-i)=tmp_read
                	READF,azm_lun,tmp_read
                	azm_map(n,*,pxl-1-i)=tmp_read
        	ENDFOR
        	FREE_LUN,azm_lun
        	FREE_LUN,ele_lun

; Read glat and glon info
        	found_file=FILE_SEARCH(dir_psa+'fovd/gla_glo/C'+STRING(camid,FORMAT='(I1.1)')+'*'+ $
                	STRING(alt,FORMAT='(I3.3)')+'_lat.txt',COUNT=num_found_file)
        	OPENR,gla_lun,found_file(num_found_file-1),/GET_LUN
		;PRINT,'Latitude FOV file:',found_file(num_found_file-1)
        	found_file=FILE_SEARCH(dir_psa+'fovd/gla_glo/C'+STRING(camid,FORMAT='(I1.1)')+'*'+ $
                	STRING(alt,FORMAT='(I3.3)')+'_long.txt',COUNT=num_found_file)
        	OPENR,glo_lun,found_file(num_found_file-1),/GET_LUN
		;PRINT,'Longitude FOV file:',found_file(num_found_file-1)
        	FOR i=0,pxl-1 DO BEGIN
                	tmp_read=FLTARR(pxl)
                	READF,gla_lun,tmp_read
			IF camid EQ 6 THEN gla_map(n,*,i)=tmp_read ELSE gla_map(n,*,pxl-1-i)=tmp_read
                	READF,glo_lun,tmp_read
			IF camid EQ 6 THEN glo_map(n,*,i)=tmp_read ELSE glo_map(n,*,pxl-1-i)=tmp_read
        	ENDFOR
        	FREE_LUN,gla_lun
        	FREE_LUN,glo_lun

		ENDIF

        ENDFOR

; Plot area
	IF KEYWORD_SET(position) THEN !P.POSITION=position
	pos=!P.POSITION

; Set up the elevation angle limitation
	IF NOT KEYWORD_SET(ele_lim) THEN ele_lim=15.0

; Set up the plotting area
	IF NOT KEYWORD_SET(xrange) THEN BEGIN
		xrange=[5,40]
	ENDIF
	IF NOT KEYWORD_SET(yrange) THEN BEGIN
		yrange=[63,74]
	ENDIF

; Plot frame
        PLOT,[0],[0],XRANGE=xrange,YRANGE=yrange,/XSTYLE,/YSTYLE,XTITLE='Geo. Lon.',YTITLE='Geo. Lat'

	FOR n=0,num_cam-1 DO BEGIN

		IF avl_map(n) EQ 1 THEN BEGIN

		camid=camids(n)

; Prepare array for contouring
		gla_psa=REFORM(gla_map(n,*,*))
		glo_psa=REFORM(glo_map(n,*,*))
		azm_psa=REFORM(azm_map(n,*,*))
		ele_psa=REFORM(ele_map(n,*,*))
		tmp_img=REFORM(tif_map(n,*,*))
        	tmp_bin=WHERE(gla_psa NE 99999.9 AND gla_psa NE 99999.9 AND ele_psa GE ele_lim,no_tmp_bin)
        	x_cont=MAKE_ARRAY(no_tmp_bin,/FLOAT,VALUE=0.0)
        	y_cont=MAKE_ARRAY(no_tmp_bin,/FLOAT,VALUE=0.0)
        	z_cont=MAKE_ARRAY(no_tmp_bin,/FLOAT,VALUE=0.0)

; Put values into the arrays for contouring
        	num_cnt=0L

		IF camid EQ 1 THEN BEGIN
        	FOR x=0,255 DO BEGIN
                	FOR y=0,255 DO BEGIN
                        	IF gla_psa(x,y) NE 99999.9 AND glo_psa(x,y) NE 99999.9 AND ele_psa(x,y) GE ele_lim $
					AND glo_psa(x,y) LE 22.6 THEN BEGIN

                                	z_cont(num_cnt)=FLOAT(tmp_img(x,y));*COS(zan_psa(x,y)*!DTOR)
                                	IF z_cont(num_cnt) GT max_psa THEN z_cont(num_cnt)=max_psa
                                	IF z_cont(num_cnt) LT min_psa THEN z_cont(num_cnt)=min_psa

                                	x_cont(num_cnt)=glo_psa(x,y)
                                	y_cont(num_cnt)=gla_psa(x,y)

                                	num_cnt++
                        	ENDIF
                	ENDFOR
        	ENDFOR
		x_cont=x_cont(0:num_cnt-1)
		y_cont=y_cont(0:num_cnt-1)
		z_cont=z_cont(0:num_cnt-1)
		ENDIF

		IF camid EQ 2 THEN BEGIN
        	FOR x=0,255 DO BEGIN
                	FOR y=0,255 DO BEGIN
                        	IF gla_psa(x,y) NE 99999.9 AND glo_psa(x,y) NE 99999.9 AND ele_psa(x,y) GE ele_lim $
					AND 0.748*(glo_psa(x,y)-19.96)+67.06 GE gla_psa(x,y) $
					AND -0.0757*(glo_psa(x,y)-32.78)+68.28 GE gla_psa(x,y) THEN BEGIN

                                	z_cont(num_cnt)=FLOAT(tmp_img(x,y));*COS(zan_psa(x,y)*!DTOR)
                                	IF z_cont(num_cnt) GT max_psa THEN z_cont(num_cnt)=max_psa
                                	IF z_cont(num_cnt) LT min_psa THEN z_cont(num_cnt)=min_psa

                                	x_cont(num_cnt)=glo_psa(x,y)
                                	y_cont(num_cnt)=gla_psa(x,y)

                                	num_cnt++
                        	ENDIF
                	ENDFOR
        	ENDFOR
		x_cont=x_cont(0:num_cnt-1)
		y_cont=y_cont(0:num_cnt-1)
		z_cont=z_cont(0:num_cnt-1)
		ENDIF

		IF camid EQ 6 THEN BEGIN
        	FOR x=0,255 DO BEGIN
                	FOR y=0,255 DO BEGIN
                        	IF gla_psa(x,y) NE 99999.9 AND glo_psa(x,y) NE 99999.9 AND ele_psa(x,y) GE ele_lim $
					AND glo_psa(x,y) GT 22.55 THEN BEGIN

                                	z_cont(num_cnt)=FLOAT(tmp_img(x,y));*COS(zan_psa(x,y)*!DTOR)
                                	IF z_cont(num_cnt) GT max_psa THEN z_cont(num_cnt)=max_psa
                                	IF z_cont(num_cnt) LT min_psa THEN z_cont(num_cnt)=min_psa

                                	x_cont(num_cnt)=glo_psa(x,y)
                                	y_cont(num_cnt)=gla_psa(x,y)

                                	num_cnt++
                        	ENDIF
                	ENDFOR
        	ENDFOR
		x_cont=x_cont(0:num_cnt-1)
		y_cont=y_cont(0:num_cnt-1)
		z_cont=z_cont(0:num_cnt-1)
		ENDIF

; Contouring
        	contour_lines=30
        	levels_set=(max_psa-min_psa)*FINDGEN(contour_lines)/FLOAT(contour_lines)+min_psa
        	colour_set=!D.TABLE_SIZE*INDGEN(contour_lines)/contour_lines
        	CONTOUR,z_cont,x_cont,y_cont,/OVERPLOT,/FILL,LEVELS=levels_set,C_COLORS=colour_set,/IRREGULAR

		ENDIF
	ENDFOR

        IF KEYWORD_SET(erg) THEN BEGIN

                fname='/radar01/work/Aurora/Kiban-S/Ground/Campaign/orb_txt/'+date+'.txt'
                SPAWN,'grep '+erg_tim+' '+fname+' > tmp_erg_orb.txt'

                la_fot=0.0
                ma_fot=0.0
                lo_fot=0.0
                ml_fot=0.0
                hm_fot='abc'

                OPENR,inp_lun,'tmp_erg_orb.txt',/GET_LUN

                tmp_string='abc'
                READF,inp_lun,tmp_string
                hh=FIX(STRMID(tmp_string,11,2))
                mm=FIX(STRMID(tmp_string,14,2))
                ss=FIX(STRMID(tmp_string,17,2))
                hm_fot=STRING(hh,FORMAT='(I2.2)')+STRING(mm,FORMAT='(I2.2)')+' UT'
                tmp_read=FLTARR(4)
                READS,STRMID(tmp_string,21,173),tmp_read
                la_fot=tmp_read(0)
                lo_fot=tmp_read(1)
                LOADCT,39,/SILENT
                psym8,/FILL
                OPLOT,[lo_fot],[la_fot],PSYM=8,COL=200,SYMSIZE=1.0
                psym8
                OPLOT,[lo_fot],[la_fot],PSYM=8,COL=255,SYMSIZE=1.0
                XYOUTS,lo_fot,la_fot+0.2,'ARASE',/DATA,NOCLIP=0,COL=200,ALI=0.5,CHARSIZE=!P.CHARSIZE
                LOADCT,0,/SILENT
                FREE_LUN,inp_lun
        ENDIF

; Plot coastlines
	overlay_coast

; Overlay stations
        gla_sta=[69.58,69.76,67.37,67.28,67.84,69.07]
        glo_sta=[19.23,27.01,26.63,20.72,20.42,20.76]
        let_sta=['TRO','KEV','SOD','TJA','KRN','KIL']
        PSYM8,/FILL
        LOADCT,39,/SILENT
        FOR i=0,N_ELEMENTS(let_sta)-1 DO BEGIN
                OPLOT,[glo_sta(i)],[gla_sta(i)],PSYM=8,SYMSIZE=0.5,COL=250
                OPLOT,[glo_sta(i)],[gla_sta(i)],PSYM=1,SYMSIZE=1.0,COL=250
                XYOUTS,glo_sta(i),gla_sta(i)+0.2,let_sta(i),COL=252,ALI=0.5,CHARSIZE=0.6*!P.CHARSIZE
        ENDFOR
        LOADCT,0,/SILENT

; Plot messages
	pos=!P.POSITION
	XYOUTS,pos(0),pos(3)+0.01,date,/NORMAL
	XYOUTS,pos(2),pos(3)+0.01,STRING(hhmmss,FORMAT='(I6.6)')+' UT',/NORMAL,ALI=1

; Plot colour bar
	typ_raw=0
        IF NOT KEYWORD_SET(no_bar) THEN BEGIN
		position=[!P.POSITION(2)+0.01,!P.POSITION(1),!P.POSITION(2)+0.02,!P.POSITION(3)]
		IF typ_raw EQ 1 THEN legend='Absolute Optical Intensity (R)' ELSE legend='Raw Count'
		plot_psa_colour_bar,position=position,legend=legend
	ENDIF

	skip0:

	END

;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	PLOT_PSA_TIF_MULTI_PANEL_FIN
;

	PRO plot_psa_tif_multi_panel_fin,date,hhmmss,abs=abs,alt=alt,position=position,ele_lim=ele_lim,check=check,erg=erg

	COMMON prf_psa

        IF NOT KEYWORD_SET(alt) THEN alt=110

	IF KEYWORD_SET(check) THEN BEGIN

       		SET_PLOT,'X'
		WINDOW,XSIZE=800,YSIZE=800
		dpanel_psa,1,1,0,0,/square
		set_scale_psa,2000,3000
	ENDIF

; The number of cameras and pixel size
	num_cam=2
	pxl=256
	
; Break date
	yr_read=STRMID(date,0,4)
	mo_read=STRMID(date,4,2)
	dy_read=STRMID(date,6,2)
	hms_sel=STRMID(STRING(hhmmss,FORMAT='(I6.6)'),0,5)
	hh_read=STRMID(STRING(hhmmss,FORMAT='(I6.6)'),0,2)
	mm_read=STRMID(STRING(hhmmss,FORMAT='(I6.6)'),2,2)
	erg_tim=yr_read+'-'+mo_read+'-'+dy_read+'/'+hh_read+':'+mm_read+':00'

; Print message
	PRINT,date+': '+STRING(hhmmss,FORMAT='(I6.6)')+' UT'

; Make arrays for images
	avl_map=INTARR(num_cam)
	tif_map=FLTARR(num_cam,pxl,pxl)
        azm_map=FLTARR(num_cam,pxl,pxl)
        ele_map=FLTARR(num_cam,pxl,pxl)
        gla_map=FLTARR(num_cam,pxl,pxl)
        glo_map=FLTARR(num_cam,pxl,pxl)

; Path for the data
	camids=[6,2]
	FOR n=0,num_cam-1 DO BEGIN
		camid=camids(n)

; Read the data
		path=dir_psa+'cam'+STRING(camid,FORMAT='(I1.1)')+'/'+yr_read+'/'+mo_read+'/'+dy_read+'_QL'
        	fname=path+'/C'+STRING(camid,FORMAT='(I1.1)')+'_'+date+'_'+hms_sel+'*.tif'
		found_files=FILE_SEARCH(fname,COUNT=num_found_files)
		IF num_found_files EQ 1 THEN BEGIN

		avl_map(n)=1
               	img=READ_TIFF(found_files(0))
		tif_map(n,*,*)=ROTATE(img,7)

; Read azimuth and elevation
	        found_file=FILE_SEARCH(dir_psa+'fovd/azm_ele/C'+STRING(camid,FORMAT='(I1.1)')+'*az.txt',COUNT=num_found_file)
        	OPENR,azm_lun,found_file(num_found_file-1),/GET_LUN
	        found_file=FILE_SEARCH(dir_psa+'fovd/azm_ele/C'+STRING(camid,FORMAT='(I1.1)')+'*el.txt',COUNT=num_found_file)
        	OPENR,ele_lun,found_file(num_found_file-1),/GET_LUN
        	FOR i=0,pxl-1 DO BEGIN
                	tmp_read=FLTARR(pxl)
                	READF,ele_lun,tmp_read
                	ele_map(n,*,pxl-1-i)=tmp_read
                	READF,azm_lun,tmp_read
                	azm_map(n,*,pxl-1-i)=tmp_read
        	ENDFOR
        	FREE_LUN,azm_lun
        	FREE_LUN,ele_lun

; Read glat and glon info
        	found_file=FILE_SEARCH(dir_psa+'fovd/gla_glo/C'+STRING(camid,FORMAT='(I1.1)')+'*'+ $
                	STRING(alt,FORMAT='(I3.3)')+'_lat.txt',COUNT=num_found_file)
        	OPENR,gla_lun,found_file(num_found_file-1),/GET_LUN
        	found_file=FILE_SEARCH(dir_psa+'fovd/gla_glo/C'+STRING(camid,FORMAT='(I1.1)')+'*'+ $
                	STRING(alt,FORMAT='(I3.3)')+'_long.txt',COUNT=num_found_file)
        	OPENR,glo_lun,found_file(num_found_file-1),/GET_LUN
        	FOR i=0,pxl-1 DO BEGIN
                	tmp_read=FLTARR(pxl)
                	READF,gla_lun,tmp_read
			IF camid EQ 6 THEN gla_map(n,*,i)=tmp_read ELSE gla_map(n,*,pxl-1-i)=tmp_read
                	READF,glo_lun,tmp_read
			IF camid EQ 6 THEN glo_map(n,*,i)=tmp_read ELSE glo_map(n,*,pxl-1-i)=tmp_read
        	ENDFOR
        	FREE_LUN,gla_lun
        	FREE_LUN,glo_lun

		ENDIF

        ENDFOR

; Plot area
	IF KEYWORD_SET(position) THEN !P.POSITION=position
	pos=!P.POSITION

; Set up the elevation angle limitation
	IF NOT KEYWORD_SET(ele_lim) THEN ele_lim=15.0

; Set up the plotting area
	IF NOT KEYWORD_SET(xrange) THEN BEGIN
		xrange=[10,40]
		;IF cid_raw EQ 2 THEN xrange=[14,38]
	ENDIF
	IF NOT KEYWORD_SET(yrange) THEN BEGIN
		yrange=[63,74]
		;IF cid_raw EQ 2 THEN yrange=[63,72]
	ENDIF

; Plot frame
        PLOT,[0],[0],XRANGE=xrange,YRANGE=yrange,/XSTYLE,/YSTYLE,XTITLE='Geo. Lon.',YTITLE='Geo. Lat'

	FOR n=0,num_cam-1 DO BEGIN

		IF avl_map(n) EQ 1 THEN BEGIN

		camid=camids(n)

; Prepare array for contouring
		gla_psa=REFORM(gla_map(n,*,*))
		glo_psa=REFORM(glo_map(n,*,*))
		azm_psa=REFORM(azm_map(n,*,*))
		ele_psa=REFORM(ele_map(n,*,*))
		tmp_img=REFORM(tif_map(n,*,*))
        	tmp_bin=WHERE(gla_psa NE 99999.9 AND gla_psa NE 99999.9 AND ele_psa GE ele_lim,no_tmp_bin)
        	x_cont=MAKE_ARRAY(no_tmp_bin,/FLOAT,VALUE=0.0)
        	y_cont=MAKE_ARRAY(no_tmp_bin,/FLOAT,VALUE=0.0)
        	z_cont=MAKE_ARRAY(no_tmp_bin,/FLOAT,VALUE=0.0)

; Put values into the arrays for contouring
        	num_cnt=0L

		IF camid EQ 2 THEN BEGIN
        	FOR x=0,255 DO BEGIN
                	FOR y=0,255 DO BEGIN
                        	IF gla_psa(x,y) NE 99999.9 AND glo_psa(x,y) NE 99999.9 AND ele_psa(x,y) GE ele_lim $
					AND -0.0350*(glo_psa(x,y)-32.72)+68.26 GE gla_psa(x,y) $
					THEN BEGIN

                                	z_cont(num_cnt)=FLOAT(tmp_img(x,y));*COS(zan_psa(x,y)*!DTOR)
                                	IF z_cont(num_cnt) GT max_psa THEN z_cont(num_cnt)=max_psa
                                	IF z_cont(num_cnt) LT min_psa THEN z_cont(num_cnt)=min_psa

                                	x_cont(num_cnt)=glo_psa(x,y)
                                	y_cont(num_cnt)=gla_psa(x,y)

                                	num_cnt++
                        	ENDIF
                	ENDFOR
        	ENDFOR
		x_cont=x_cont(0:num_cnt-1)
		y_cont=y_cont(0:num_cnt-1)
		z_cont=z_cont(0:num_cnt-1)
		ENDIF

		IF camid EQ 6 THEN BEGIN
        	FOR x=0,255 DO BEGIN
                	FOR y=0,255 DO BEGIN
                        	IF gla_psa(x,y) NE 99999.9 AND glo_psa(x,y) NE 99999.9 AND ele_psa(x,y) GE ele_lim $
					THEN BEGIN
					;AND glo_psa(x,y) GT 22.55 THEN BEGIN

                                	z_cont(num_cnt)=FLOAT(tmp_img(x,y));*COS(zan_psa(x,y)*!DTOR)
                                	IF z_cont(num_cnt) GT max_psa THEN z_cont(num_cnt)=max_psa
                                	IF z_cont(num_cnt) LT min_psa THEN z_cont(num_cnt)=min_psa

                                	x_cont(num_cnt)=glo_psa(x,y)
                                	y_cont(num_cnt)=gla_psa(x,y)

                                	num_cnt++
                        	ENDIF
                	ENDFOR
        	ENDFOR
		x_cont=x_cont(0:num_cnt-1)
		y_cont=y_cont(0:num_cnt-1)
		z_cont=z_cont(0:num_cnt-1)
		ENDIF

; Contouring
        	contour_lines=30
        	levels_set=(max_psa-min_psa)*FINDGEN(contour_lines)/FLOAT(contour_lines)+min_psa
        	colour_set=!D.TABLE_SIZE*INDGEN(contour_lines)/contour_lines
        	CONTOUR,z_cont,x_cont,y_cont,/OVERPLOT,/FILL,LEVELS=levels_set,C_COLORS=colour_set,/IRREGULAR

		ENDIF
	ENDFOR

	IF KEYWORD_SET(erg) THEN BEGIN

		fname='/radar01/work/Aurora/Kiban-S/Ground/Campaign/orb_txt/'+date+'.txt'
		SPAWN,'grep '+erg_tim+' '+fname+' > tmp_erg_orb.txt'

        	la_fot=0.0
        	ma_fot=0.0
        	lo_fot=0.0
        	ml_fot=0.0
        	hm_fot='abc'

        	OPENR,inp_lun,'tmp_erg_orb.txt',/GET_LUN

               	tmp_string='abc'
               	READF,inp_lun,tmp_string
               	hh=FIX(STRMID(tmp_string,11,2))
               	mm=FIX(STRMID(tmp_string,14,2))
               	ss=FIX(STRMID(tmp_string,17,2))
               	hm_fot=STRING(hh,FORMAT='(I2.2)')+STRING(mm,FORMAT='(I2.2)')+' UT'
               	tmp_read=FLTARR(4)
               	READS,STRMID(tmp_string,21,173),tmp_read
               	la_fot=tmp_read(0)
               	lo_fot=tmp_read(1)
		LOADCT,39,/SILENT
                psym8,/FILL
                OPLOT,[lo_fot],[la_fot],PSYM=8,COL=200,SYMSIZE=1.0
                psym8
                OPLOT,[lo_fot],[la_fot],PSYM=8,COL=255,SYMSIZE=1.0
                XYOUTS,lo_fot,la_fot+0.2,'ARASE',/DATA,NOCLIP=0,COL=200,ALI=0.5,CHARSIZE=!P.CHARSIZE
		LOADCT,0,/SILENT
        	FREE_LUN,inp_lun
	ENDIF

; Plot coastlines
	overlay_coast

; Overlay stations
        gla_sta=[69.58,69.76,67.37,67.28,67.84,69.07]
        glo_sta=[19.23,27.01,26.63,20.72,20.42,20.76]
        let_sta=['TRO','KEV','SOD','TJA','KRN','KIL']
        PSYM8,/FILL
	LOADCT,39,/SILENT
        FOR i=0,N_ELEMENTS(let_sta)-1 DO BEGIN
                OPLOT,[glo_sta(i)],[gla_sta(i)],PSYM=8,SYMSIZE=0.5,COL=250
                OPLOT,[glo_sta(i)],[gla_sta(i)],PSYM=1,SYMSIZE=1.0,COL=250
                XYOUTS,glo_sta(i),gla_sta(i)+0.2,let_sta(i),COL=252,ALI=0.5,CHARSIZE=0.6*!P.CHARSIZE
        ENDFOR
	LOADCT,0,/SILENT

; Plot messages
	pos=!P.POSITION
	XYOUTS,pos(0),pos(3)+0.01,date,/NORMAL
	XYOUTS,pos(2),pos(3)+0.01,STRING(hhmmss,FORMAT='(I6.6)')+' UT',/NORMAL,ALI=1

; Plot colour bar
	typ_raw=0
        IF NOT KEYWORD_SET(no_bar) THEN BEGIN
		position=[!P.POSITION(2)+0.01,!P.POSITION(1),!P.POSITION(2)+0.02,!P.POSITION(3)]
		IF typ_raw EQ 1 THEN legend='Absolute Optical Intensity (R)' ELSE legend='Raw Count'
		plot_psa_colour_bar,position=position,legend=legend
	ENDIF

	skip0:

	END

;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	PLOT_PSA_TIF_MULTI_PANEL_ASK
;

	PRO plot_psa_tif_multi_panel_ask,date,hhmmss,abs=abs,alt=alt,position=position,ele_lim=ele_lim,check=check,erg=erg

	COMMON prf_psa

        IF NOT KEYWORD_SET(alt) THEN alt=110

	IF KEYWORD_SET(check) THEN BEGIN

       		SET_PLOT,'X'
		WINDOW,XSIZE=800,YSIZE=800
		dpanel_psa,1,1,0,0,/square
		set_scale_psa,2000,3000
	ENDIF

; The number of cameras and pixel size
	num_cam=1
	pxl=256
	
; Break date
	yr_read=STRMID(date,0,4)
	mo_read=STRMID(date,4,2)
	dy_read=STRMID(date,6,2)
	hms_sel=STRMID(STRING(hhmmss,FORMAT='(I6.6)'),0,5)
	hh_read=STRMID(STRING(hhmmss,FORMAT='(I6.6)'),0,2)
	mm_read=STRMID(STRING(hhmmss,FORMAT='(I6.6)'),2,2)
	erg_tim=yr_read+'-'+mo_read+'-'+dy_read+'/'+hh_read+':'+mm_read+':00'

; Print message
	PRINT,date+': '+STRING(hhmmss,FORMAT='(I6.6)')+' UT'

; Make arrays for images
	avl_map=INTARR(num_cam)
	tif_map=FLTARR(num_cam,pxl,pxl)
        azm_map=FLTARR(num_cam,pxl,pxl)
        ele_map=FLTARR(num_cam,pxl,pxl)
        gla_map=FLTARR(num_cam,pxl,pxl)
        glo_map=FLTARR(num_cam,pxl,pxl)

; Path for the data
	camids=[7]
	FOR n=0,num_cam-1 DO BEGIN
		camid=camids(n)

; Read the data
		path=dir_psa+'cam'+STRING(camid,FORMAT='(I1.1)')+'/'+yr_read+'/'+mo_read+'/'+dy_read+'_QL'
        	fname=path+'/C'+STRING(camid,FORMAT='(I1.1)')+'_'+date+'_'+hms_sel+'*.tif'
		found_files=FILE_SEARCH(fname,COUNT=num_found_files)
		IF num_found_files EQ 1 THEN BEGIN

		avl_map(n)=1
               	img=READ_TIFF(found_files(0))
		tif_map(n,*,*)=ROTATE(img,7)

; Read azimuth and elevation
	        found_file=FILE_SEARCH(dir_psa+'fovd/azm_ele/C'+STRING(camid,FORMAT='(I1.1)')+'*az.txt',COUNT=num_found_file)
        	OPENR,azm_lun,found_file(num_found_file-1),/GET_LUN
	        found_file=FILE_SEARCH(dir_psa+'fovd/azm_ele/C'+STRING(camid,FORMAT='(I1.1)')+'*el.txt',COUNT=num_found_file)
        	OPENR,ele_lun,found_file(num_found_file-1),/GET_LUN
        	FOR i=0,pxl*2-1 DO BEGIN
                	tmp_read=FLTARR(pxl*2)
                	READF,ele_lun,tmp_read
			;ele_map(n,*,pxl-1-i)=tmp_read
			IF i MOD 2 EQ 0 THEN FOR j=0,pxl-1 DO ele_map(n,j,pxl-1-i/2)=tmp_read(2*j)
                	READF,azm_lun,tmp_read
			;azm_map(n,*,pxl-1-i)=tmp_read
			IF i MOD 2 EQ 0 THEN FOR j=0,pxl-1 DO azm_map(n,j,pxl-1-i/2)=tmp_read(2*j)
        	ENDFOR
        	FREE_LUN,azm_lun
        	FREE_LUN,ele_lun

; Read glat and glon info
        	found_file=FILE_SEARCH(dir_psa+'fovd/gla_glo/C'+STRING(camid,FORMAT='(I1.1)')+'*'+ $
                	STRING(alt,FORMAT='(I3.3)')+'_lat.txt',COUNT=num_found_file)
        	OPENR,gla_lun,found_file(num_found_file-1),/GET_LUN
        	found_file=FILE_SEARCH(dir_psa+'fovd/gla_glo/C'+STRING(camid,FORMAT='(I1.1)')+'*'+ $
                	STRING(alt,FORMAT='(I3.3)')+'_long.txt',COUNT=num_found_file)
        	OPENR,glo_lun,found_file(num_found_file-1),/GET_LUN
        	FOR i=0,pxl*2-1 DO BEGIN
                	tmp_read=FLTARR(pxl*2)
                	READF,gla_lun,tmp_read
			;gla_map(n,*,pxl-1-i)=tmp_read
			IF i MOD 2 EQ 0 THEN FOR j=0,pxl-1 DO gla_map(n,j,pxl-1-i/2)=tmp_read(2*j)
			;IF camid EQ 6 THEN gla_map(n,*,i)=tmp_read ELSE gla_map(n,*,pxl-1-i)=tmp_read
                	READF,glo_lun,tmp_read
			;glo_map(n,*,pxl-1-i)=tmp_read
			IF i MOD 2 EQ 0 THEN FOR j=0,pxl-1 DO glo_map(n,j,pxl-1-i/2)=tmp_read(2*j)
			;IF camid EQ 6 THEN glo_map(n,*,i)=tmp_read ELSE glo_map(n,*,pxl-1-i)=tmp_read
        	ENDFOR
        	FREE_LUN,gla_lun
        	FREE_LUN,glo_lun

		ENDIF

        ENDFOR

; Plot area
	IF KEYWORD_SET(position) THEN !P.POSITION=position
	pos=!P.POSITION

; Set up the elevation angle limitation
	IF NOT KEYWORD_SET(ele_lim) THEN ele_lim=15.0

; Set up the plotting area
	IF NOT KEYWORD_SET(xrange) THEN BEGIN
		;xrange=[200,230]
		xrange=[205,225]
		;IF cid_raw EQ 2 THEN xrange=[14,38]
	ENDIF
	IF NOT KEYWORD_SET(yrange) THEN BEGIN
		;yrange=[57,70]
		yrange=[58,67]
		;IF cid_raw EQ 2 THEN yrange=[63,72]
	ENDIF

; Plot frame
        PLOT,[0],[0],XRANGE=xrange,YRANGE=yrange,/XSTYLE,/YSTYLE,XTITLE='Geo. Lon.',YTITLE='Geo. Lat'

	FOR n=0,num_cam-1 DO BEGIN

		IF avl_map(n) EQ 1 THEN BEGIN

		camid=camids(n)

; Prepare array for contouring
		gla_psa=REFORM(gla_map(n,*,*))
		glo_psa=REFORM(glo_map(n,*,*))
		azm_psa=REFORM(azm_map(n,*,*))
		ele_psa=REFORM(ele_map(n,*,*))
		tmp_img=REFORM(tif_map(n,*,*))
        	tmp_bin=WHERE(gla_psa NE 99999.9 AND gla_psa NE 99999.9 AND ele_psa GE ele_lim,no_tmp_bin)
        	x_cont=MAKE_ARRAY(no_tmp_bin,/FLOAT,VALUE=0.0)
        	y_cont=MAKE_ARRAY(no_tmp_bin,/FLOAT,VALUE=0.0)
        	z_cont=MAKE_ARRAY(no_tmp_bin,/FLOAT,VALUE=0.0)

; Put values into the arrays for contouring
        	num_cnt=0L

		IF camid EQ 7 THEN BEGIN
        	FOR x=0,255 DO BEGIN
                	FOR y=0,255 DO BEGIN
                        	IF gla_psa(x,y) NE 99999.9 AND glo_psa(x,y) NE 99999.9 AND ele_psa(x,y) GE ele_lim $
					;AND -0.0350*(glo_psa(x,y)-32.72)+68.26 GE gla_psa(x,y) $
					THEN BEGIN

                                	z_cont(num_cnt)=FLOAT(tmp_img(x,y));*COS(zan_psa(x,y)*!DTOR)
                                	IF z_cont(num_cnt) GT max_psa THEN z_cont(num_cnt)=max_psa
                                	IF z_cont(num_cnt) LT min_psa THEN z_cont(num_cnt)=min_psa

                                	x_cont(num_cnt)=glo_psa(x,y)
                                	y_cont(num_cnt)=gla_psa(x,y)

                                	num_cnt++
                        	ENDIF
                	ENDFOR
        	ENDFOR
		x_cont=x_cont(0:num_cnt-1)
		y_cont=y_cont(0:num_cnt-1)
		z_cont=z_cont(0:num_cnt-1)
		ENDIF

; Contouring
        	contour_lines=30
        	levels_set=(max_psa-min_psa)*FINDGEN(contour_lines)/FLOAT(contour_lines)+min_psa
        	colour_set=!D.TABLE_SIZE*INDGEN(contour_lines)/contour_lines
        	CONTOUR,z_cont,x_cont,y_cont,/OVERPLOT,/FILL,LEVELS=levels_set,C_COLORS=colour_set,/IRREGULAR

		ENDIF
	ENDFOR

	IF KEYWORD_SET(erg) THEN BEGIN

		fname='/radar01/work/Aurora/Kiban-S/Ground/Campaign/orb_txt/'+date+'.txt'
		SPAWN,'grep '+erg_tim+' '+fname+' > tmp_erg_orb.txt'

        	la_fot=0.0
        	ma_fot=0.0
        	lo_fot=0.0
        	ml_fot=0.0
        	hm_fot='abc'

        	OPENR,inp_lun,'tmp_erg_orb.txt',/GET_LUN

               	tmp_string='abc'
               	READF,inp_lun,tmp_string
               	hh=FIX(STRMID(tmp_string,11,2))
               	mm=FIX(STRMID(tmp_string,14,2))
               	ss=FIX(STRMID(tmp_string,17,2))
               	hm_fot=STRING(hh,FORMAT='(I2.2)')+STRING(mm,FORMAT='(I2.2)')+' UT'
               	tmp_read=FLTARR(4)
               	READS,STRMID(tmp_string,21,173),tmp_read
               	la_fot=tmp_read(0)
               	lo_fot=tmp_read(1)
		LOADCT,39,/SILENT
                psym8,/FILL
                OPLOT,[lo_fot],[la_fot],PSYM=8,COL=200,SYMSIZE=1.0
                psym8
                OPLOT,[lo_fot],[la_fot],PSYM=8,COL=255,SYMSIZE=1.0
                XYOUTS,lo_fot,la_fot+0.2,'ARASE',/DATA,NOCLIP=0,COL=200,ALI=0.5,CHARSIZE=!P.CHARSIZE
		LOADCT,0,/SILENT
        	FREE_LUN,inp_lun
	ENDIF

; Plot coastlines
	overlay_coast

; Overlay stations
        gla_sta=[69.58,69.76,67.37,67.28,66.2,64.7,54.7,65.13,-69.0,62.34,70.03,66.78,62.39,49.41,56.6,60.05,67.84]
        glo_sta=[19.23,27.01,26.63,20.72,342.9,338.9,246.7,212.51,39.59,25.51,88.01,123.37,214.78,277.50,298.3,150.73,20.42]
        let_sta=['TRO','KEV','SOD','TJA','TJO','HUS','ATH','PFRR','SYO','NYR','ITK','ZGN','GAK','KAP','NAI','MGD','KIR']
        PSYM8,/FILL
	LOADCT,39,/SILENT
        FOR i=0,N_ELEMENTS(let_sta)-1 DO BEGIN
                OPLOT,[glo_sta(i)],[gla_sta(i)],PSYM=8,SYMSIZE=0.5,COL=250
                OPLOT,[glo_sta(i)],[gla_sta(i)],PSYM=1,SYMSIZE=1.0,COL=250
                XYOUTS,glo_sta(i),gla_sta(i)+0.2,let_sta(i),COL=252,ALI=0.5,CHARSIZE=0.6*!P.CHARSIZE
        ENDFOR
	LOADCT,0,/SILENT

; Plot messages
	pos=!P.POSITION
	XYOUTS,pos(0),pos(3)+0.01,date,/NORMAL
	XYOUTS,pos(2),pos(3)+0.01,STRING(hhmmss,FORMAT='(I6.6)')+' UT',/NORMAL,ALI=1

; Plot colour bar
	typ_raw=0
        IF NOT KEYWORD_SET(no_bar) THEN BEGIN
		position=[!P.POSITION(2)+0.01,!P.POSITION(1),!P.POSITION(2)+0.02,!P.POSITION(3)]
		IF typ_raw EQ 1 THEN legend='Absolute Optical Intensity (R)' ELSE legend='Raw Count'
		plot_psa_colour_bar,position=position,legend=legend
	ENDIF

	skip0:

	END
;------------------------------------------------------------------------------
;
;  NAME:
;
;        PLOT_GEOGRA_MAP
;
;   
;
;  THEMIS > PLOT_GEOGRA_MAP
;  THEMIS > plot_geogra_map,'20161004'(day:YYYYMMDD),'2300'(start
;  time:HHMM),2(duration:hour),10E6(scale of map)


        PRO  PLOT_GEOGRA_MAP,day,hhmmss,dura,scale,col,count=count,no_bar=no_bar,no_change=no_change,no_plot=no_plot
  
   	COMMON prf_psa
	COMMON fov_psa
        COMMON raw_psa
        COMMON cor_psa
        COMMON struct_psa
        COMMON apa_psa
        COMMON apa_fov

; set time span
     IF NOT KEYWORD_SET(day) THEN BEGIN
        day=''
        hhmmss=''
     
        READ,day,PROMPT='Enter day(YYYYMMDD):'
        READ,hhmmss,PROMPT='Enter start time(HHMMSS):'
        READ,dura,PROMPT='Enter duration(hour):'
        
     ENDIF

     
     IF NOT KEYWORD_SET(count) THEN count=0L
     IF NOT KEYWORD_SET(scale) THEN scale=10E6
     IF NOT KEYWORD_SET(col) THEN col='1'
     IF NOT KEYWORD_SET(YYYY_orb) THEN BEGIN

        YYYY_orb=STRMID(day,0,4)
        MM_orb=STRMID(day,4,2)
        DD_orb=STRMID(day,6,2)
        hour_orb=STRMID(hhmmss,0,2)
        min_orb=STRMID(hhmmss,2,2)
        dura_orb=dura
     ENDIF
     
     YYYY=STRMID(day,0,4)
     MM=STRMID(day,4,2)
     DD=STRMID(day,6,2)
     hhmm=STRMID(hhmmss,0,4)
     hour=STRMID(hhmmss,0,2)
     min=STRMID(hhmmss,2,2)
     second=STRMID(hhmmss,4,2)

     
     timespan,YYYY+'-'+MM+'-'+DD+'/'+hour+':'+min+':00',dura,/hour

     IF NOT KEYWORD_SET(img_raw) THEN file_psa_raw,2,day,hhmm,1
     
     GO_TIME_PSA,LONG(hhmmss)
     
     num_sod=0L
     IF NOT KEYWORD_SET(no_plot) THEN BEGIN
        !P.POSITION=[0,0,0,0.97]

        map2d_init
        map2d_coord,'geo'
        window,0,xsize=800,ysize=640 & erase
        map2d_set,glatc=67.58,glonc=33.39,scale=scale

        gla_sta=[67.25]
        glo_sta=[26.35]
        let_sta=['SOD']
        gla_sta_apa=[67.58]
        glo_sta_apa=[33.39]
        let_sta_apa=['APA']
     ENDIF

     IF NOT KEYWORD_SET(ele_lim) THEN ele_lim=15.0
     ele_lim=13.0
        tmp_bin=WHERE(gla_apa NE 99999.9 AND gla_apa NE 99999.9 AND ele_apa GE ele_lim,no_tmp_bin)
        x_cont=MAKE_ARRAY(no_tmp_bin,/FLOAT,VALUE=0.0)
        y_cont=MAKE_ARRAY(no_tmp_bin,/FLOAT,VALUE=0.0)
        z_cont=MAKE_ARRAY(no_tmp_bin,/FLOAT,VALUE=0.0)

; Put values into the arrays for contouring
        num_cnt=0L
        max_apa=3400
        min_apa=500
        FOR x=0,382 DO BEGIN
                FOR y=0,286 DO BEGIN
                        IF gla_apa(x,y) NE 99999.9 AND glo_apa(x,y) NE 99999.9 AND ele_apa(x,y) GE ele_lim THEN BEGIN

                                z_cont(num_cnt)=FLOAT(img_apa(now_apa,x,y));*COS(zan_psa(x,y)*!DTOR)
                                IF z_cont(num_cnt) GT max_apa THEN z_cont(num_cnt)=max_apa
                                IF z_cont(num_cnt) LT min_apa THEN z_cont(num_cnt)=min_apa
                                IF glo_apa(x,y) LE 30 THEN z_cont(num_cnt)='NaN'
                                
                                x_cont(num_cnt)=glo_apa(x,y)
                                y_cont(num_cnt)=gla_apa(x,y)

                                num_cnt++
                        ENDIF
                ENDFOR
        ENDFOR

        
; Contouring

        IF col EQ '0' THEN loadct,[0]
        IF col EQ '1' THEN loadct,[39]
        contour_lines=30
        levels_set=(max_apa-min_apa)*FINDGEN(contour_lines)/FLOAT(contour_lines)+min_apa
        colour_set=!D.TABLE_SIZE*INDGEN(contour_lines)/contour_lines
        CONTOUR,z_cont,x_cont,y_cont,/OVERPLOT,/CELL_FILL,LEVELS=levels_set,C_COLORS=colour_set,/IRREGULAR
        
; Prepare array for contouring
        IF NOT KEYWORD_SET(ele_lim) THEN ele_lim=15.0
        ele_lim=13.0
        tmp_bin=WHERE(gla_psa NE 99999.9 AND gla_psa NE 99999.9 AND ele_psa GE ele_lim,no_tmp_bin)
        x_cont=MAKE_ARRAY(no_tmp_bin,/FLOAT,VALUE=0.0)
        y_cont=MAKE_ARRAY(no_tmp_bin,/FLOAT,VALUE=0.0)
        z_cont=MAKE_ARRAY(no_tmp_bin,/FLOAT,VALUE=0.0)

; Put values into the arrays for contouring
        num_cnt=0L
        ;img_raw_smo=FLTARR(256,256)
        
        ;FOR k=0,tr_img-1 DO BEGIN
           
        ;  IF now_raw+k GE num_raw THEN BEGIN
        ;     k=k-1
        ;     BREAK
        ;   ENDIF

        ;   img_raw_smo(*,*)+=img_raw(now_raw+k,*,*)
        ;ENDFOR
        ;img_raw_smo=img_raw_smo/k

        IF (count EQ 0) AND (NOT KEYWORD_SET(no_change)) THEN BEGIN
           img_changed=CHANGE_TIME_RES_IMG(img_raw)
        ENDIF  
           img_changed_smo=REFORM(img_changed(count,*,*))
        
        FOR x=0,255 DO BEGIN
                FOR y=0,255 DO BEGIN
                        IF gla_psa(x,y) NE 99999.9 AND glo_psa(x,y) NE 99999.9 AND ele_psa(x,y) GE ele_lim THEN BEGIN

                                z_cont(num_cnt)=FLOAT(img_changed_smo(x,y));*COS(zan_psa(x,y)*!DTOR)
                                IF z_cont(num_cnt) GT max_psa THEN z_cont(num_cnt)=max_psa
                                IF z_cont(num_cnt) LT min_psa THEN z_cont(num_cnt)=min_psa
                                IF glo_psa(x,y) GE 30.5 THEN z_cont(num_cnt)='NaN'

                                x_cont(num_cnt)=glo_psa(x,y)
                                y_cont(num_cnt)=gla_psa(x,y)

                                num_cnt++
                        ENDIF
                ENDFOR
        ENDFOR

        
; Contouring
        IF NOT KEYWORD_SET(no_plot) THEN BEGIN
        IF col EQ '0' THEN loadct,[0]
        IF col EQ '1' THEN loadct,[39]
        contour_lines=30
        levels_set=(max_psa-min_psa)*FINDGEN(contour_lines)/FLOAT(contour_lines)+min_psa
        colour_set=!D.TABLE_SIZE*INDGEN(contour_lines)/contour_lines
        CONTOUR,z_cont,x_cont,y_cont,/OVERPLOT,/CELL_FILL,LEVELS=levels_set,C_COLORS=colour_set,/IRREGULAR
        
; Get orbit data
        timespan,YYYY_orb+'-'+MM_orb+'-'+DD_orb+'/'+hour_orb+':'+min_orb+':00',dura_orb,/hour
        get_data,'rbspa_ifoot_geo_lat',data=data
        IF NOT KEYWORD_SET(data) THEN GET_ORBIT_DATA,YYYY_orb+MM+DD,hour_orb+min_orb,dura_orb
        
; overlay satelite orbit
        overlay_map_coast
        overlay_map_sc_ifoot,'rbspa_ifoot_geo_lat','rbspa_ifoot_geo_lon'
        timespan,YYYY+'-'+MM+'-'+DD+'/'+hour+':'+min+':00',dura,/hour
        
; plot station and satelite position
        
        PSYM8,/FILL
     	OPLOT,[glo_sta],[gla_sta],PSYM=3,SYMSIZE=0.5,COL=fsc_color("black")
        OPLOT,[glo_sta],[gla_sta],PSYM=1,SYMSIZE=1.0,COL=fsc_color("black")
        
        OPLOT,[glo_sta_apa],[gla_sta_apa],PSYM=3,SYMSIZE=0.5,COL=fsc_color("black")
        OPLOT,[glo_sta_apa],[gla_sta_apa],PSYM=1,SYMSIZE=1.0,COL=fsc_color("black")
        PLOT_POS_FOOTPRINT
        XYOUTS,glo_sta,gla_sta+0.2,let_sta,COL=fsc_color("black"),ALI=0.5,CHARSIZE=0.8*!P.CHARSIZE
        XYOUTS,glo_sta_apa,gla_sta_apa+0.2,let_sta_apa,COL=fsc_color("black"),ALI=0.5,CHARSIZE=0.8*!P.CHARSIZE
; Plot messages
        IF KEYWORD_SET(no_bar) THEN BEGIN
           day_mes=YYYY+'-'+MM+'-'+DD
           ti_mes=hour+':'+min+':'+second+' '+'5 min'
        ENDIF ELSE BEGIN
        
           day_mes=STRMID(time_string(img_changed_tim_ms(count)),0,10)
           ti_mes=STRMID(time_string(img_changed_tim_ms(count)),11,8)+'('+STRING(ms_changed(count),FORMAT='(I3.3)')+')'
        ENDELSE
        
        pos=!P.POSITION
        XYOUTS,pos(0)+0.05,pos(3)-0.01,day_mes+' '+ti_mes+': Cam'+STRING(cid_raw,FORMAT='(I1.1)'),/NORMAL,charsize=2.0

; Plot colour bar
        typ_raw=0
        IF NOT KEYWORD_SET(no_bar) THEN BEGIN
		position=[!P.POSITION(2)+0.01,!P.POSITION(1)+0.02,!P.POSITION(2)+0.02,!P.POSITION(3)]
		IF typ_raw EQ 1 THEN legend='Absolute Optical Intensity (R)' ELSE legend='Raw Count'
		plot_psa_colour_bar,position=position,legend=legend
	ENDIF

     ENDIF
        
        day_mes=STRMID(time_string(img_changed_tim_ms(count)),0,10)
        ti_mes=STRMID(time_string(img_changed_tim_ms(count)),11,8)+'('+STRING(ms_changed(count),FORMAT='(I3.3)')+')'
        geogra_map=CREATE_STRUCT(NAME='geogra_map','x',x_cont,'y',y_cont,'z',z_cont,'mes',day_mes+' '+ti_mes,'pos',!P.POSITION)     

   END
;------------------------------------------------------------------------------
;
;  NAME:
;
;        SAVE_MERGED_MAP
;
;   
;
;  THEMIS > PLOT_GEOGRA_MAP
;  THEMIS > plot_geogra_map,'20161004'(day:YYYYMMDD),'2300'(start
;  time:HHMM),2(duration:hour),10E6(scale of map)

        PRO save_merged_map
        
        COMMON prf_psa
	COMMON fov_psa
        COMMON raw_psa
        COMMON cor_psa
        COMMON struct_psa
        COMMON apa_psa
        COMMON apa_fov

; set time span
     IF NOT KEYWORD_SET(day) THEN BEGIN
        day=''
        hhmmss=''
     
        READ,day,PROMPT='Enter day(YYYYMMDD):'
        READ,hhmmss,PROMPT='Enter start time(HHMMSS):'
        READ,dura,PROMPT='Enter duration(hour):'
        
     ENDIF

     
     IF NOT KEYWORD_SET(count) THEN count=0L
     IF NOT KEYWORD_SET(scale) THEN scale=10E6
     IF NOT KEYWORD_SET(col) THEN col='1'
     

     ;IF NOT KEYWORD_SET(img_apa) THEN file_psa_apa,path='/home/suguru/spedas/image_data/russia/apa_txt_23/'
     min=STRMID(hhmmss,2,2)
     i=0L
     WHILE hhmmss NE '005959' DO BEGIN
        IF hhmmss EQ '000000' THEN file_psa_apa,path='/home/suguru/spedas/image_data/russia/apa_txt_00/'
        IF LONG(STRMID(hhmmss,2,2)) - FIX(min) NE 0 THEN BEGIN
           FREE_VARIABLE
           count=0L
        ENDIF
        plot_geogra_map,day,hhmmss,dura,count=count
        makepng,'/home/suguru/spedas/image_data/merged_geo/'+day+hhmmss
; Set next time
        min=STRMID(hhmmss,2,2)

        hhmmss = STRING(LONG(hhmmss)+1,FORMAT='(I6.6)')
        IF STRMID(hhmmss,4,2) EQ '60' THEN BEGIN
           hhmmss = STRING(LONG(hhmmss)+100-60,FORMAT='(I6.6)')                
        ENDIF              
        
        IF STRMID(hhmmss,2,2) EQ '60' THEN BEGIN
           hhmmss = STRING(LONG(hhmmss)+10000-6000,FORMAT='(I6.6)')
        ENDIF
        
        IF STRMID(hhmmss,0,2) EQ '24' THEN BEGIN
           day = STRING(LONG(day)+1,FORMAT='(I8.8)')
           hhmmss = '000000'
        ENDIF
        hhmm=STRMID(hhmmss,0,4)
        i++
        count=count+8
        go_apa,i
     ENDWHILE
     END
        
     
;-------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       PLOT_POS_FOOTPRINT
;
; NOTES:
;
;
;
        
        
        PRO PLOT_POS_FOOTPRINT

          COMMON raw_psa
          
           get_data,'rbspa_ifoot_geo_lat',data=lat_data
           get_data,'rbspa_ifoot_geo_lon',data=lon_data
           now_tim=STRING(yr_raw(now_raw),FORMAT='(I4.4)')+'-'+STRING(mo_raw(now_raw),FORMAT='(I2.2)')+'-'+STRING(dy_raw(now_raw),FORMAT='(I2.2)')+'/'+STRING(hh_raw(now_raw),FORMAT='(I2.2)')+':'+STRING(mm_raw(now_raw),FORMAT='(I2.2)')+':30'
           
           k=WHERE(time_string(lat_data.x) EQ now_tim)

           k=k(0)

           OPLOT,[lon_data.y(k)],[lat_data.y(k)],PSYM=8,SYMSIZE=1.0,COL=fsc_color("black")
           XYOUTS,lon_data.y(k),lat_data.y(k)+0.2,'RBSPA',COL=fsc_color("black"),ALI=0.5,CHARSIZE=0.8*!P.CHARSIZE
        END


;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	GET_ORBIT_DATA
;
         PRO get_orbit_data,day,hhmm,dura
  
           IF NOT KEYWORD_SET(day) THEN BEGIN
              day=''
              hhmm=''
                 
              READ,day,PROMPT='Enter day(YYYYMMDD):'
              READ,hhmm,PROMPT='Enter start time(HHMM):'
              READ,dura,PROMPT='Enter duration(hour):'
           ENDIF
           YYYY=STRMID(day,0,4)
           MM=STRMID(day,4,2)
           DD=STRMID(day,6,2)
           hour=STRMID(hhmm,0,2)
           min=STRMID(hhmm,2,2)
           
           ;timespan,YYYY+'-'+MM+'-'+DD+'/'+hour+':'+min+':00',dura,/hour
 
           map2d_init
           map2d_coord,'geo'
           rbsp_ifoot,probe='a',autocalc_parmod='autocalc_parmod'
           get_data,'rbspa_ifoot_geo_lat',data=lat_data
           get_data,'rbspa_ifoot_geo_lon',data=lon_data

      END
       
;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	save_png_raw
;
        PRO save_png_raw
          
           COMMON prf_psa
           COMMON fov_psa
           COMMON raw_psa
           COMMON cor_psa
           COMMON struct_psa

; set time range 
           day=''
           hhmmss=''
     
           READ,day,PROMPT='Enter day(YYYYMMDD):'
           READ,hhmmss,PROMPT='Enter start time(HHMMSS):'
           READ,dura,PROMPT='Enter duration(hour):'

           YYYY=STRMID(day,0,4)
           MM=STRMID(day,4,2)
           DD=STRMID(day,6,2)
           hhmm=STRMID(hhmmss,0,4)
           hour=STRMID(hhmmss,0,2)
           min=STRMID(hhmmss,2,2)
           tr_change=8.0
           count=0L
           
; plot map and save

           WHILE (1) DO BEGIN

              hhmm=STRMID(hhmmss,0,4)              
              file_psa_raw,2,day,hhmm,1
              num_img=N_ELEMENTS(img_raw(*,0,0))*(tr_change/tr_img)

; plot map
              FOR i=0,num_img-1 DO BEGIN
                 PLOT_GEOGRA_MAP,day,hhmmss,dura,count=i
; save img as png
                 print,img_changed_tim_ms(i)
                 tmp_yr=STRMID(img_changed_tim_ms(i),0,4)
                 tmp_mo=STRMID(img_changed_tim_ms(i),5,2)
                 tmp_dy=STRMID(img_changed_tim_ms(i),8,2)
                 tmp_tim=STRMID(img_changed_tim_ms(i),11,12)
                 makepng,'/home/suguru/spedas/image_data/coh_sod/' 
                 makepng,'/home/suguru/spedas/image_data/img_raw_8Hz/'+geogra_map.mes
              ENDFOR

              i=0L
; Set next time range
              hhmmss = STRING(LONG(hhmmss)+100,FORMAT='(I6.6)')

              IF STRMID(hhmmss,4,2) EQ '60' THEN BEGIN
                 hhmmss = STRING(LONG(hhmmss)+100-60,FORMAT='(I6.6)')
              ENDIF
              
              
              IF STRMID(hhmmss,2,2) EQ '60' THEN BEGIN
                 hhmmss = STRING(LONG(hhmmss)+10000-6000,FORMAT='(I6.6)')
              ENDIF

              IF STRMID(hhmmss,0,2) EQ '24' THEN BEGIN
                 day = STRING(LONG(day)+1,FORMAT='(I8.8)')
                 hhmmss = '000000'
              ENDIF
              
              FREE_VARIABLE
              
           ENDWHILE
           
        END
        

;-------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       RANGE_TO_COORDS
;
; NOTES:
;
;       Given a start "lat" and "lon" find the position of the point "range" in km away
;       at an azimuth of "azimuth" from north.
;

        FUNCTION range_to_coords,lat,lon,azimuth,range

        IF azimuth GT 180 THEN az=azimuth-360 ELSE az=azimuth
        Re=6370.0
        coLat=90-lat

        coLat_point=ACOS(COS(range/Re)*COS(coLat*!DTOR)+SIN(range/Re)* $
                SIN(coLat*!DTOR)*COS(az*!DTOR))*!RADEG
        lon_point=ACOS(MIN([(COS(range/Re)-COS(coLat_point*!DTOR)* $
                COS(coLat*!DTOR))/(SIN(coLat_point*!DTOR)*SIN(coLat*!DTOR)),1.0]))*!RADEG
        lat_point=90-coLat_point
        IF az GT 0 THEN lon_point=lon+lon_point ELSE lon_point=lon-lon_point

        RETURN,[lat_point,lon_point]

        END

;-------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       ANG_TO_LOC
;
; NOTE:
;
;       Given an ground position (geog coords) and azimuth and the zenith of
;       observations and the assumed emission altitude return the geographic
;       lat and lon (i.e. [ele,azm,alt] -> [lat,lon])
;

        FUNCTION ang_to_loc,lat,lon,azimuth,zenith,altitude

; Earth radius
        Re=6370.0
        
; Assume that msp look direction maps to within +/-25 degrees of the station
; and determine heights that projection maps to within this range of lats
        delta_lat=25*FINDGEN(1000)/1000
        h=Re*(TAN((90-ABS(zenith)+delta_lat)*!pi/180)*SIN(delta_lat*!pi/180)+COS(delta_lat*!pi/180)-1)

; Find h that corresponds best to desired altitude
        min_dev=MIN(ABS(altitude-h),min_pos)

; Convert position to a distance
        range=Re*delta_lat(min_pos)*!pi/180

; Determine latitude and longitude of distance from msp location at given azimuth
        IF zenith LT 0 THEN az_point=azimuth+180 ELSE az_point=azimuth
        
        RETURN,range_to_coords(lat,lon,az_point,range)

        END

;----------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	OVERLAY_IMAGINARY_FOV
;

	PRO overlay_imaginary_fov,lat,lon,zan_lim=zan_lim,height=height,col=col,thick=thick

	IF NOT KEYWORD_SET(zan_lim) THEN zan_lim=75.0
	IF NOT KEYWORD_SET(height) THEN height=110.0

	lat_plot=FLTARR(360)
	lon_plot=FLTARR(360)
	FOR i=0,359 DO BEGIN
        	tmp=ang_to_loc(lat,lon,i,zan_lim,height)
		lat_plot(i)=tmp(0)
		lon_plot(i)=tmp(1)
	ENDFOR

	IF NOT KEYWORD_SET(col) THEN col=253
	IF NOT KEYWORD_SET(thick) THEN thick=3
	OPLOT,lon_plot,lat_plot,COL=col,THICK=thick
     END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	 normalize_routine
;

        FUNCTION normalize_routine,array

          IF MAX(ABS(array)) EQ 0 THEN BEGIN
             ret_array=DBLARR(N_ELEMENTS(array))
             RETURN,ret_array
          ENDIF
          
          
           ret_array=DBLARR(N_ELEMENTS(array))
           FOR i=0,N_ELEMENTS(array)-1 DO BEGIN
              ret_array(i)=(array(i)-min(array))/(max(array)-min(array))
           ENDFOR

           RETURN,ret_array

        END
;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	 SAVE_PNG_COR
;
        PRO save_png_cor

           COMMON prf_psa
           COMMON fov_psa
           COMMON raw_psa
           COMMON cor_psa

; input time info
           IF NOT KEYWORD_SET(day) THEN BEGIN
              day=''
              hhmmss=''
              
              READ,day,PROMPT='Enter day(YYYYMMDD):'
              READ,hhmmss,PROMPT='Enter start time(HHMMSS):'
              READ,dura,PROMPT='Enter duration(hour):'
              READ,span,PROMPT='Enter calculation span of CC(second):'

           ENDIF
              
              YYYY=STRMID(day,0,4)
              MM=STRMID(day,4,2)
              DD=STRMID(day,6,2)
              hhmm=STRMID(hhmmss,0,4)
              hour=STRMID(hhmm,0,2)
              min=STRMID(hhmm,2,2)
              second=STRMID(hhmmss,4,2)
              ;count=0

              file_psa_raw,2,day,hhmm,span/60+2
              FOR i=0,(dura*60)/(cal_span/60)-1 DO BEGIN

                 IF FIX(STRMID(hhmmss,2,2))-FIX(min) EQ 1 THEN BEGIN
                    FREE_VARIABLE
                 
                    file_psa_raw,2,day,hhmm,span/60+2
                 ENDIF
                 

                 PRINT,'Current time: '+day+'-'+hhmmss
              
; Plot correlation coefficient with contour
                 PLOT_COR_CONTOUR,day,hhmmss,dura,span              
                
; Set next time range
                 min=STRMID(hhmmss,2,2)
                 hhmmss = STRING(LONG(hhmmss)+cal_span,FORMAT='(I6.6)')
                 IF STRMID(hhmmss,4,2) EQ '60' THEN BEGIN
                    hhmmss = STRING(LONG(hhmmss)+100-60,FORMAT='(I6.6)')
                 ENDIF
                 
              
                    IF STRMID(hhmmss,2,2) EQ '60' THEN BEGIN
                       hhmmss = STRING(LONG(hhmmss)+10000-6000,FORMAT='(I6.6)')
                    ENDIF

                    IF STRMID(hhmmss,0,2) EQ '24' THEN BEGIN
                       day = STRING(LONG(day)+1,FORMAT='(I8.8)')
                       hhmmss = '000000'
                    ENDIF

                    hhmm=STRMID(hhmmss,0,4)
                    ;count=count+tr_img
                    ;GO_PSA_RAW,count
                 ENDFOR
        END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	 SAVE_PNG_COR_LINE
;
        PRO save_png_cor_line

           COMMON prf_psa
           COMMON fov_psa
           COMMON raw_psa
           COMMON cor_psa

; input time info
           IF NOT KEYWORD_SET(day) THEN BEGIN
              day=''
              hhmmss=''
              
              READ,day,PROMPT='Enter day(YYYYMMDD):'
              READ,hhmmss,PROMPT='Enter start time(HHMMSS):'
              READ,dura,PROMPT='Enter duration(hour):'
              READ,span,PROMPT='Enter calculation span of CC(second):'

           ENDIF
              
              YYYY=STRMID(day,0,4)
              MM=STRMID(day,4,2)
              DD=STRMID(day,6,2)
              hhmm=STRMID(hhmmss,0,4)
              hour=STRMID(hhmm,0,2)
              min=STRMID(hhmm,2,2)
              second=STRMID(hhmmss,4,2)
              ;count=0

              file_psa_raw,2,day,hhmm,span/60+2
              FOR i=0,(dura*60)/(cal_span/60)-1 DO BEGIN

                 IF FIX(STRMID(hhmmss,2,2))-FIX(min) EQ 1 THEN BEGIN
                    FREE_VARIABLE
                 
                    file_psa_raw,2,day,hhmm,span/60+2
                 ENDIF
                 

                 PRINT,'Current time: '+day+'-'+hhmmss
              
; Plot correlation coefficient with contour
                 PLOT_COR_LINE,day,hhmmss,dura,span
                 makepng,'/home/suguru/spedas/image_data/cor_line/'+day+hhmmss+'(5 min)'
                
; Set next time range
                 min=STRMID(hhmmss,2,2)
                 hhmmss = STRING(LONG(hhmmss)+cal_span,FORMAT='(I6.6)')
                 IF STRMID(hhmmss,4,2) EQ '60' THEN BEGIN
                    hhmmss = STRING(LONG(hhmmss)+100-60,FORMAT='(I6.6)')
                 ENDIF
                 
              
                    IF STRMID(hhmmss,2,2) EQ '60' THEN BEGIN
                       hhmmss = STRING(LONG(hhmmss)+10000-6000,FORMAT='(I6.6)')
                    ENDIF

                    IF STRMID(hhmmss,0,2) EQ '24' THEN BEGIN
                       day = STRING(LONG(day)+1,FORMAT='(I8.8)')
                       hhmmss = '000000'
                    ENDIF

                    hhmm=STRMID(hhmmss,0,4)
                    ;count=count+tr_img
                    ;GO_PSA_RAW,count
                 ENDFOR
        END
       
              
;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	 ROLLING_THUNDER
;
; example:  PsA > file_psa_raw,2,day,hhmm,dura
;           PsA > cor=ROLLING_THUNDER(day,hhmmss,dura,span)
;
         FUNCTION rolling_thunder,day,hhmmss,dura,span

           COMMON prf_psa
           COMMON fov_psa
           COMMON raw_psa
           COMMON cor_psa


; input time info
           IF NOT KEYWORD_SET(day) THEN BEGIN
              day=''
              hhmmss=''
              
              READ,day,PROMPT='Enter day(YYYYMMDD):'
              READ,hhmmss,PROMPT='Enter start time(HHMMSS):'
              READ,dura,PROMPT='Enter duration(hour):'
              READ,span,PROMPT='Enter calculation span of CC(second):'
              
           ENDIF
           
              YYYY=STRMID(day,0,4)
              MM=STRMID(day,4,2)
              DD=STRMID(day,6,2)
              hhmm=STRMID(hhmmss,0,4)
              hour=STRMID(hhmm,0,2)
              min=STRMID(hhmm,2,2)
              second=STRMID(hhmmss,4,2)

           timespan,YYYY+'-'+MM+'-'+DD+'/'+hour+':'+min+':'+second,dura,/hour

           IF NOT KEYWORD_SET(num_raw) THEN file_psa_raw,2,day,hhmm,span/60+2
           print,N_ELEMENTS(img_raw(*,0,0))*(8.0/tr_img)-2*6000*(8.0/tr_img)
; preparing array           
           cor=FLTARR(256,256)
           img_changed=CHANGE_TIME_RES_IMG(img_raw)
           ;img_changed=DBLARR(span+120,256,256)
           ;img_changed_tim=STRARR(span+120)
           ;img_changed_tim_ms=STRARR(span+120)
           ;ts_raw_changed=DBLARR(span)
           tr_change=8.0
           num_img=N_ELEMENTS(img_raw(*,0,0))*(tr_change/tr_img)-2*6000*(tr_change/tr_img)
           lag=[-2,-1,0,1,2]
           max_cor=1.0
           min_cor=0.0
           max_pixel=INTARR(2)
           max_pixel_cor=0.0
           xlim_pxl=256
           ylim_pxl=256
           
; get fbk data of rbspa
           get_data,'rbspa_efw_fbk_7_fb2_pk',data=fbk_data
           IF NOT KEYWORD_SET(fbk_data) THEN BEGIN
              timespan,YYYY+'-'+MM+'-'+DD+'/'+hour+':'+min+':'+second,3,/h
              rbsp_load_efw_fbk,probe='a',type='calibrated'             
           ENDIF
           get_data,'rbspa_efw_fbk_7_fb2_pk_5',data=data
           IF NOT KEYWORD_SET(data) THEN split_vec,'rbspa_efw_fbk_7_fb2_pk' 
           tclip_start=YYYY+'-'+MM+'-'+DD+'/'+hour+':'+min+':'+second
           tclip_end_min=FIX(min)+span/60
           tclip_end_hour=hour
           tclip_end_day=DD
           IF tclip_end_min GE 60 THEN BEGIN
              tclip_end_min=STRING(tclip_end_min-60,FORMAT='(I2.2)')
              tclip_end_hour=STRING(FIX(tclip_end_hour)+1,FORMAT='(I2.2)')
              IF tclip_end_hour EQ '24' THEN BEGIN
                 tclip_end_hour='00'
                 tclip_end_day=STRING(FIX(DD)+1,FORMAT='(I2.2)')
              ENDIF              
           ENDIF
           tclip_end_min=STRING(tclip_end_min,FORMAT='(I2.2)')
           
           tclip_end=YYYY+'-'+MM+'-'+tclip_end_day+'/'+tclip_end_hour+':'+tclip_end_min+':'+second
           time_clip,'rbspa_efw_fbk_7_fb2_pk_5',tclip_start,tclip_end
           get_data,'rbspa_efw_fbk_7_fb2_pk_5_tclip',data=fbk_data
           
              fbk_tim_changed=fbk_data.x;STRARR(N_ELEMENTS(fbk_tim)-1)
              fbk_data_changed=fbk_data.y;DBLARR(N_ELEMENTS(fbk_tim)-1)
              timespan,YYYY+'-'+MM+'-'+DD+'/'+hour+':'+min+':'+second,dura,/hour

; cut data
              start_time = YYYY+'-'+MM+'-'+DD+'/'+hour+':'+min+':'+second
              s_range=WHERE(img_changed_tim EQ start_time)
              s_range=s_range(0)
              
              img_changed_cut=img_changed(s_range : s_range+num_img-1,*,*)
              img_changed_tim_ms_cut=img_changed_tim_ms(s_range : s_range+num_img-1)

              IF N_ELEMENTS(fbk_tim_changed) LT N_ELEMENTS(img_changed_tim_ms_cut) THEN BEGIN
                 img_changed_tim_ms_cut=img_changed_tim_ms_cut(0:N_ELEMENTS(fbk_tim_changed)-1)
                 img_changed_cut=img_changed_cut(0:N_ELEMENTS(fbk_tim_changed)-1,*,*)
              ENDIF ELSE BEGIN
                 fbk_data_changed=fbk_data_changed(0 : num_img-1)
                 fbk_tim_changed=fbk_tim_changed(0 : num_img-1)
              ENDELSE 

              print,'--------------------------------------------------------------------------------------'
              print,'Start time(fbk): '+time_string(fbk_tim_changed(0))
              print,'End time(fbk): '+time_string(fbk_tim_changed(N_ELEMENTS(fbk_tim_changed)-1))

              print,'Start time(img): '+time_string(img_changed_tim_ms_cut(0))
              print,'End time(img): '+time_string(img_changed_tim_ms_cut(N_ELEMENTS(img_changed_tim_ms_cut)-1))
              print,'--------------------------------------------------------------------------------------'
;================================================================================================================
           
; change time resolution to 1s about chorus data
           ;i=0L
           ;FOR k=0,N_ELEMENTS(fbk_tim)-1 DO BEGIN
           ;   IF (k MOD tr_sat) EQ 0 THEN BEGIN
           ;      FOR s=0,tr_sat-1 DO BEGIN
           ;         fbk_data_changed(i)+=fbk_data(k+s)
           ;      ENDFOR

           ;      fbk_tim_changed(i)=fbk_tim(k)
           ;      i++
           ;   ENDIF
           ;ENDFOR
           ;fbk_data_changed=fbk_data_changed/tr_sat
; cut fbk data
           ;tmp1=DD & tmp2=hour & tmp3=min & tmp4=second
           
           ;start_time = YYYY+'-'+MM+'-'+DD+'/'+hour+':'+min+':'+second
           ;min=STRING(FIX(min)+span/60,FORMAT='(I2.2)')
           ;IF second EQ '60' THEN BEGIN
           ;   min = STRING(FIX(min)+1,FORMAT='(I2.2)')
           ;   second='00'
           ;ENDIF           
           
           ;IF min EQ '60' THEN BEGIN
           ;   hour = STRING(FIX(hour)+1,FORMAT='(I2.2)')
           ;   min='00'
           ;ENDIF
           
           ;IF hour EQ '24' THEN BEGIN
           ;   DD = STRING(FIX(DD)+1,FORMAT='(I2.2)')
           ;   hour='00'
           ;   min='00'
           ;ENDIF
           
          ;end_time = YYYY+'-'+MM+'-'+DD+'/'+hour+':'+min+':'+second

          ;s_range=where(fbk_tim_changed EQ start_time)
          ;s_range=s_range(0)
          ;e_range=where(fbk_tim_changed EQ end_time)
          ;e_range=e_range(N_ELEMENTS(e_range)-1)
          ;fbk_tim_changed=fbk_tim_changed(s_range : e_range-1)
          ;fbk_data_changed=fbk_data_changed(s_range : e_range-1)

          ;DD=tmp1 & hour=tmp2 & min=tmp3 & second=tmp4

          ;fbk_data_changed = normalize_routine(fbk_data_changed)

          ;PRINT,'Complete loading fbk data'
           
; change time resolution to 1s about img data
           ;i=0L
           ;FOR k=0,num_raw-1 DO BEGIN

           ;   IF (k MOD tr_img) EQ 0 THEN BEGIN
                 
           ;      FOR j=0,tr_img-1 DO BEGIN
           ;         IF j+k GE num_raw THEN BEGIN
           ;            j=j-1
           ;            BREAK
           ;        ENDIF
                    
           ;         img_changed(i,*,*)+=img_raw(j+k,*,*)
           ;      ENDFOR

           ;      img_changed(i,*,*)=img_changed(i,*,*)/(j+1)
                 
           ;         img_changed_tim_ms(i)=STRING(yr_raw(k),FORMAT='(I4.4)')+'-'+STRING(mo_raw(k),FORMAT='(I2.2)')+'-'+STRING(dy_raw(k),FORMAT='(I2.2)')+'/'+STRING(hh_raw(k),FORMAT='(I2.2)')+':'+STRING(mm_raw(k),FORMAT='(I2.2)')+':'+STRING(ss_raw(k),FORMAT='(I2.2)') +':'+STRING(ms_raw(k),FORMAT='(I3.3)')

           ;         IF ms_raw(k) GE 998 THEN BEGIN
           ;            img_changed_tim(i)=STRING(yr_raw(k),FORMAT='(I4.4)')+'-'+STRING(mo_raw(k),FORMAT='(I2.2)')+'-'+STRING(dy_raw(k),FORMAT='(I2.2)')+'/'+STRING(hh_raw(k),FORMAT='(I2.2)')+':'+STRING(mm_raw(k),FORMAT='(I2.2)')+':'+STRING(ss_raw(k)+1,FORMAT='(I2.2)')
           ;         ENDIF ELSE BEGIN
           ;            img_changed_tim(i)=STRING(yr_raw(k),FORMAT='(I4.4)')+'-'+STRING(mo_raw(k),FORMAT='(I2.2)')+'-'+STRING(dy_raw(k),FORMAT='(I2.2)')+'/'+STRING(hh_raw(k),FORMAT='(I2.2)')+':'+STRING(mm_raw(k),FORMAT='(I2.2)')+':'+STRING(ss_raw(k),FORMAT='(I2.2)')
           ;         ENDELSE
                                     
           ;      PRINT,'changing time resolution of img:'+STRING(i)+'/'+STRING(FIX(span+120-1),FORMAT='(I3.3)');+STRING(span-1,FORMAT='(I3.3)')
                 
           ;      i++
           ;   ENDIF
              
           ;ENDFOR

           ;smo_data=smooth(img_changed,120)

; Cut img data           
          ; tmp1=DD & tmp2=hour & tmp3=min & tmp4=second
           
          ; start_time = YYYY+'-'+MM+'-'+DD+'/'+hour+':'+min+':'+second
          ; min=STRING(FIX(min)+span/60,FORMAT='(I2.2)')
          ; IF second EQ '60' THEN BEGIN
          ;    min = STRING(FIX(min)+1,FORMAT='(I2.2)')
          ;    second='00'
          ; ENDIF           
           
          ; IF min EQ '60' THEN BEGIN
          ;    hour = STRING(FIX(hour)+1,FORMAT='(I2.2)')
          ;    min='00'
          ; ENDIF
           
          ; IF hour EQ '24' THEN BEGIN
          ;    DD = STRING(FIX(DD)+1,FORMAT='(I2.2)')
          ;    hour='00'
          ;    min='00'
          ; ENDIF
           
          ;end_time = YYYY+'-'+MM+'-'+DD+'/'+hour+':'+min+':'+second

          ;s_range=where(img_changed_tim EQ start_time)
          ;s_range=s_range(0)
          ;e_range=where(img_changed_tim EQ end_time)
          ;e_range=e_range(N_ELEMENTS(e_range)-1)
          ;img_changed_tim=img_changed_tim(s_range : e_range-1)
          ;img_changed_tim_ms=img_changed_tim_ms(s_range : e_range-1)
          ;img_changed=img_changed(s_range : e_range-1,*,*)
          ;smo_data=smo_data(s_range : e_range-1,*,*)

          ;DD=tmp1 & hour=tmp2 & min=tmp3 & second=tmp4
           
           
           ;img_changed_tim=time_double(img_changed_tim)

           ;IF WHERE(img_changed EQ 0) NE -1 THEN STOP
          
           ;PRINT,'Complete loading img data'

           ;IF N_ELEMENTS(img_changed(*,0,0)) LT N_ELEMENTS(fbk_data_changed) THEN BEGIN
           ;   fbk_data_changed=fbk_data_changed(0:N_ELEMENTS(img_changed(*,0,0))-1)
           ;   fbk_tim_changed=fbk_tim_changed(0:N_ELEMENTS(img_changed(*,0,0))-1)
           ;ENDIF

           ;IF N_ELEMENTS(img_changed(*,0,0)) GT N_ELEMENTS(fbk_data_changed) THEN BEGIN
           ;   img_changed=img_changed(0:N_ELEMENTS(fbk_data_changed)-1,*,*)
           ;   smo_data=smo_data(0:N_ELEMENTS(fbk_data_changed)-1,*,*)
           ;   img_changed_tim=img_changed_tim(0:N_ELEMENTS(fbk_data_changed)-1)
           ;   img_changed_tim_ms=img_changed_tim_ms(0:N_ELEMENTS(fbk_data_changed)-1)
           ;ENDIF
              
; ============================================================================================================

              
; calculate correlation coefficient about each pixel
           j=0L
           FOR x=0,255 DO BEGIN
              
                 PRINT,'calculating correlate coefficient about each pixel:'+STRING(x)+'/255'
                 
                 FOR y=0,255 DO BEGIN

                    IF (x LT xlim_pxl) AND (y LT ylim_pxl) THEN BEGIN 
                    
                       img_tmp = img_changed_cut(*,x,y)-smooth(img_changed_cut(*,x,y),480,/EDGE_MIRROR)

                       ;img_tmp=NORMALIZE_routine(img_tmp)
                         
                       cor_tmp=C_CORRELATE(img_tmp,fbk_data_changed,lag)
                       
                       cor_tmp=MAX(cor_tmp)
                    
                       cor(x,y)=cor_tmp
                    ENDIF ELSE BEGIN

                       cor(x,y)=0.0

                    ENDELSE
                                                            
                ; search max pixel
                    IF (cor(x,y) GT max_pixel_cor) AND (ele_psa(x,y) GT 15.0)  THEN BEGIN
                       max_pixel(0)=x
                       max_pixel(1)=y
                       max_pixel_cor=cor(x,y)
                    ENDIF
                                      
                 ENDFOR
              ENDFOR           

           RETURN,cor
        END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	 PLOT_COR_CONTOUR
;
        PRO plot_cor_contour,day,hhmmss,dura,span,no_plot=no_plot

           COMMON prf_psa
           COMMON fov_psa
           COMMON raw_psa
           COMMON cor_psa
           COMMON struct_psa
           
; input time info
           IF NOT KEYWORD_SET(day) THEN BEGIN
              day=''
              hhmmss=''
              
              READ,day,PROMPT='Enter day(YYYYMMDD):'
              READ,hhmmss,PROMPT='Enter start time(HHMMSS):'
              READ,dura,PROMPT='Enter duration(hour):'
              READ,span,PROMPT='Enter calculation span of CC(second):'

           ENDIF
           
              YYYY=STRMID(day,0,4)
              MM=STRMID(day,4,2)
              DD=STRMID(day,6,2)
              hhmm=STRMID(hhmmss,0,4)
              hour=STRMID(hhmm,0,2)
              min=STRMID(hhmm,2,2)
              second=STRMID(hhmmss,4,2)

              n=0L
              max_cor=1.0
              min_cor=0

           timespan,YYYY+'-'+MM+'-'+DD+'/'+hour+':'+min+':'+second,dura,/hour

           ele_lim=10.0

           IF NOT KEYWORD_SET(num_raw) THEN file_psa_raw,2,day,hhmm,span/60+2
           
           
; calculate correlation coefficient about each pixel
            IF NOT KEYWORD_SET(no_plot) THEN cor = ROLLING_THUNDER(day,hhmmss,dura,span)

; over plot contour on the map about each time span

           ; plot map
              IF KEYWORD_SET(no_plot) THEN BEGIN
                 PLOT_GEOGRA_MAP,day,hhmmss,dura,10E5*6.8,'0',no_bar=1,no_change=1,no_plot=1
              ENDIF ELSE BEGIN
                 PLOT_GEOGRA_MAP,day,hhmmss,dura,10E5*6.8,'0',no_bar=1,no_change=1
              ENDELSE
              
           ; Prepare array for contouring
              IF NOT KEYWORD_SET(ele_lim) THEN ele_lim=10.0
              tmp_bin=WHERE(gla_psa NE 99999.9 AND gla_psa NE 99999.9 AND ele_psa GE ele_lim,no_tmp_bin)
              x_cont=MAKE_ARRAY(no_tmp_bin,/FLOAT,VALUE=0.0)
              y_cont=MAKE_ARRAY(no_tmp_bin,/FLOAT,VALUE=0.0)
              z_cont=MAKE_ARRAY(no_tmp_bin,/FLOAT,VALUE=0.0)

              num_cnt=0L
              
           ; Put values into the arrays for contouring
              FOR x=0,255 DO BEGIN
                 FOR y=0,255 DO BEGIN
                    IF gla_psa(x,y) NE 99999.9 AND glo_psa(x,y) NE 99999.9 AND ele_psa(x,y) GE ele_lim THEN BEGIN

                       z_cont(num_cnt)=FLOAT(cor(x,y)) ;*COS(zan_psa(x,y)*!DTOR)
                       IF z_cont(num_cnt) GT max_cor THEN z_cont(num_cnt)=max_cor
                       IF z_cont(num_cnt) LT min_cor THEN z_cont(num_cnt)=min_cor
                    
                       x_cont(num_cnt)=glo_psa(x,y)
                       y_cont(num_cnt)=gla_psa(x,y)

                       num_cnt++
                    ENDIF
                 ENDFOR
              ENDFOR
              
            ; Contouring
              IF NOT KEYWORD_SET(no_plot) THEN BEGIN
              loadct,[39]
              contour_lines=30
              levels_set=(max_cor-min_cor)*FINDGEN(contour_lines)/FLOAT(contour_lines)+min_cor
              colour_set=!D.TABLE_SIZE*INDGEN(contour_lines)/contour_lines
              CONTOUR,z_cont,x_cont,y_cont,/OVERPLOT,/CELL_FILL,LEVELS=levels_set,C_COLORS=colour_set,/IRREGULAR,MIN_VALUE=0.25

           
; Plot colour bar
          IF KEYWORD_SET(no_bar) THEN  DELVAR,no_bar
           
           position=[!P.POSITION(2)+0.01,!P.POSITION(1),!P.POSITION(2)+0.02,!P.POSITION(3)]

           plot_colour_bar,position=position,legend=legend,bar_title='cor',max_val=max_cor,min_val=min_cor
           
; Plot Max pixel and Max correlation coefficient
           PSYM8,/fill
           print,glo_psa(max_pixel(0),max_pixel(1))
           print,gla_psa(max_pixel(0),max_pixel(1))
           print,max_pixel
           OPLOT,[glo_psa(max_pixel(0),max_pixel(1))],[gla_psa(max_pixel(0),max_pixel(1))],PSYM=8,SYMSIZE=SYMSIZE,COL=fsc_color("black")

           XYOUTS,0.4,0.01,'Max Correlate Coefficient:'+STRING(max_pixel_cor),CHARSIZE=2.0,/NORMAL
           
; Save as png
           IF FINITE(MAX(cor),/NAN) NE 1 THEN BEGIN
              makepng,'/home/suguru/spedas/image_data/cor_img/'+day+'-'+hhmmss+'-'+STRING(span/60,FORMAT='(I1.1)')+'min'
           ENDIF

           ENDIF
           cor_contour=CREATE_STRUCT(NAME='cor_contour','x',x_cont,'y',y_cont,'z',z_cont,'mes',max_pixel_cor,'pos',!P.POSITION)
           END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	 PLOT_COR_LINE
;
        PRO plot_cor_line,day,hhmmss,dura,span,merged=merged

           COMMON prf_psa
           COMMON fov_psa
           COMMON raw_psa
           COMMON cor_psa

; input time info
           IF NOT KEYWORD_SET(day) THEN BEGIN
              day=''
              hhmmss=''
              
              READ,day,PROMPT='Enter day(YYYYMMDD):'
              READ,hhmmss,PROMPT='Enter start time(HHMMSS):'
              READ,dura,PROMPT='Enter duration(hour):'
              READ,span,PROMPT='Enter calculation span of CC(second):'

           ENDIF
              
              YYYY=STRMID(day,0,4)
              MM=STRMID(day,4,2)
              DD=STRMID(day,6,2)
              hhmm=STRMID(hhmmss,0,4)
              hour=STRMID(hhmm,0,2)
              min=STRMID(hhmm,2,2)
              second=STRMID(hhmmss,4,2)
              
              IF NOT KEYWORD_SET(img_raw) THEN file_psa_raw,2,day,hhmm,span/60+2
              
              ;FOR i=0,(dura*60)/(cal_span/60)-1 DO BEGIN
                 
                 ;IF FIX(STRMID(hhmmss,2,2))-FIX(min) EQ 1 THEN BEGIN
                 ;   FREE_VARIABLE
                 ;   file_psa_raw,2,day,hhmm,span/60+1
                 ;ENDIF
              
              PRINT,'Current time: '+day+'-'+hhmmss

                 
; Calculate Correlate Coefficient                 
                 cor=rolling_thunder(day,hhmmss,dura,span)
                 timespan,YYYY+'-'+MM+'-'+DD+'/'+STRMID(hhmmss,0,2)+':'+STRMID(hhmmss,2,2)+':'+STRMID(hhmmss,4,2),span/60,/min
; Normalize data
                 img_normed=normalize_routine(img_changed(*,max_pixel(0),max_pixel(1))-smooth(img_changed(*,max_pixel(0),max_pixel(1)),480,/EDGE_MIRROR))
                 fbk_normed=normalize_routine(fbk_data_changed-smooth(fbk_data_changed,480,/EDGE_MIRROR))

; Store img data and fbk data to tplot value
                 store_data,'normed_fbk_data_changed',data={x:fbk_tim_changed,y:fbk_normed}
                 store_data,'normed_img_data_changed',data={x:img_changed_tim_ms,y:img_normed}
                 ;ylim,'normed_fbk_data_changed',0,1,1
                 ;ylim,'normed_img_data_changed',0,1,0
; Merge img plot and fbk plot
                 store_data,'merged_plot',data=['normed_fbk_data_changed','normed_img_data_changed']

; Set plot option and plot
                 options,'merged_plot','ytitle','Normalized Intensity[Arbitrary Unit]'
                 options,'merged_plot','xtitle',''
                 loadct,39
                 IF KEYWORD_SET(merged) THEN tplot_options,'region',[0.0,0.05,1.0,0.45] ;[0.5,0.6,0.95,0.99]
                 thm_init,/reset
                 options,'merged_plot','colors',[6,4]
                 tplot,'merged_plot'
              
                 XYOUTS,0.4,0.01,'Correlate Coefficient:'+STRING(max_pixel_cor),CHARSIZE=2.0,/NORMAL

; Save img as png
                 ;makepng,'/home/suguru/spedas/image_data/merged_cor_line/'+day+'-'+hhmmss+'-'+STRING(span/60,FORMAT='(I1.1)')+'min(line)'



                 
; Plot fbk data
              ;options,'normed_fbk_data_changed','ytitle','Chorus Intensity [Arbitrary Unit]'
              ;options,'normed_fbk_data_changed','xtitle',''
              ;tplot,'normed_fbk_data_changed'

; Save fbk line plot as png
              ;makepng,'line_fbk_raw/'+day+'-'+hhmm+'-'+STRING(span/60,FORMAT='(I1.1)')+'min(fbk plot)'

; Plot img data at max pixel
              ;options,'normed_img_data_changed','ytitle','Auroral Intensity [Arbitrary Unit]'
              ;options,'normed_img_data_changed','xtitle',''
              ;tplot,'normed_img_data_changed'

; Save img line plot as png
              ;makepng,'line_img_raw/'+day+'-'+hhmm+'-'+STRING(span/60,FORMAT='(I1.1)')+'min(img
              ;plot)'


                 
              
; Set next time
                 ;min=STRMID(hhmmss,2,2)
                 ;hhmmss = STRING(LONG(hhmmss)+cal_span,FORMAT='(I6.6)')
                 ;IF STRMID(hhmmss,4,2) EQ '60' THEN BEGIN
                 ;   hhmmss = STRING(LONG(hhmmss)+100-60,FORMAT='(I6.6)')                
                 ;ENDIF
              
              
                 ;IF STRMID(hhmmss,2,2) EQ '60' THEN BEGIN
                 ;   hhmmss = STRING(LONG(hhmmss)+10000-6000,FORMAT='(I6.6)')
                 ;ENDIF

                 ;IF STRMID(hhmmss,0,2) EQ '24' THEN BEGIN
                 ;   day = STRING(LONG(day)+1,FORMAT='(I8.8)')
                 ;   hhmmss = '000000'
                 ;ENDIF
                 ;hhmm=STRMID(hhmmss,0,4)

              ;ENDFOR
        END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	 PLOT_FFT_CONTOUR
;
        PRO plot_fft_contour,day,hhmmss,dura,span,no_plot=no_plot

           COMMON prf_psa
           COMMON fov_psa
           COMMON raw_psa
           COMMON cor_psa
           COMMON struct_psa
           
           IF NOT KEYWORD_SET(day) THEN BEGIN
              day=''
              hhmmss=''
              
              READ,day,PROMPT='Enter day(YYYYMMDD):'
              READ,hhmmss,PROMPT='Enter start time(HHMMSS):'
              READ,dura,PROMPT='Enter duration(hour):'
              READ,span,PROMPT='Enter calculation span of FFT(second):'

           ENDIF

                      
              YYYY=STRMID(day,0,4)
              MM=STRMID(day,4,2)
              DD=STRMID(day,6,2)
              hhmm=STRMID(hhmmss,0,4)
              hour=STRMID(hhmm,0,2)
              min=STRMID(hhmm,2,2)
              second=STRMID(hhmmss,4,2)

              n=0L

              timespan,YYYY+'-'+MM+'-'+DD+'/'+hour+':'+min+':'+second,dura,/hour

              ele_lim=10.0
           
              IF NOT KEYWORD_SET(num_raw) THEN file_psa_raw,2,day,hhmm,span/60+2

; Preparing arrays
           IF NOT KEYWORD_SET(img_changed) THEN img_changed=CHANGE_TIME_RES_IMG(img_raw);DBLARR(span+120,256,256)
           ;img_changed_tim=STRARR(span+120)
           ;img_changed_tim_ms=STRARR(span+120)
           fft_peak=FLTARR(256,256)
           tr_change=8.0
           num_img=N_ELEMENTS(img_raw(*,0,0))*(tr_change/tr_img)
           max_peak=40.0
           min_peak=0.0
           smo_data=smooth(img_changed,120)
              
; Cut img data           
           start_time = YYYY+'-'+MM+'-'+DD+'/'+hour+':'+min+':'+second
           s_range=WHERE(img_changed_tim EQ start_time)
           s_range=s_range(0)
              
           img_changed_cut=img_changed(s_range : s_range+num_img-1,*,*)
           img_changed_tim_ms_cut=img_changed_tim_ms(s_range : s_range+num_img-1)
        
;;;;;;;;;;;;;;;;;;;;;;;;;; FFT about each pixel ;;;;;;;;;;;;;;;;;;;;;;;;;;;

           FOR x=0,255 DO BEGIN
              PRINT,'calculating FFT about each pixel:'+STRING(x)+'/255'
                 FOR y=0,255 DO BEGIN
                 
                 ;Set time range for FFT
                 
                    num_data=N_ELEMENTS(img_changed_cut(*,0,0))
                    ts=0.125
                    total_time=ts*(num_data-1)
                    fft_results=DBLARR(num_data/2)
                    period=INDGEN(num_data/2)*ts
                    vib_num=INDGEN(num_data/2)
                    tr_data=img_changed_cut(*,x,y)
                    tr_data=tr_data-normalize_routine(smo_data(*,x,y))
                    
                    fft_results=ABS(FFT(tr_data))
                    fft_results=fft_results[0:num_data/2-1]
                    fft_data=REVERSE(fft_results)

                    FOR i=0,(num_data/2)-1 DO BEGIN
                       period(i)=total_time/vib_num(i)
                    ENDFOR

                    period[0]=period(1)*1.2
                    period=REVERSE(period)

                    tmp=WHERE(period GT 30)
                    IF N_ELEMENTS(tmp) GT 1 THEN tmp=tmp(0)

                    period=period(0:tmp-1)
                    fft_data=fft_data(0:tmp-1)
                    tmp=MAX(fft_data,max_sub)
                    fft_peak(x,y)=period(max_sub)
                                ;fft_data=smooth(fft_data,2)
                                ;store_data,'test',data={x:time_double(img_changed_tim_ms),y:tr_data}
                                ;timespan,YYYY+'-'+MM+'-'+DD+'/'+hour+':'+min+':'+second,2,/min
                                ;tplot,'test'
                                ;PLOT,period,fft_data
                    
                 ENDFOR
              ENDFOR
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              
; plot map           
              PLOT_GEOGRA_MAP,day,hhmmss,dura,10E5*6.8,'0',no_bar=1,no_change=1
              
; Prepare array for contouring
              IF NOT KEYWORD_SET(ele_lim) THEN ele_lim=10.0
              tmp_bin=WHERE(gla_psa NE 99999.9 AND gla_psa NE 99999.9 AND ele_psa GE ele_lim,no_tmp_bin)
              x_cont=MAKE_ARRAY(no_tmp_bin,/FLOAT,VALUE=0.0)
              y_cont=MAKE_ARRAY(no_tmp_bin,/FLOAT,VALUE=0.0)
              z_cont=MAKE_ARRAY(no_tmp_bin,/FLOAT,VALUE=0.0)

              num_cnt=0L
              
; Put values into the arrays for contouring

              FOR x=0,255 DO BEGIN
                 FOR y=0,255 DO BEGIN
                    IF gla_psa(x,y) NE 99999.9 AND glo_psa(x,y) NE 99999.9 AND ele_psa(x,y) GE ele_lim THEN BEGIN

                       z_cont(num_cnt)=FLOAT(fft_peak(x,y)) ;*COS(zan_psa(x,y)*!DTOR)
                       IF z_cont(num_cnt) GT max_peak THEN z_cont(num_cnt)=max_peak
                       IF z_cont(num_cnt) LT min_peak THEN z_cont(num_cnt)=min_peak
                    
                       x_cont(num_cnt)=glo_psa(x,y)
                       y_cont(num_cnt)=gla_psa(x,y)

                       num_cnt++
                    ENDIF
                 ENDFOR
              ENDFOR
              
; Contouring
              IF NOT KEYWORD_SET(no_plot) THEN BEGIN
              loadct,[39]
              contour_lines=30
              levels_set=(max_peak-min_peak)*FINDGEN(contour_lines)/FLOAT(contour_lines)+min_peak
              colour_set=!D.TABLE_SIZE*INDGEN(contour_lines)/contour_lines
              CONTOUR,z_cont,x_cont,y_cont,/OVERPLOT,/CELL_FILL,LEVELS=levels_set,C_COLORS=colour_set,/IRREGULAR
        
; overlay satelite orbit
              timespan,YYYY_orb+'-'+MM_orb+'-'+DD_orb+'/'+hour_orb+':'+min_orb+':00',dura_orb,/hour
              overlay_map_sc_ifoot,'rbspa_ifoot_geo_lat','rbspa_ifoot_geo_lon'
              timespan,YYYY+'-'+MM+'-'+DD+'/'+hour+':'+min+':00',dura,/hour
              PLOT_POS_FOOTPRINT
              
; Plot colour bar
             IF KEYWORD_SET(no_bar) THEN  DELVAR,no_bar
           
             position=[!P.POSITION(2)+0.01,!P.POSITION(1),!P.POSITION(2)+0.02,!P.POSITION(3)]

             plot_colour_bar,position=position,legend=legend,bar_title='fft',max_val=max_peak,min_val=min_peak

             ENDIF
              
             fft_contour=CREATE_STRUCT(NAME='fft_contour','x',x_cont,'y',y_cont,'z',z_cont,'pos',!P.POSITION)
              
              END
        
                 
;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	 SAVE_PNG_FFT
;
        PRO save_png_fft,day,hhmmss,dura,span

           COMMON prf_psa
           COMMON fov_psa
           COMMON raw_psa
           COMMON cor_psa

           IF NOT KEYWORD_SET(day) THEN BEGIN
              day=''
              hhmmss=''
              
              READ,day,PROMPT='Enter day(YYYYMMDD):'
              READ,hhmmss,PROMPT='Enter start time(HHMMSS):'
              READ,dura,PROMPT='Enter duration(hour):'
              READ,span,PROMPT='Enter calculation span of FFT(second):'

           ENDIF

                      
           YYYY=STRMID(day,0,4)
           MM=STRMID(day,4,2)
           DD=STRMID(day,6,2)
           hhmm=STRMID(hhmmss,0,4)
           hour=STRMID(hhmm,0,2)
           min=STRMID(hhmm,2,2)
           second=STRMID(hhmmss,4,2)
           
           n=0L
           
           timespan,YYYY+'-'+MM+'-'+DD+'/'+hour+':'+min+':'+second,dura,/hour
           
           ele_lim=10.0

           file_psa_raw,2,day,hhmm,span/60+2
           FOR i=0,(dura*60)/(cal_span/60)-1 DO BEGIN

              IF FIX(STRMID(hhmmss,2,2))-FIX(min) EQ 1 THEN BEGIN
                 FREE_VARIABLE
                 file_psa_raw,2,day,hhmm,span/60+2
              ENDIF

              PRINT,'Current time: '+day+'-'+hhmmss

; Plot FFT peak with contour
              plot_fft_contour,day,hhmmss,dura,span

; Save image as png
              makepng,'/home/suguru/spedas/image_data/fft_peak_contour/'+day+'-'+hhmmss+'-'+STRING(span/60,FORMAT='(I1.1)')+'min'
              
; Set next time
           min=STRMID(hhmmss,2,2)
           hhmmss = STRING(LONG(hhmmss)+cal_span,FORMAT='(I6.6)')
           IF STRMID(hhmmss,4,2) EQ '60' THEN BEGIN
              hhmmss = STRING(LONG(hhmmss)+100-60,FORMAT='(I6.6)')                
           ENDIF
              
              
           IF STRMID(hhmmss,2,2) EQ '60' THEN BEGIN
              hhmmss = STRING(LONG(hhmmss)+10000-6000,FORMAT='(I6.6)')
           ENDIF

           IF STRMID(hhmmss,0,2) EQ '24' THEN BEGIN
              day = STRING(LONG(day)+1,FORMAT='(I8.8)')
              hhmmss = '000000'
           ENDIF
           hhmm=STRMID(hhmmss,0,4)
           
        ENDFOR

        END
;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	 PLOT_FFT_LINE
;
        PRO plot_fft_line,day,hhmmss,span,no_img_line=no_img_line,fft_hist=fft_hist,count=count

           COMMON prf_psa
           COMMON fov_psa
           COMMON raw_psa
           COMMON cor_psa

; Set time info
           IF NOT KEYWORD_SET(day) THEN BEGIN
              day=''
              hhmmss=''
              
              READ,day,PROMPT='Enter day(YYYYMMDD):'
              READ,hhmmss,PROMPT='Enter start time(HHMMSS):'
              READ,span,PROMPT='Enter calculation span of FFT(second):'

           ENDIF

                      
              YYYY=STRMID(day,0,4)
              MM=STRMID(day,4,2)
              DD=STRMID(day,6,2)
              hhmm=STRMID(hhmmss,0,4)
              hour=STRMID(hhmm,0,2)
              min=STRMID(hhmm,2,2)
              second=STRMID(hhmmss,4,2)

              n=0L

              timespan,YYYY+'-'+MM+'-'+DD+'/'+hour+':'+min+':'+second,dura,/hour

              ele_lim=10.0
              tr_change=8.0

              IF NOT KEYWORD_SET(img_raw) THEN file_psa_raw,2,day,hhmm,span/60+2
              num_img=N_ELEMENTS(img_raw(*,0,0))*(tr_change/tr_img)-2*6000*(tr_change/tr_img)
              
              
; Get highest correlate coeficient pixel
              
              IF NOT KEYWORD_SET(cor) THEN cor=ROLLING_THUNDER(day,hhmmss,2,span)

              IF NOT KEYWORD_SET(img_changed) THEN img_changed=CHANGE_TIME_RES_IMG(img_raw)

; Cut img data           
              start_time = YYYY+'-'+MM+'-'+DD+'/'+hour+':'+min+':'+second
              s_range=WHERE(img_changed_tim EQ start_time)
              s_range=s_range(0)
              
              img_changed_cut=img_changed(s_range : s_range+num_img-1,*,*)
              img_changed_tim_ms_cut=img_changed_tim_ms(s_range : s_range+num_img-1)
              img_changed_cut=normalize_routine(img_changed_cut(*,max_pixel(0),max_pixel(1)))
              
              
;;;;;;;;;;;;;;;;;;;;;;;;;; FFT about img ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                 
                 ;Set time range for FFT
                 
                    num_data=N_ELEMENTS(img_changed_cut)
                    ts=0.125
                    total_time=ts*(num_data-1)
                    fft_results=DBLARR(num_data/2)
                    period_img=INDGEN(num_data/2)*ts
                    vib_num=INDGEN(num_data/2)
                    smo_data=smooth(img_changed_cut,240,/EDGE_MIRROR)
                    tr_data=img_changed_cut-smo_data
                    
                    fft_results=ABS(FFT(tr_data))
                    fft_results=fft_results[0:num_data/2-1]
                    fft_data_img=REVERSE(fft_results)

                    FOR i=0,(num_data/2)-1 DO BEGIN
                       period_img(i)=total_time/vib_num(i)
                    ENDFOR

                    period_img[0]=period_img(1)*1.2
                    period_img=REVERSE(period_img)
                    print,period_img
                    tmp=WHERE(period_img GT 30)
                    IF N_ELEMENTS(tmp) GT 1 THEN tmp=tmp(0)

                    period_img=period_img(0:tmp-1)
                    print,period_img
                    fft_data_img=fft_data_img(0:tmp-1)
                    tmp_max=MAX(fft_data_img,max_sub)
                    max_period_img=period_img(max_sub)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                    IF KEYWORD_SET(no_img_line) OR KEYWORD_SET(fft_hist) THEN BEGIN
                       
;;;;;;;;;;;;;;;;;;;;;;;;;; FFT about fbk ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                 
                 ;Set time range for FFT
                 
                       num_data=N_ELEMENTS(fbk_data_changed)
                       smo_data=smooth(fbk_data_changed,240,/EDGE_MIRROR)
                       ts=0.125
                       total_time=ts*(num_data-1)
                       fft_results=DBLARR(num_data/2)
                       period_fbk=INDGEN(num_data/2)*ts
                       vib_num=INDGEN(num_data/2)
                       tr_data=fbk_data_changed-smo_data
                       
                       fft_results=ABS(FFT(tr_data))
                       fft_results=fft_results[0:num_data/2-1]
                       fft_data_fbk=REVERSE(fft_results)

                       FOR i=0,(num_data/2)-1 DO BEGIN
                          period_fbk(i)=total_time/vib_num(i)
                       ENDFOR

                       period_fbk[0]=period_fbk(1)*1.2
                       period_fbk=REVERSE(period_fbk)

                       tmp=WHERE(period_fbk GT 30)
                       IF N_ELEMENTS(tmp) GT 1 THEN tmp=tmp(0)

                       period_fbk=period_fbk(0:tmp-1)
                       fft_data_fbk=fft_data_fbk(0:tmp-1)
                       tmp_max=MAX(fft_data_fbk,max_sub)
                       max_period_fbk=period_fbk(max_sub)

                    
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                       IF NOT KEYWORD_SET(fft_hist) THEN BEGIN
                          PLOT,period_img,fft_data_img,xtitle='Period [s]',ytitle='Amplitude (IMG)',CHARSIZE=2.0,position=[0.06,0.79,0.5,0.94]
                          PLOT,period_fbk,fft_data_fbk,xtitle='Period [s]',ytitle='Amplitude (FBK)',CHARSIZE=2.0,position=[0.06,0.54,0.5,0.69]
                       ENDIF

                       IF KEYWORD_SET(fft_hist) THEN BEGIN
                          fft_hist_img(count)=max_period_img
                          fft_hist_fbk(count)=max_period_fbk
                       ENDIF
                       
                       ;printf,lun,day,STRING(hhmmss,FORMAT='(I6.6)'),FORMAT='(4(A,TR1))'
                       ;printf,lun,max_period_img,max_period_fbk,FORMAT='(4(F,TR1))'          
                    ENDIF ELSE BEGIN                       
                    
; store data to tplot valuable
                       store_data,'img',data={x:time_double(img_changed_tim_ms_cut),y:img_changed_cut}

; Set plot options
                       timespan,STRMID(day,0,4)+'-'+STRMID(day,4,2)+'-'+STRMID(day,6,2)+'/'+STRMID(hhmmss,0,2)+':'+STRMID(hhmmss,2,2)+':'+STRMID(hhmmss,4,2),span/60,/min

                       options,'img','ytitle','Normalized Intensity [Arbitrary Unit]'
                       options,'img','xtitle','UT'

                       TPLOT_OPTIONS,'region',[0.0,0.5,1.0,1.0]

                       !P.MULTI=[0,1,2]
              
                       TPLOT,'img'
                       PLOT,period_img,fft_data_img,position=[0.045,0.1,0.955,0.5],xtitle='Period [s]',ytitle='Amplitude [Arbitrary Unit]'                  
              
                       FREE_VARIABLE
                    ENDELSE
              END


;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	 PLOT_DYNA_SPEC
;
;               plot dynamic spectrum
;
        PRO plot_dyna_spec,day,hhmmss,dura,span,key=key

           COMMON prf_psa
           COMMON fov_psa
           COMMON raw_psa
           COMMON cor_psa
           
           IF NOT KEYWORD_SET(day) THEN BEGIN
              day=''
              hhmmss=''
              
              READ,day,PROMPT='Enter day(YYYYMMDD):'
              READ,hhmmss,PROMPT='Enter start time(HHMMSS):'
              READ,dura,PROMPT='Enter duration(hour):'
              READ,span,PROMPT='Enter calculation span of FFT(second):'

           ENDIF

                      
              YYYY=STRMID(day,0,4)
              MM=STRMID(day,4,2)
              DD=STRMID(day,6,2)
              hhmm=STRMID(hhmmss,0,4)
              hour=STRMID(hhmm,0,2)
              min=STRMID(hhmm,2,2)
              second=STRMID(hhmmss,4,2)
              s_hhmmss=hhmmss

              n=0L

                 timespan,YYYY+'-'+MM+'-'+DD+'/'+hour+':'+min+':'+second,span/60,/min
                 get_timespan,time
                 tmp=time_string(time)
                 
                 ;file_psa_raw,2,day,hhmm,span/60+2
                 
                 ;img_changed=change_time_res_img(img_raw)
                 ;cor=ROLLING_THUNDER(day,hhmmss,2,span)

                 FOR k=0,60*dura/(span/60)-1 DO BEGIN
                    timespan,STRMID(day,0,4)+'-'+STRMID(day,4,2)+'-'+STRMID(day,6,2)+'/'+STRMID(hhmmss,0,2)+':'+STRMID(hhmmss,2,2)+':'+STRMID(hhmmss,4,2),span/60,/min
                    file_psa_raw,2,day,hhmm,span/60
                    img_changed=CHANGE_TIME_RES_IMG(img_raw)
                    img=img_changed(*,34:94,34:94) ;center (64,64) plus minus 30
                    img_mean=DBLARR(N_ELEMENTS(img(*,0,0)))
                 
                    FOR i=0,N_ELEMENTS(img(*,0,0))-1 DO BEGIN
                       img_mean(i)=MEAN(img(i,*,*))
                    ENDFOR
                    IF k EQ 0 THEN BEGIN
                       store_data,'img',data={x:img_changed_tim_ms,y:img_mean}
                    ENDIF
                    get_data,'img',data=data
                    img_tim=[data.x,img_changed_tim_ms]
                    help,img_tim
                    img_data=[data.y,img_mean]
                    store_data,'img',data={x:img_tim,y:img_data}
                 
                 ; Set next time
                    min=STRMID(hhmmss,2,2)
                    hhmmss = STRING(LONG(hhmmss)+span/60*100,FORMAT='(I6.6)')
                    IF STRMID(hhmmss,4,2) EQ '60' THEN BEGIN
                       hhmmss = STRING(LONG(hhmmss)+100-60,FORMAT='(I6.6)')                
                    ENDIF
              
              
                    IF STRMID(hhmmss,2,2) EQ '60' THEN BEGIN
                       hhmmss = STRING(LONG(hhmmss)+10000-6000,FORMAT='(I6.6)')
                    ENDIF

                    IF STRMID(hhmmss,0,2) EQ '24' THEN BEGIN
                       day = STRING(LONG(day)+1,FORMAT='(I8.8)')
                       hhmmss = '000000'
                    ENDIF
                    hhmm=STRMID(hhmmss,0,4)
                    FREE_VARIABLE
                 ENDFOR  
                       
                 ;store_data,'img',data={x:img_changed_tim_ms,y:img_mean}
                 ;time_clip,'img',tmp(0),tmp(1)
                 tdpwrspc,'img',nbox=2400,nshift=10
                 ;wav_data,'img'
           
              IF KEYWORD_SET(key) THEN BEGIN
                 tplot,'img_dpwrspc',/oplot
              ENDIF ELSE BEGIN
                 timespan,YYYY+'-'+MM+'-'+DD+'/'+hour+':'+min+':'+second,dura,/h
                 tplot,'img_dpwrspc'
              ENDELSE
              
           END
;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	 PLOT_FFT_HIST
;
        PRO plot_fft_hist,day,hhmmss,span,no_img_line=no_img_line,fft_hist=fft_hist,count=count

           COMMON prf_psa
           COMMON fov_psa
           COMMON raw_psa
           COMMON cor_psa

; Set time info
           IF NOT KEYWORD_SET(day) THEN BEGIN
              day=''
              hhmmss=''
              
              READ,day,PROMPT='Enter day(YYYYMMDD):'
              READ,hhmmss,PROMPT='Enter start time(HHMMSS):'
              READ,span,PROMPT='Enter calculation span of FFT(second):'

           ENDIF

                      
              YYYY=STRMID(day,0,4)
              MM=STRMID(day,4,2)
              DD=STRMID(day,6,2)
              hhmm=STRMID(hhmmss,0,4)
              hour=STRMID(hhmm,0,2)
              min=STRMID(hhmm,2,2)
              second=STRMID(hhmmss,4,2)
              fft_hist_img=FLTARR(10000)
              fft_hist_fbk=FLTARR(10000)
              i=0L
              
              WHILE (hhmmss NE '005450') DO BEGIN
                 IF FIX(STRMID(hhmmss,2,2))-FIX(min) GE 1 THEN BEGIN
                    FREE_VARIABLE
                    file_psa_raw,2,day,hhmm,span/60+2
                 ENDIF
                 
                 PLOT_FFT_LINE,day,hhmmss,span,fft_hist=1,count=i
; Set next time
                 min=STRMID(hhmmss,2,2)
                 hhmmss = STRING(LONG(hhmmss)+cal_span,FORMAT='(I6.6)')
                 IF STRMID(hhmmss,4,2) EQ '60' THEN BEGIN
                    hhmmss = STRING(LONG(hhmmss)+100-60,FORMAT='(I6.6)')                
                 ENDIF
              
                 
                 IF STRMID(hhmmss,2,2) EQ '60' THEN BEGIN
                    hhmmss = STRING(LONG(hhmmss)+10000-6000,FORMAT='(I6.6)')
                 ENDIF

                 IF STRMID(hhmmss,0,2) EQ '24' THEN BEGIN
                    day = STRING(LONG(day)+1,FORMAT='(I8.8)')
                    hhmmss = '000000'
                 ENDIF
                 hhmm=STRMID(hhmmss,0,4)
                 i++
              ENDWHILE

              sub_img_tmp=WHERE(fft_hist_img EQ 0)
              sub_fbk_tmp=WHERE(fft_hist_fbk EQ 0)
              fft_hist_img=fft_hist_img(0:sub_img_tmp(0))
              fft_hist_fbk=fft_hist_fbk(0:sub_fbk_tmp(0))

              hist_img=HISTOGRAM(fft_hist_img,BINSIZE=1)
              hist_fbk=HISTOGRAM(fft_hist_fbk,BINSIZE=1)
              xbin=INDGEN(31)
              xbin=xbin(1:30)

              PLOT,xbin,hist_img,PSYM=10,TITLE='Image',ytitle='Occurence',xtitle='Period [s]'
              makepng,'/home/suguru/spedas/image_data/histgram/hist_img'
              PLOT,xbin,hist_fbk,PSYM=10,TITLE='Filter Bank Data',ytitle='Occurence',xtitle='Period [s]'
              makepng,'/home/suguru/spedas/image_data/histgram/hist_fbk'
              
              
           END
        
;--------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	CHANGE_TIME_RES_IMG
;
; PURPOSE:
;
;	change time resolution of img 
;

        FUNCTION change_time_res_img,arr_img

          COMMON prf_psa
          COMMON raw_psa
          COMMON cor_psa

; time resolution you want to change
          tr_change=8.0

; preparing arrays and valiables
          ave_img_val=FLOAT(FIX(100/tr_change))     
          tmp=size(arr_img)
          dim_array=tmp(0)
          
          IF dim_array EQ 3 THEN BEGIN
             num_img=N_ELEMENTS(arr_img(*,0,0))*(tr_change/tr_img)
             img_changed_tim=STRARR(num_img)
             ms_changed=DBLARR(num_img)
             img=DBLARR(num_img,256,256)             
             count=0L
             i=0L

             WHILE count LT num_img DO BEGIN
                PRINT,'changing time resolution of img:'+STRING(count)+'/'+STRING(num_img-1,FORMAT='(I4.4)')

                IF i LT ave_img_val/2 THEN key=1
                IF (i GE ave_img_val/2) AND (i+ave_img_val/2 LT N_ELEMENTS(arr_img(*,0,0))) THEN key=2
                IF i+ave_img_val/2 GE N_ELEMENTS(arr_img(*,0,0)) THEN key=3

                CASE key OF
                   1 : BEGIN                      
                      img(count,*,*)=MEAN(arr_img(i : i+ave_img_val-1,*,*),DIMENSION=1)
                   END

                   2 : BEGIN
                      s = ave_img_val/2
                      img(count,*,*)=MEAN(arr_img(i-s : i+s-1,*,*),DIMENSION=1)
                   END

                   3 : BEGIN
                      s = ave_img_val
                      img(count,*,*)=MEAN(arr_img(i-s : i-1,*,*),DIMENSION=1)
                   END
                ENDCASE

                img_changed_tim(count)=STRING(yr_raw(i),FORMAT='(I4.4)')+'-'+STRING(mo_raw(i),FORMAT='(I2.2)')+'-'+STRING(dy_raw(i),FORMAT='(I2.2)')+'/'+STRING(hh_raw(i),FORMAT='(I2.2)')+':'+STRING(mm_raw(i),FORMAT='(I2.2)')+':'+STRING(ss_raw(i),FORMAT='(I2.2)')

                ms_changed(count)=ms_raw(i)
                
                count++
                IF (count MOD 2) EQ 0 THEN i=i+ave_img_val ELSE i=i+(ave_img_val+1)
        
             ENDWHILE

             img_changed_tim_ms=time_double(img_changed_tim)+DOUBLE(ms_changed/1e3)
             
          ENDIF ELSE BEGIN
             num_img=N_ELEMENTS(arr_img)*(tr_change/tr_img)
             img_changed_tim_ms=STRARR(num_img)
             img=DBLARR(num_img)
             count=0L
             i=0L

             WHILE count LT num_img DO BEGIN
                PRINT,'changing time resolution of img:'+STRING(count)+'/'+STRING(num_img-1,FORMAT='(I3.3)')
                
                SWITCH 1 OF
                   (i LT ave_img_val/2): BEGIN                      
                      FOR k=0,ave_img_val-1 DO BEGIN
                         img(count)+=arr_img(i+k)
                      ENDFOR
                      img(count)=img(count)/(ave_img_val+1)
                   END

                   ((i GE ave_img_val/2) AND (i+ave_img_val/2 LT N_ELEMENTS(arr_img))) : BEGIN
                      s = ave_img_val/2
                      FOR k=0,ave_img_val-1 DO BEGIN
                         img(count)+=arr_img(i-s+k)
                      ENDFOR
                      img(count)=img(count)/(ave_img_val+1)
                   END

                   (i+ave_img_val/2 GE N_ELEMENTS(arr_img)) : BEGIN
                      s= ave_img_val
                      FOR k=0,ave_img_val-1 DO BEGIN
                         img(count)+=arr_img(i-s+k)
                      ENDFOR
                      img(count)=img(count)/(ave_img_val+1)
                   END
                ENDSWITCH

                img_changed_tim_ms(count)=STRING(yr_raw(i),FORMAT='(I4.4)')+'-'+STRING(mo_raw(i),FORMAT='(I2.2)')+'-'+STRING(dy_raw(i),FORMAT='(I2.2)')+'/'+STRING(hh_raw(i),FORMAT='(I2.2)')+':'+STRING(mm_raw(i),FORMAT='(I2.2)')+':'+STRING(ss_raw(i),FORMAT='(I2.2)') +':'+STRING(ms_raw(i),FORMAT='(I3.3)')
                
                count++
                i=i+ave_img_val
             ENDWHILE                        
          ENDELSE
          RETURN,img
       END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	 MERGED_PLOT
;
        PRO merged_plot,day,hhmmss,dura,span

           COMMON prf_psa
           COMMON fov_psa
           COMMON raw_psa
           COMMON cor_psa
           COMMON struct_psa

           IF NOT KEYWORD_SET(day) THEN BEGIN
              day=''
              hhmmss=''
              
              READ,day,PROMPT='Enter day(YYYYMMDD):'
              READ,hhmmss,PROMPT='Enter start time(HHMMSS):'
              READ,dura,PROMPT='Enter duration(hour):'
              READ,span,PROMPT='Enter calculation span of FFT(second):'

           ENDIF
                      
              YYYY=STRMID(day,0,4)
              MM=STRMID(day,4,2)
              DD=STRMID(day,6,2)
              hhmm=STRMID(hhmmss,0,4)
              hour=STRMID(hhmm,0,2)
              min=STRMID(hhmm,2,2)
              second=STRMID(hhmmss,4,2)

              ;window_psa,/landscape
              !P.MULTI=[0,3,2]

              PLOT_COR_LINE,day,hhmmss,dura,span,merged=1
              

              ;!P.position=[0.55,0.3,0.9,0.55]

              PLOT_FFT_LINE,day,hhmmss,span,no_img_line=1
              
              !P.position=[0.55,0.54,0.87,0.94]
              PLOT_COR_CONTOUR,day,hhmmss,dura,span,no_plot=1

              map2d_init
              map2d_coord,'geo'
              ;scale=10E6*6.8
              map2d_set,glatc=67.25,glonc=26.35,scale=10E5*4.4,position=!P.POSITION

              gla_sta=[67.25]
              glo_sta=[26.35]
              let_sta=['SOD']

              loadct,[0]
              contour_lines=30
              levels_set=(max_psa-min_psa)*FINDGEN(contour_lines)/FLOAT(contour_lines)+min_psa
              colour_set=!D.TABLE_SIZE*INDGEN(contour_lines)/contour_lines
              x_cont=geogra_map.x
              y_cont=geogra_map.y
              z_cont=geogra_map.z
              
              CONTOUR,z_cont,x_cont,y_cont,/OVERPLOT,/FILL,LEVELS=levels_set,C_COLORS=colour_set,/IRREGULAR

; Contouring
              max_cor=1.0
              min_cor=0.0
              loadct,[39]
              contour_lines=30
              levels_set=(max_cor-min_cor)*FINDGEN(contour_lines)/FLOAT(contour_lines)+min_cor
              colour_set=!D.TABLE_SIZE*INDGEN(contour_lines)/contour_lines
              x_cont=cor_contour.x
              y_cont=cor_contour.y
              z_cont=cor_contour.z
              
              CONTOUR,z_cont,x_cont,y_cont,/OVERPLOT,/CELL_FILL,LEVELS=levels_set,C_COLORS=colour_set,/IRREGULAR,MIN_VALUE=0.1
              
; Get orbit data
              timespan,YYYY_orb+'-'+MM_orb+'-'+DD_orb+'/'+hour_orb+':'+min_orb+':00',dura_orb,/hour
              get_data,'rbspa_ifoot_geo_lat',data=data
              IF NOT KEYWORD_SET(data) THEN GET_ORBIT_DATA,YYYY_orb+MM+DD,hour_orb+min_orb,dura_orb
        
; overlay satelite orbit
              overlay_map_coast
              overlay_map_sc_ifoot,'rbspa_ifoot_geo_lat','rbspa_ifoot_geo_lon'
              ;timespan,YYYY+'-'+MM+'-'+DD+'/'+hour+':'+min+':00',dura,/hour
        
; plot station and satelite position
        
              PSYM8,/FILL
              OPLOT,[glo_sta],[gla_sta],PSYM=3,SYMSIZE=0.5,COL=fsc_color("black")
              OPLOT,[glo_sta],[gla_sta],PSYM=1,SYMSIZE=1.0,COL=fsc_color("black")
              PLOT_POS_FOOTPRINT
              XYOUTS,glo_sta,gla_sta+0.2,let_sta,COL=fsc_color("black"),ALI=0.5,CHARSIZE=0.8*!P.CHARSIZE

              OPLOT,[glo_psa(max_pixel(0),max_pixel(1))],[gla_psa(max_pixel(0),max_pixel(1))],PSYM=8,SYMSIZE=1.0,COL=fsc_color("Red")

              ;XYOUTS,0.4,0.01,'Max Correlate Coefficient:'+STRING(max_pixel_cor),CHARSIZE=2.0,/NORMAL
        
; Plot messages
              pos=[0.05,0,0,0.97]
              XYOUTS,pos(0),pos(3),STRMID(day,0,4)+'-'+STRMID(day,4,2)+'-'+STRMID(day,6,2)+' '+STRMID(hhmmss,0,2)+':'+STRMID(hhmmss,2,2)+':'+STRMID(hhmmss,4,2)+' 5 min : Cam'+STRING(cid_raw,FORMAT='(I1.1)'),/NORMAL,charsize=2.0

; Plot colour bar
              position=[!P.POSITION(2)+0.01,!P.POSITION(1),!P.POSITION(2)+0.02,!P.POSITION(3)]

              plot_colour_bar,position=position,legend=legend,bar_title='cor',max_val=max_cor,min_val=min_cor


              !P.position=[0.1,0.1,0.5,0.5]
; plot geogra map              
              ;map2d_init
              ;map2d_coord,'geo'
              ;map2d_set,glatc=67.25,glonc=26.35,scale=scale

              ;gla_sta=[67.25]
              ;glo_sta=[26.35]
              ;let_sta=['SOD']

              ;loadct,[0]
              ;contour_lines=30
              ;levels_set=(max_psa-min_psa)*FINDGEN(contour_lines)/FLOAT(contour_lines)+min_psa
              ;colour_set=!D.TABLE_SIZE*INDGEN(contour_lines)/contour_lines
              ;CONTOUR,geogra_map.z,geogra_map.x,geogra_map.y,/OVERPLOT,/FILL,LEVELS=levels_set,C_COLORS=colour_set,/IRREGULAR

              ;PLOT_FFT_CONTOUR,day,hhmmss,dura,span
           END
;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	 SAVE_MERGED_PLOT
;
        PRO save_merged_plot,day,hhmmss,dura,span

           COMMON prf_psa
           COMMON fov_psa
           COMMON raw_psa
           COMMON cor_psa
           COMMON struct_psa

           IF NOT KEYWORD_SET(day) THEN BEGIN
              day=''
              hhmmss=''
              
              READ,day,PROMPT='Enter day(YYYYMMDD):'
              READ,hhmmss,PROMPT='Enter start time(HHMMSS):'
              READ,dura,PROMPT='Enter duration(hour):'
              READ,span,PROMPT='Enter calculation span of FFT(second):'

           ENDIF
                      
              YYYY=STRMID(day,0,4)
              MM=STRMID(day,4,2)
              DD=STRMID(day,6,2)
              hhmm=STRMID(hhmmss,0,4)
              hour=STRMID(hhmm,0,2)
              min=STRMID(hhmm,2,2)
              second=STRMID(hhmmss,4,2)

              openw,lun,'/home/suguru/spedas/image_data/txt_period.txt',/get_lun              

              FOR i=0,60/cal_span*dura*60-1 DO BEGIN

                 IF LONG(STRMID(hhmmss,2,2))-min EQ 1 THEN BEGIN
                    FREE_VARIABLE
                    file_psa_raw,2,day,hhmm,span/60+2
                 ENDIF
                 
                 MERGED_PLOT,day,hhmmss,dura,span

                 makepng,'/home/suguru/spedas/image_data/merged_plot_8Hz/'+day+'-'+hhmmss+'(5 min)'

                 ; Set next time
                 min=STRMID(hhmmss,2,2)
                 hhmmss = STRING(LONG(hhmmss)+cal_span,FORMAT='(I6.6)')
                 IF STRMID(hhmmss,4,2) EQ '60' THEN BEGIN
                    hhmmss = STRING(LONG(hhmmss)+100-60,FORMAT='(I6.6)')                
                 ENDIF              
              
                 IF STRMID(hhmmss,2,2) EQ '60' THEN BEGIN
                    hhmmss = STRING(LONG(hhmmss)+10000-6000,FORMAT='(I6.6)')
                 ENDIF

                 IF STRMID(hhmmss,0,2) EQ '24' THEN BEGIN
                    day = STRING(LONG(day)+1,FORMAT='(I8.8)')
                    hhmmss = '000000'
                 ENDIF
                 hhmm=STRMID(hhmmss,0,4)

              ENDFOR
              free_lun,lun

           END
        
;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	 FREE VARIABLE
;
        PRO free_variable

           COMMON prf_psa
           COMMON fov_psa
           COMMON raw_psa
           COMMON cor_psa

           img_raw=!NULL
           num_raw=!NULL
           day_raw=!NULL
           yr_raw=!NULL
           mo_raw=!NULL
           dy_raw=!NULL
           hh_raw=!NULL
           mm_raw=!NULL
           ss_raw=!NULL
           ms_raw=!NULL
           ts_raw=!NULL
           ti_raw=!NULL
           IF KEYWORD_SET(cor) THEN cor=!NULL
           IF KEYWORD_SET(img_changed) THEN img_changed=!NULL
           IF KEYWORD_SET(img_changed_tim) THEN img_changed_tim=!NULL
           IF KEYWORD_SET(img_changed_tim_ms) THEN img_changed_tim_ms=!NULL
           IF KEYWORD_SET(fbk_data_changed) THEN fbk_data_changed=!NULL
           IF KEYWORD_SET(fbk_tim_changed) THEN fbk_tim_changed=!NULL

        END
        
;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	 REPLICATE_VECTOR
;
         FUNCTION REPLICATE_VECTOR, vector, numberReplicates, COLUMNS=columns

           IF KEYWORD_SET( columns ) THEN BEGIN
              RETURN, vector ## REPLICATE( 1, numberReplicates, 1 )
           ENDIF ELSE BEGIN
              RETURN, REPLICATE( 1, numberReplicates ) ## vector
           ENDELSE

        END
;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	 NORMALIZE
;
         FUNCTION NORMALIZE, vector

;+
; Determine the dimensionality of the vector(s) provided
;-
           dimensions = SIZE( vector, /DIMENSIONS )
           numberElements = dimensions[0]

;+
; Create a double-precision local copy of the provided vector(s) for use in the
; normalization computation
;-
           v = DOUBLE( vector )

;+
; Replicate the magnitude in a vector the same size as the provided vector dimensions
;-
           m = REPLICATE_VECTOR( MAGNITUDE( v ), numberElements, /COLUMNS )

;+
; Return the normalized vector(s)
;-
           RETURN, v / m

        END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	 MAGNITUDE
;
         FUNCTION MAGNITUDE, vector

;+
; Create a double-precision local copy of the provided vector(s) for use in the
; computation of the magnitude
;-
           v = DOUBLE( vector )
   
;+
; Return the magnitude(s) of the provided vector(s)
;-
           RETURN, SQRT( TOTAL( v * v, 1 ) )

        END
     
;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	 COHERENCE
;
        PRO coherence,day,hhmmss,dura,span
          
           COMMON prf_psa
           COMMON fov_psa
           COMMON raw_psa
           COMMON cor_psa

           IF NOT KEYWORD_SET(day) THEN BEGIN
              day=''
              hhmmss=''
              
              READ,day,PROMPT='Enter day(YYYYMMDD):'
              READ,hhmmss,PROMPT='Enter start time(HHMMSS):'
              READ,dura,PROMPT='Enter duration(hour):'
              READ,span,PROMPT='Enter calculation span of FFT(second):'

           ENDIF
                      
              YYYY=STRMID(day,0,4)
              MM=STRMID(day,4,2)
              DD=STRMID(day,6,2)
              hhmm=STRMID(hhmmss,0,4)
              hour=STRMID(hhmm,0,2)
              min=STRMID(hhmm,2,2)
              second=STRMID(hhmmss,4,2)
              
              tmp1=ROLLING_THUNDER(day,hhmmss,dura,span)
; Cut img data
           num_img=span*8
           start_time = YYYY+'-'+MM+'-'+DD+'/'+hour+':'+min+':'+second
           s_range=WHERE(img_changed_tim EQ start_time)
           s_range=s_range(0)
              
           img_changed_cut=img_changed(s_range : s_range+num_img-1,*,*)
           img_changed_tim_ms_cut=img_changed_tim_ms(s_range : s_range+num_img-1)
           
           store_data,'fbk_data',data={x:fbk_tim_changed,y:fbk_data_changed}
           store_data,'img_data',data={x:img_changed_tim_ms_cut,y:REFORM(img_changed_cut(*,max_pixel(0),max_pixel(1)))}
           ;thigh_path_filter,'fbk_data',30
           ;thigh_path_filter,'img_data',30
           
           timespan,img_changed_tim_ms_cut(0),5,/min

           uspec_coh_sod,'img_data','fbk_data',deltat=0.125

	   XYOUTS,0.1,0.97,YYYY+'-'+MM+'-'+DD+'/'+hour+':'+min+':'+second+' 5 min',charsize=2.0,/NORMAL
           
        END
;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	 SAVE_COHERENCE
;
        PRO save_coherence
          
           COMMON prf_psa
           COMMON fov_psa
           COMMON raw_psa
           COMMON cor_psa

           IF NOT KEYWORD_SET(day) THEN BEGIN
              day=''
              hhmmss=''
              
              READ,day,PROMPT='Enter day(YYYYMMDD):'
              READ,hhmmss,PROMPT='Enter start time(HHMMSS):'
              READ,dura,PROMPT='Enter duration(hour):'
              READ,span,PROMPT='Enter calculation span of FFT(second):'

           ENDIF
                      
              YYYY=STRMID(day,0,4)
              MM=STRMID(day,4,2)
              DD=STRMID(day,6,2)
              hhmm=STRMID(hhmmss,0,4)
              hour=STRMID(hhmm,0,2)
              min=STRMID(hhmm,2,2)
              second=STRMID(hhmmss,4,2)

              FOR i=0,(dura*60)/(cal_span/60)-1 DO BEGIN

                 IF FIX(STRMID(hhmmss,2,2))-FIX(min) NE 0 THEN BEGIN
                    FREE_VARIABLE
                 
                    file_psa_raw,2,day,hhmm,span/60+2
                 ENDIF
                 

                 PRINT,'Current time: '+day+'-'+hhmmss
                 
		 COHERENCE,day,hhmmss,dura,span
                 makepng,'/home/suguru/spedas/image_data/coh_sod/'+day+hhmmss+'(5 min)'
; Set next time range
                 min=STRMID(hhmmss,2,2)
                 hhmmss = STRING(LONG(hhmmss)+cal_span,FORMAT='(I6.6)')
                 IF STRMID(hhmmss,4,2) EQ '60' THEN BEGIN
                    hhmmss = STRING(LONG(hhmmss)+100-60,FORMAT='(I6.6)')
                 ENDIF
                 
              
                    IF STRMID(hhmmss,2,2) EQ '60' THEN BEGIN
                       hhmmss = STRING(LONG(hhmmss)+10000-6000,FORMAT='(I6.6)')
                    ENDIF

                    IF STRMID(hhmmss,0,2) EQ '24' THEN BEGIN
                       day = STRING(LONG(day)+1,FORMAT='(I8.8)')
                       hhmmss = '000000'
                    ENDIF

                    hhmm=STRMID(hhmmss,0,4)
                    ;count=count+tr_img
                    ;GO_PSA_RAW,count                 
		ENDFOR
        END
       
              
;------------------------------------------------------------------------------------------------------------------
