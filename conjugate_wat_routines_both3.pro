;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	COMMON BLOCK
;

	COMMON	syowa,	        img_wat,hms_wat,tsc_wat,num_wat,now_wat,keo_wat, $
				max_wat,min_wat,siz_wat,stm_wat,etm_wat,tim_wat, $
				day_wat,a_v_wat,rti_wat,rog_wat,rom_wat,gla_wat, $
				glo_wat,mla_wat,mlo_wat,yer_wat,jdy_wat,yrs_wat, $
				cx1_wat

        COMMON	tjornes,	img_wat2,hms_wat2,tsc_wat2,num_wat2,now_wat2,keo_wat2, $
				max_wat2,min_wat2,siz_wat2,stm_wat2,etm_wat2,tim_wat2, $
				day_wat2,a_v_wat2,rti_wat2,rog_wat2,rom_wat2,gla_wat2, $
				glo_wat2,mla_wat2,mlo_wat2,yer_wat2,jdy_wat2,yrs_wat2, $
                                zan_wat2

        COMMON  fft,            num_data,fft_tjo,fft_syo, $
                                plot_wat,plot_wat2,x_plot,x_plot2


;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       RADAR_TO_GROUND_RANGE
;
; NOTE:
;
;       Function for converting the radar range to the ground range
;

        FUNCTION radar_to_ground_range,radar_range,altitude

; Earth radius
        Re=6370.0

; Just using "YOGEN-TEIRI"
        sin_el=((altitude+Re)^2-Re^2-radar_range^2)/(2*Re*radar_range)
        ground_range=Re*ATAN(SQRT(1-sin_el^2)/(Re/radar_range+sin_el))
        
        RETURN,ground_range

        END

;-------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       FIND_CELL
;
; NOTE:
;
;       Function for finding radar range gates without rbpos
;

        FUNCTION find_cell,radarlat,radarlon,nobeams,boresite,beamsep,gatelen,firstrange,beamno,gateno

; Earth radius
        Re=6370.0

        lon=radarlon*!DTOR
        colat=(90-radarlat)*!DTOR
        azimuth=((((boresite+beamsep*(beamno-0.5*(nobeams-1)))+900) MOD 360)-180)*!DTOR
        range=firstrange+gateno*gatelen
        
; Uncomment this line to take account of ground range being less than radar range
        range=radar_to_ground_range(range,110.0)

; Just using "KYUMEN-SANKAKU"
        ccolat=ACOS(COS(range/Re)*COS(colat)+SIN(range/Re)*SIN(colat)*COS(azimuth))
        clon=ACOS(MIN([(COS(range/Re)-COS(ccolat)*COS(colat))/(SIN(ccolat)*SIN(colat)), $
                1.0]))
        clat=90-ccolat/!DTOR
        IF azimuth GT 0 THEN clon=(lon+clon)/!DTOR ELSE clon=(lon-clon)/!DTOR
        
        RETURN,[clat,clon]

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

;-------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       GET_MAPPED_POINT
;
; PURPOSE:
;
;       Calculate ionospheric point in lat and lon for given ele and azm
;       assuming altitude of the ionosphere.
;

        FUNCTION get_mapped_point,ele,azm,alt

        lat_tjo=66.20
        lon_tjo=342.88

        lat_and_lon=ang_to_loc(lat_tjo,lon_tjo,azm,90-ele,alt)

        RETURN,lat_and_lon

        END

;-------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	GET_GROUND_DISTANCE_FROM_TJORNES
; 
;       Get distance (on the surface of the Earth) in km from Tjornes to given point.
;

	FUNCTION get_ground_distance_from_tjornes,lat_giv,lon_giv

; Location of Tjornes
	lat_tjo= 66.20
	lon_tjo=342.88

; Radius of the Earth
        Re=6370.

; Co-latitude from pole to Tjornes
        colata=(90.0-lat_tjo)/180.*!PI

; Co-latitude from pole to point B
        colatb=(90.0-lat_giv)/180.*!PI

; Angle between the lines from pole to A and from pole to B
        inpro=ABS(lon_giv-lon_tjo)/180.0*!PI

; Distance between point A and point B
        ground_dist=Re*ACOS(COS(colata)*COS(colatb)+SIN(colata)*SIN(colatb)*COS(inpro))

; Return the results
        RETURN,ground_dist

	END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	GET_DISTANCE_FROM_TJORNES
;

	FUNCTION get_distance_from_tjornes,lat_giv,lon_giv,height

; Radius of the Earth
        Re=6370.

; Get ground distance from Tjornes
	ground_distance=get_ground_distance_from_tjornes(lat_giv,lon_giv)

; Get offset angle on the big circle including Tjornes and given point
	ang_offset=360.0*ground_distance/(2*!PI*Re)

; Get straigh line distance
	distance=SQRT(Re^2+(Re+height)^2-2*Re*(Re+height)*COS(!PI*ang_offset/180.0))

	RETURN,distance

	END

;-------------------------------------------------------------------------------------------------------------------
; NAME:
;
; 	GET_ZENITH_ANGLE_FROM_TJORNES
;

	FUNCTION get_zenith_angle_from_tjornes,lat_giv,lon_giv,height

; Radius of the Earth
        Re=6370.

; Get straight line distance from Tjornes to the point specified
	distance=get_distance_from_tjornes(lat_giv,lon_giv,height)

; Get elevation angle
	angle1=ACOS((distance^2-height^2-2*Re*height)/(2*distance*Re))*180.0/!PI
	zenith_angle=180.0-angle1

	RETURN,zenith_angle

	END

;-------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	GET_AZIMUTH_FROM_TJORNES
;

	FUNCTION get_azimuth_from_tjornes,lat_giv,lon_giv

; Location of Tjornes
	lat_tjo= 66.20
	lon_tjo=-17.12

; Radius of the Earth
        Re=6370.

; Angle between the lines from pole to A and from pole to B
;        small_c=!PI*(lon_giv-lon_tjo)/180.0
        large_c=!PI*(lon_giv-lon_tjo)/180.0

; Get colatitude of Tjornes and given point
	colat_tjo=90.0-lat_tjo
	colat_giv=90.0-lat_giv
	small_b=!PI*colat_tjo/180.0
	small_a=!PI*colat_giv/180.0

; Get ground distance from Tjornes
	cos_small_c=COS(small_a)*COS(small_b)+SIN(small_a)*SIN(small_b)*COS(large_c)
	small_c=ACOS(cos_small_c)

; Get azimuth angle
	cos_large_a=(COS(small_a)-COS(small_b)*COS(small_c))/(SIN(small_b)*SIN(small_c))
	azimuth=180.0*ACOS(cos_large_a)/!PI
	IF lon_giv LT lon_tjo THEN azimuth=-azimuth
	IF azimuth LT 0 THEN azimuth=azimuth+360.0

	RETURN,azimuth

	END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	GET_PIXEL
;
; PURPOSE:
;
;	Get x and y of pixel of the specified zenith and azimuth angle.
;

	FUNCTION get_pixel,zan,azm

	COMMON syowa

; Values of fish-eye lens define with define_zenith procedure
	X0=siz_wat/2
	Y0=siz_wat/2

	azm2=360-azm
	azm3=azm2 MOD 360
	azm4=azm3-rog_wat
	x=FIX(X0+(0.5*siz_wat*zan/75.0)*SIN(azm4*!PI/180.0))
	y=FIX(Y0+(0.5*siz_wat*zan/75.0)*COS(azm4*!PI/180.0))
	
	RETURN,[x,y]

	END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	GET_ZAN_AZM
;
; PURPOSE:
;
;	Get zenith and azimuth angle of specified x and y pixel on the image.
;

	FUNCTION get_zan_azm,x,y

	COMMON syowa

; Values of fish-eye lens define with define_zenith procedure
	X0=siz_wat/2
	Y0=siz_wat/2

; Displacement from zenith
	dx=FLOAT(x-X0)
	dy=FLOAT(y-Y0)

; Convert (dx,dy) into zenith and azimuth
	zan=75.0*SQRT(dx^2+dy^2)/FLOAT(siz_wat*0.5)
	azm0=90.0-180.0*ATAN(dy,dx)/!PI
	azm1=azm0+rog_wat
	azm2=360.0-azm1
	IF azm2 LT   0.0 THEN azm2=azm2+360.0
	IF azm2 GT 360.0 THEN azm2=azm2-360.0

	RETURN,[zan,azm2]
	;RETURN,[ROUND(zan),ROUND(azm2)]

	END

;-------------------------------------------------------------------------------------------------------------------
; NAME
;
;	FILE_WAT
;

	PRO file_wat,path=path,retry=retry

	COMMON syowa

; Imomusi data directory
	IF NOT KEYWORD_SET(path) THEN path='/raid/data/Aurora/Tjorwat/20110930_htr_tks/syo/'
	;IF NOT KEYWORD_SET(path) THEN path='/raid/data/Aurora/Tjorwat/20110930_htr_all/'

; Search the data files
        found_file=FILE_SEARCH(path+'*.jpg',COUNT=num_found_file)
        num_wat=num_found_file
        print,num_wat

        
; Values of fish-eye lens define with define_zenith procedure (see alignment.txt in detail)
        X0=318         ; X of central pixel 320
        Y0=248         ; Y of central pixel 240
	a_v_wat=130    ; A-value (diameter/3.1415)
        rog_wat=24.0   ; GN is rotated to the East by this angle in degree
	rom_wat=-16.0 ; Declination at SYOWA (taken from 元場さん) in 2009

; Grid adjustment
	cnt_wat=[X0,Y0]
	siz_wat=FIX(a_v_wat*!PI*75.0/90.0)            ; zenith angle limit of 75 deg
	IF siz_wat MOD 2 EQ 1 THEN siz_wat=siz_wat+1  ; siz_wat must be an even number

; Initialize arrays

	img_wat1=INTARR(num_wat,siz_wat+1,siz_wat+1)
        img_wat=INTARR(num_wat,siz_wat+1,siz_wat+1)
        keo_wat=INTARR(10,num_wat,siz_wat+1)
	day_wat=STRARR(num_wat)
	yer_wat=STRARR(num_wat)
	yrs_wat=STRARR(num_wat)
	jdy_wat=STRARR(num_wat)
	tsc_wat=LONARR(num_wat)
	tim_wat=STRARR(num_wat)
	hms_wat=STRARR(num_wat)

; Read the data
	FOR i=0,num_wat-1 DO BEGIN

                hh=STRMID(found_file(i),STRLEN(path)+17,2)
                IF hh EQ 00 THEN hh=STRING(24)
                mm=STRMID(found_file(i),STRLEN(path)+20,2)
                ss=STRMID(found_file(i),STRLEN(path)+23,2)
                tsc_wat(i)=3600L*FIX(hh)+60L*FIX(mm)+FIX(ss) ; 通し秒 of 一日
                hms_wat(i)=hh+mm+ss
                READ_JPEG,found_file(i),tmp_read,/TRUE
                img_tmp=REFORM(tmp_read(0,*,*))
                print,SIZE(img_tmp,/DIM)

; Get time info
                yr=STRMID(found_file(i),STRLEN(path)+08,2)
                mo=STRMID(found_file(i),STRLEN(path)+11,2)
                dy=STRMID(found_file(i),STRLEN(path)+14,2)
                yer_wat(i)=2000+FIX(yr)
                day_wat(i)='20'+yr+mo+dy
                tim_wat(i)=hh+mm+' '+ss+'s UT'
                
                jdy_wat(i)=JULDAY(FIX(mo),FIX(dy),yer_wat(i))-JULDAY(1,0,yer_wat(i))
		yrs_wat(i)=86400L*(jdy_wat(i)-1)+tsc_wat(i)
                PRINT,'Scan'+STRING(i,FORMAT='(I5)')+': '+day_wat(i)+' '+tim_wat(i)+ $
			' ('+STRMID(found_file(i),STRLEN(path),32)+')'

; Get 2D image ,decide cutting circle
                img_wat(i,*,*)=img_tmp(cnt_wat(0)-siz_wat/2:cnt_wat(0)+siz_wat/2, $
                                        cnt_wat(1)-siz_wat/2:cnt_wat(1)+siz_wat/2)
                ;img_wat(i,*,*)=ROT(img_wat1,24,/INTERP)
                ;print,SIZE(img_wat,/DIM)
                ;loadct,[0]
                ;TV,img_wat

; Get keogram from GS to GN (Line ID 0)
		ang_for_keo=rog_wat
		offset=FIX(TAN(ang_for_keo*!PI/180.0)*(siz_wat/2))
		FOR k=0,siz_wat DO BEGIN
			x_keo=FIX(siz_wat/2+offset-k*TAN(ang_for_keo*!PI/180.0))
			y_keo=k
			keo_wat(0,i,k)=img_wat(i,x_keo,y_keo)
		ENDFOR

; Get keogram from GE to GW (Line ID 1)
		ang_for_keo=rog_wat
		offset=FIX(TAN(ang_for_keo*!PI/180.0)*(siz_wat/2))
		FOR k=0,siz_wat DO BEGIN
			x_keo=k
			y_keo=FIX(siz_wat/2-offset+k*TAN(ang_for_keo*!PI/180.0))
			keo_wat(1,i,k)=img_wat(i,x_keo,y_keo)
		ENDFOR

; Get keogram from MS to MN (Line ID 2)
		ang_for_keo=rog_wat+rom_wat
		offset=FIX(TAN(ang_for_keo*!PI/180.0)*(siz_wat/2))
		FOR k=0,siz_wat DO BEGIN
			x_keo=FIX(siz_wat/2+offset-k*TAN(ang_for_keo*!PI/180.0))
			y_keo=k
			keo_wat(2,i,k)=img_wat(i,x_keo,y_keo)
		ENDFOR

; Get keogram from ME to MW (Line ID 3)
		ang_for_keo=rog_wat+rom_wat
		offset=FIX(TAN(ang_for_keo*!PI/180.0)*(siz_wat/2))
		FOR k=0,siz_wat DO BEGIN
			x_keo=k
			y_keo=FIX(siz_wat/2-offset+k*TAN(ang_for_keo*!PI/180.0))
			keo_wat(3,i,k)=img_wat(i,x_keo,y_keo)
		ENDFOR

; Get keogram (2: beam 7 align)
                boresite_from_tjornes=32.0 ; NOTE: This is not the radar boresite! Azimuth of the beam 7 from TJO.
              ;  ang_for_keo=rog_wat+boresite_from_tjornes-3.24*0.5
              ;  offset=FIX(TAN(ang_for_keo*!PI/180.0)*(siz_wat/2))
	;	FOR k=0,siz_wat DO BEGIN
	;		x_keo=FIX(siz_wat/2+offset-k*TAN(ang_for_keo*!PI/180.0))
	;		y_keo=k
	;		keo_wat(4,i,k)=img_wat(i,x_keo,y_keo)
	;	ENDFOR

; Get Kishiyama's keogram (for syowa)
              ;  ang_for_keo=rog_wat+rom_wat
               ; offset=FIX(TAN(ang_for_keo*!PI/180.0)*(siz_wat/2))
                ;FOR k=0,siz_wat DO BEGIN
                 ;  x_keo=FIX(siz_wat/2+offset-k*TAN(ang_for_keo*!PI/180.0))
                  ; y_keo=k
                   ;keo_wat(5,i,k)=img_wat(i,x_keo,y_keo)
                                ; ENDFOR
             ENDFOR
        

; Summary of settings
	PRINT,''
	PRINT,'Summary of current settings'
	PRINT,'------------------------------------------'

; Current image
	go_wat,0L

; Start and end time
	time_wat,hms_wat(0),hms_wat(num_wat-1)
	IF STRMID(hms_wat(num_wat-1),3,1) EQ '9' THEN time_wat,hms_wat(0),hms_wat(num_wat-1)

; Set scale
	set_scale_wat,0,200
	PRINT,'------------------------------------------'
	PRINT,''

; Set default charsize
	!P.CHARSIZE=1.0

	END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;
;       GET_BOTH_WAT

        PRO get_both_wat,time;,xmaps,ymaps
          

          COMMON syowa
          COMMON tjornes

          ;FOR t=time,235959 DO BEGIN
          print,time
          GO_TIME_WAT,time
          GO_TIME_WAT2,time
          mm=STRMID(time,0,4)
          ss=STRMID(time,4,2)
          
          
          XYOUTS,240,siz_wat2+50,20111001,$
                 /DATA,COL=foreground,CHARSIZE=1.4 ;charsize kaetayo
          XYOUTS,siz_wat2+250,siz_wat2+50,mm+' '+ss+'s UT',$
                 /DATA,COL=foreground,ALI=1,CHARSIZE=1.4
         
          PLOT_BOTH             ;,xmaps,ymaps

          ;ENDFOR




          END

        
;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       SAVE_PNG_WAT
;

        PRO save_png_wat,stime,etime

          COMMON tjornes
          COMMON syowa

          
          t=stime+10
          print,t
          FOR i=0,10 DO BEGIN
           get_both_wat,'00000'+i
           fname='pic_movie1/'+STRING(i, FORMAT='(I6.6)')+'.png'
           write_png,fname,tvrd(true=1)
           clear_page           
        ENDFOR   
           
 
       END
  




;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       SAVE_ONE
;

        PRO save_one,time

          COMMON tjornes
          COMMON syowa

  

         ; FOR i=time,time2 DO BEGIN
           get_both_wat,time
           print,time
           fname='pic_movie1/'+STRING(time, FORMAT='(I6.6)')+'.png'
           write_png,fname,tvrd(true=1)
           clear_page
           
         ; ENDFOR   
           
 
       END
  






        







        

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	MAKE_FOV_WAT
;
	
	PRO make_fov_wat

	COMMON syowa

	height=110.0

; Initialize array
	gla_wat=FLTARR(siz_wat+1,siz_wat+1)
	glo_wat=FLTARR(siz_wat+1,siz_wat+1)
	mla_wat=FLTARR(siz_wat+1,siz_wat+1)
	mlo_wat=FLTARR(siz_wat+1,siz_wat+1)
	zan_wat=FLTARR(siz_wat+1,siz_wat+1)

	;PLOT,[0],[0],XRANGE=[330,360],YRANGE=[60,75]
	FOR xx=0,siz_wat DO BEGIN
		FOR yy=0,siz_wat DO BEGIN

			tmp=get_zan_azm(xx,yy)
			zan=tmp(0)
			azm=tmp(1)
			IF zan LT 0.0 OR zan GT 75.0 THEN BEGIN
				zan=99999.9
				azm=99999.9
			ENDIF
			IF FINITE(zan,/NAN) THEN BEGIN
				zan=99999.9
				azm=99999.9
			ENDIF
			zan_wat(xx,yy)=zan
			IF zan NE 99999.9 AND azm NE 99999.9 THEN BEGIN
        			tmp=get_mapped_point(90-zan,azm,height)
				;!P.SYMSIZE=1
				;OPLOT,[tmp(1)],[tmp(0)],COL=252,SYMSIZE=0.2*!P.SYMSIZE,PSYM=1
				gla_wat(xx,yy)=tmp(0)
				glo_wat(xx,yy)=tmp(1)
				tmp=cnvcoord(gla_wat(xx,yy),glo_wat(xx,yy),height)
				mla_wat(xx,yy)=tmp(0)
				mlo_wat(xx,yy)=tmp(1)
			ENDIF ELSE BEGIN
				gla_wat(xx,yy)=99999.9
				glo_wat(xx,yy)=99999.9
				mla_wat(xx,yy)=99999.9
				mlo_wat(xx,yy)=99999.9
			ENDELSE

		ENDFOR
	ENDFOR

	END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	GO_WAT
;
; EXAMPLE:
;
;       Go > go_wat,10
;

	PRO go_wat,img_to_go

	COMMON syowa

; Set new current image
	now_wat=0L
	now_wat=img_to_go

; Verbose
        current_hms=STRING(hms_wat(now_wat),FORMAT='(I6.6)')
        PRINT,'Time of current image: '+current_hms+' UT'

	END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	NEXT_IMAGE
;

	PRO next_image

	COMMON syowa

	go_wat,now_wat+1

	END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       GO_TIME_WAT
;
; EXAMPLE:
;
;       Go > go_time_wat,025800
;

        PRO go_time_wat,image_time

        COMMON syowa

; Finding image
        found_image=WHERE(hms_wat EQ STRING(image_time,FORMAT='(I6.6)'),num_found_image)
       ; print,found_image
       ; print,num_found_image
; Jump image
        IF num_found_image EQ 0 THEN BEGIN
                PRINT,'scan out of range!'
        ENDIF ELSE BEGIN
                go_wat,found_image
        ENDELSE

	END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	TIME_WAT
;
; EXAMPLE:
;
;	Go > time_wat,030000,031000
;
; NOTE:
;
;	Time should be specified in 6-digits (hhmmss).
;

        PRO time_wat,stime,etime ;stime=start time, etime=end time

        COMMON syowa

; Store the data
        stm_wat=stime
        etm_wat=etime

; Verbose
	PRINT,'Start time: '+STRMID(STRING(stime,FORMAT='(I6.6)'),0,4)+' '+ $
		STRMID(STRING(stime,FORMAT='(I6.6)'),4,2)+'s UT'
	PRINT,'End   time: '+STRMID(STRING(etime,FORMAT='(I6.6)'),0,4)+' '+ $
		STRMID(STRING(etime,FORMAT='(I6.6)'),4,2)+'s UT'

; Change time info on Go
        time,stime/100L,etime/100L

        END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       SET_SCALE_WAT
;
; EXAMPLE:
;
;	Go > set_scale_wat,0,300
;

        PRO set_scale_wat,min_val,max_val,quiet=quiet

        COMMON syowa

; Store the data in common block
        min_wat=min_val
        max_wat=max_val

; Verbose
	IF NOT KEYWORD_SET(quiet) THEN BEGIN
		PRINT,'Min intensity: '+STRING(min_val,FORMAT='(I4)')
		PRINT,'Max intensity: '+STRING(max_val,FORMAT='(I4)')
	ENDIF

	END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       PLOT_BOTH
;
; EXAMPLE:
;
;	Go > plot_wat,2,3
;

        ;PRO plot_BOTH,xmaps,ymaps,skip=skip,ctab=ctab,_extra=e

        ;COMMON syowa
        ;COMMON tjornes
	;COMMON ps_info

; Keep current charsize
	;charsize_keep=!P.CHARSIZE

; Refresh window
	;clear_page

; Setting panel positions
        ;IF N_PARAMS() NE 2 THEN BEGIN & xmaps=2 & ymaps=1 & ENDIF
        ;IF NOT KEYWORD_SET(skip) THEN skip=1

; Change charsize
	;!P.CHARSIZE=MAX([MIN([1.0/xmaps,1.0/ymaps]),0.24])
	;IF psinfo.open EQ 0 THEN !P.CHARSIZE=1.5*!P.CHARSIZE
	;IF psinfo.open EQ 1 THEN !P.CHARSIZE=1.1*!P.CHARSIZE



       ; plot_wat_panel,2,1,1,0 ;南極昭和のやつ  ;xmaps,ymaps,xmap,ymap,ctab=ctab,_extra=e     
      ; now_wat+=skip
       ; plot_wat_panel2,2,1,0,0 ;チョルネスのやつ
      ; now_wat2+=skip   

                        
       

; Restore old charsize
	;!P.CHARSIZE=charsize_keep

        ;END

;----------------------------------------------------------------------------------
;
;
;        PLOT_BOTH   ;koretsukattemasu


         PRO plot_both
           COMMON tjornes
           COMMON syowa

           
           plot_wat_panel,2,1,1,0
           plot_wat_panel2,2,1,0,0

      END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       PLOT_WAT
;
; EXAMPLE:
;
;	Go > plot_wat,2,3
;

        PRO plot_wat,xmaps,ymaps,skip=skip,ctab=ctab,_extra=e

        COMMON syowa
        ;COMMON tjornes
	COMMON ps_info

; Keep current charsize
	charsize_keep=!P.CHARSIZE

; Refresh window
	clear_page

; Setting panel positions
        IF N_PARAMS() NE 2 THEN BEGIN & xmaps=1 & ymaps=1 & ENDIF
        IF NOT KEYWORD_SET(skip) THEN skip=1

; Change charsize
	!P.CHARSIZE=MAX([MIN([1.0/xmaps,1.0/ymaps]),0.24])
	IF psinfo.open EQ 0 THEN !P.CHARSIZE=1.5*!P.CHARSIZE
	IF psinfo.open EQ 1 THEN !P.CHARSIZE=1.1*!P.CHARSIZE

; Plot panels
        FOR ymap=0,ymaps-1 DO BEGIN
               FOR xmap=0,xmaps-1 DO BEGIN

                        plot_wat_panel,xmaps,ymaps,xmap,ymap,ctab=ctab,_extra=e     
                        now_wat+=skip
              
               ENDFOR
        ENDFOR
	IF xmaps EQ 1 AND ymaps EQ 1 THEN now_wat-=skip

                        
                        
; Plot colour bar
       plot_wat_colour_bar,ctab=ctab

; Restore old charsize
	!P.CHARSIZE=charsize_keep

        END

















        
;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       PLOT_WAT_PANEL
;
; NOTE:
;
;	1. if ctab is not specified, cut_col_tab is used for plot.
;	2. if cont is specified, the data will be plotted with CONTOUR. Otherwise, plotted with TV.
;	3. if irev is specified, the image is reversed in the east-west direction.
;

PRO plot_wat_panel,xmaps,ymaps,xmap,ymap,irev=irev,position=position,ctab=ctab,cont=cont, $
  myopic_data=myopic_data,myopic_fov=myopic_fov,force_char=force_char,white_grid=white_grid, $
  mnorth=mnorth,no_title=no_title,no_grid=no_grid,special_grid=special_grid,no_ew=no_ew

  COMMON colour_info
  COMMON syowa
	COMMON data
	COMMON prefs

; Update colour info
  IF N_ELEMENTS(ctab) EQ 0 THEN BEGIN
    set_colour_table,/default
  ENDIF ELSE BEGIN
    set_colour_table,ctab
  ENDELSE
  update_colour_info

; Charsize
	IF KEYWORD_SET(force_char) THEN BEGIN
		old_charsize=!P.CHARSIZE
		!P.CHARSIZE=force_char
	ENDIF

; Define plotting region. Default is one panel on screen
  IF N_PARAMS() NE 4 THEN BEGIN & xmaps=1 & ymaps=1 & xmap=0 & ymap=0 & ENDIF
  IF NOT KEYWORD_SET(position) THEN BEGIN
    define_panel,xmaps,ymaps,xmap,ymap,/bar,/square,/no_charset
  	position=!P.POSITION
  ENDIF
  !P.POSITION=position

; Plot frame
  PLOT,[0],[0],XRANGE=[0,siz_wat],YRANGE=[0,siz_wat],/XSTYLE,/YSTYLE, $
          XTICKFORMAT='no_ticks',YTICKFORMAT='no_ticks',XTICKS=1,YTICKS=1

; Plot image with contour
  IF KEYWORD_SET(cont) THEN BEGIN
    c_lines=30
    levels_set=(max_wat-min_wat)*FINDGEN(c_lines)/FLOAT(c_lines)+min_wat
    colour_set=252.0*FINDGEN(c_lines)/FLOAT(c_lines)
    colour_set(0)=0
    x_cont=INDGEN(siz_wat+1)
    y_cont=INDGEN(siz_wat+1)
		IF KEYWORD_SET(irev) THEN BEGIN
      z_cont=REVERSE(REFORM(img_wat(now_wat,*,*)))
		ENDIF ELSE BEGIN
      z_cont=REFORM(img_wat(now_wat,*,*))
		ENDELSE
    CONTOUR,z_cont,x_cont,y_cont,/FILL,LEVELS=levels_set, $
		C_COLORS=colour_set,/OVERPLOT
  ENDIF ELSE BEGIN

; Plot the data with TV
		IF KEYWORD_SET(irev) THEN BEGIN
      img1=REVERSE(REFORM(img_wat(now_wat,*,*)))
		ENDIF ELSE BEGIN
      img1=REFORM(img_wat(now_wat,*,*))
		ENDELSE
    img2=CONGRID(img1,FIX((position(2)-position(0))*!D.X_Size), $
		  FIX((position(3)-position(1))*!D.Y_Size))
    img3=BYTSCL(img2,MAX=max_wat,MIN=min_wat,TOP=252.0)
    TV,img3, position(0)*!D.X_Size,position(1)*!D.Y_Size, XSIZE=(position(2)-position(0))*!D.X_Size, $
			YSIZE=(position(3)-position(1))*!D.Y_Size
  ENDELSE

; Plot info of DATA time
	IF NOT KEYWORD_SET(no_title) THEN BEGIN
                ;XYOUTS,0,siz_wat+5,day_wat(now_wat),/DATA,COL=foreground,CHARSIZE=1.1
        	;XYOUTS,siz_wat-1,siz_wat+5,tim_wat(now_wat),/DATA,COL=foreground,ALI=1,CHARSIZE=1.3
	ENDIF
        PRINT,'Scan'+STRCOMPRESS(now_wat)+' '+day_wat(now_wat)+' '+tim_wat(now_wat)

; Plot the geographic grid (axis is inclined by rog_wat)
	IF NOT KEYWORD_SET(no_grid) THEN BEGIN
		ofs_wid=FLOAT(siz_wat/2)*TAN(rog_wat*!PI/180.0)
		IF KEYWORD_SET(irev) THEN sign=-1 ELSE sign=+1
        	OPLOT,0.5*(siz_wat)*[1,1]+sign*[ofs_wid,-ofs_wid],[0,siz_wat],COL=white,THICK=1,LINE=2
        
	ENDIF

; Plot the geomagnetic S-N (axis is inclined by rog_wat+rom_wat)
	ofs_wid=FLOAT(siz_wat/2)*TAN((rom_wat)*!PI/180.0)
	IF KEYWORD_SET(irev) THEN sign=-1 ELSE sign=+1
       OPLOT,0.5*(siz_wat)*[1,1]+sign*[ofs_wid,-ofs_wid],[0,siz_wat],COL=white,THICK=1,LINE=1
       OPLOT,[0,siz_wat],0.5*(siz_wat)*[1,1]-sign*[ofs_wid,-ofs_wid],COL=white,THICK=1,LINE=1
; Plot direction
	IF NOT KEYWORD_SET(no_grid) THEN BEGIN
		IF NOT KEYWORD_SET(irev) THEN BEGIN
			XYOUTS,siz_wat+5,siz_wat/2+(siz_wat/2)*TAN(rog_wat*!DTOR)-135,'N',ALI=0.0,/DATA
			XYOUTS,-5,siz_wat/2-(siz_wat/2)*TAN(rog_wat*!DTOR)+125,'S',ALI=1.0,/DATA
			XYOUTS,siz_wat/2-(siz_wat/2)*TAN(rog_wat*!DTOR),siz_wat+5,'MS',ALI=0.5,/DATA
			XYOUTS,siz_wat/2+(siz_wat/2)*TAN(rog_wat*!DTOR)-125,-15,'E',ALI=0.5,/DATA
			XYOUTS,siz_wat/2-(siz_wat/2)*TAN((rog_wat+rom_wat)*!DTOR)+70,siz_wat+5,'W',ALI=0.5,/DATA
		ENDIF ELSE BEGIN
			XYOUTS,-150,siz_wat/2+(siz_wat/2)*TAN(rog_wat*!DTOR)-5,'E',ALI=1.0,/DATA
			XYOUTS,siz_wat+5,siz_wat/2-(siz_wat/2)*TAN(rog_wat*!DTOR)-5,'W',ALI=0.0,/DATA
			;XYOUTS,siz_wat/2+(siz_wat/2)*TAN(rog_wat*!DTOR),siz_wat+5,'MS',ALI=0.5,/DATA
		;	XYOUTS,siz_wat/2-(siz_wat/2)*TAN(rog_wat*!DTOR),-15,'N',ALI=0.5,/DATA
			;XYOUTS,siz_wat/2+(siz_wat/2)*TAN((rog_wat+rom_wat)*!DTOR),siz_wat+5,'MS',ALI=0.5,/DATA
		ENDELSE
	ENDIF

; Mask edge and overlay zenith
	mask_fov
	IF NOT KEYWORD_SET(no_grid) THEN overlay_zenith
	IF KEYWORD_SET(special_grid) THEN overlay_zenith_special

; Overplot SuperDARN stuff?
	IF KEYWORD_SET(myopic_fov) THEN overlay_myopic_fov,irev=irev
	IF KEYWORD_SET(myopic_data) THEN overlay_myopic_data,irev=irev

; Overplot frame again
        PLOT,[0],[0],XRANGE=[0,siz_wat],YRANGE=[0,siz_wat],/XSTYLE,/YSTYLE, $
                XTICKFORMAT='no_ticks',YTICKFORMAT='no_ticks',XTICKS=1,YTICKS=1
	IF NOT KEYWORD_SET(no_title) THEN BEGIN
	;	XYOUTS,10,10,'Outermost circle is 75!9'+STRING("260B")+'!3 zenith angle',CHARSIZE=0.6*!P.CHARSIZE,COL=white
		IF KEYWORD_SET(irev) THEN $
			XYOUTS,siz_wat+3,siz_wat,'East-West Reversed',ORI=270,CHARSIZE=0.5*!P.CHARSIZE
	ENDIF

; Plot SuperDARN scan time
	IF KEYWORD_SET(myopic_data) THEN BEGIN
		scan_beams=WHERE(beam_scan EQ scan_no,no_scan_beams)
		pos=!P.POSITION
		XYOUTS,pos(0)+0.01,pos(1)+0.02,str_time(beam_time(scan_beams(0)),1)+' UT',/NORMAL,CHARSIZE=1.4
	ENDIF

; Restore old charsize
	IF KEYWORD_SET(force_char) THEN !P.CHARSIZE=old_charsize

	END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       OVERLAY_WAT_CONT
;

        PRO overlay_wat_cont,irev=irev

        COMMON colour_info
        COMMON syowa

	COMMON data
	COMMON prefs

; Update colour info
        IF N_ELEMENTS(ctab) EQ 0 THEN BEGIN
                set_colour_table,/default
        ENDIF ELSE BEGIN
                set_colour_table,ctab
        ENDELSE
        update_colour_info

; Plot image with contour
        ;c_lines=30
        ;levels_set=(max_wat-min_wat)*FINDGEN(c_lines)/FLOAT(c_lines)+min_wat
        ;colour_set=252.0*FINDGEN(c_lines)/FLOAT(c_lines)
        ;colour_set(0)=0
	levels_set=[80,160]
        x_cont=INDGEN(siz_wat+1)
        y_cont=INDGEN(siz_wat+1)
	IF KEYWORD_SET(irev) THEN BEGIN
      		z_cont=REVERSE(REFORM(img_wat(now_wat,*,*)))
	ENDIF ELSE BEGIN
        	z_cont=REFORM(img_wat(now_wat,*,*))
	ENDELSE
        CONTOUR,z_cont,x_cont,y_cont,LEVELS=levels_set,/OVERPLOT,COL=0,THICK=6
        CONTOUR,z_cont,x_cont,y_cont,LEVELS=levels_set,/OVERPLOT,COL=254,THICK=2

	END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	PLOT_WAT_KEO
;

	PRO plot_wat_keo,_extra=e

	COMMON syowa
	COMMON ps_info

; Refresh window
	clear_page

; Check time info
        time_info=STRING(stm_wat,FORMAT='(I6.6)')+' to '+STRING(etm_wat,FORMAT='(I6.6)')+' UT'
        print,time_info

; Keep current charsize
	charsize_keep=!P.CHARSIZE

; Change charsize
	IF psinfo.open EQ 0 THEN !P.CHARSIZE=1.0*!P.CHARSIZE ;タイトルとかの文字サイズ 
	IF psinfo.open EQ 1 THEN !P.CHARSIZE=0.7*!P.CHARSIZE

; Plot two keograms
	plot_wat_keo_panel,1,2,0,0,line_id=0,_extra=e
	XYOUTS,!P.POSITION(0),!P.POSITION(3)+0.02,'Keogram from MN to MS in SYOWA', $
		/NORMAL
	XYOUTS,!P.POSITION(2),!P.POSITION(3)+0.01,time_info,/NORMAL,ALI=1.0
	plot_wat_keo_panel2,1,2,0,1,line_id=2,_extra=e
	XYOUTS,!P.POSITION(0),!P.POSITION(3)+0.01,'Keogram from MS to MN in TJORNES', $
		/NORMAL
	XYOUTS,!P.POSITION(2),!P.POSITION(3)+0.01,time_info,/NORMAL,ALI=1.0
;	plot_wat_keo_panel,1,3,0,2,line_id=4,_extra=e
;	XYOUTS,!P.POSITION(0),!P.POSITION(3)+0.01,'Kishiyamas Keogram',$;'Beam 7 aligned Keogram', $
;		/NORMAL
;	XYOUTS,!P.POSITION(2),!P.POSITION(3)+0.01,time_info,/NORMAL,ALI=1.0

; Plot colour bar
	plot_wat_colour_bar,_extra=e

; Restore old charsize
	!P.CHARSIZE=charsize_keep

     END
;------------------------------------------------------------------------------------------------------------------
; NAME: 
;
;      TIME_BOTH_WAT
;

       PRO time_both_wat,time,time2

         COMMON syowa
         COMMON tjornes
         COMMON ps_info


         time_wat,time,time2
         time_wat2,time,time2


         END
;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	SET_SCALE_BOTH
;    
        PRO set_scale_both,sc1,sc2

          COMMON syowa
          COMMON tjornes


          set_scale_Wat,sc1,sc2
          set_scale_wat2,sc1,sc2

       END
        
  

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	PLOT_WAT_KEO_BOTH
;

	PRO plot_wat_keo_both,time,time2,sc1=sc1,sc2=sc2,_extra=e

        COMMON syowa
        COMMON tjornes  		

	COMMON ps_info

; Refresh window
	clear_page


        time_both_wat,time,time2
        set_scale_both,sc1,sc2
; Check time info
        time_info=STRING(stm_wat,FORMAT='(I6.6)')+' to '+STRING(etm_wat,FORMAT='(I6.6)')+' UT'
        print,time_info

; Keep current charsize
	charsize_keep=!P.CHARSIZE

; Change charsize
	IF psinfo.open EQ 0 THEN !P.CHARSIZE=0.9*!P.CHARSIZE   ;タイトルとかの文字サイズ 
	IF psinfo.open EQ 1 THEN !P.CHARSIZE=0.9*!P.CHARSIZE

; Plot two keograms
	plot_wat_keo_panel,1,2,0,1,line_id=0,_extra=e
	XYOUTS,!P.POSITION(0),!P.POSITION(3)+0.01,'Keogram from MN to MS in SYOWA', $
		/NORMAL
	XYOUTS,!P.POSITION(2),!P.POSITION(3)+0.01,time_info,/NORMAL,ALI=1.0
	plot_wat_keo_panel2,1,2,0,0,line_id=2,_extra=e
	XYOUTS,!P.POSITION(0),!P.POSITION(3)+0.01,'Keogram from MS to MN in TJORNES', $
		/NORMAL
	XYOUTS,!P.POSITION(2),!P.POSITION(3)+0.01,time_info,/NORMAL,ALI=1.0

; Plot colour bar
	plot_wat_colour_bar,_extra=e

; Restore old charsize
	!P.CHARSIZE=charsize_keep

	END


        
;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       PLOT_WAT_KEO_PANEL
;
; EXAMPLE:
;
;	Keogram from S to N cut
;	Go > plot_wat_keo_panel,1,3,0,0,line_id=0
;

        PRO plot_wat_keo_panel,xmaps,ymaps,xmap,ymap,position=position, $
                zan=zan,azm=azm,line_id=line_id,ctab=ctab,no_x=no_x,no_y=no_y,no_data=no_data

        COMMON syowa

        COMMON colour_info
        update_colour_info

; Line id (line_id=0; N-S cut is default)
        IF NOT KEYWORD_SET(line_id) THEN line_id=0
	IF line_id GT 4 THEN BEGIN
		PRINT,'line_id must be'
		PRINT,'0: GS to GN'
		PRINT,'1: GE to GW'
		PRINT,'2: MS to MN'
		PRINT,'3: ME to MW'
		PRINT,'4: Beam 7 aligned'
		RETURN
	ENDIF

; Tick mark setting
        IF NOT KEYWORD_SET(no_x) THEN no_x=0

; Update colour info
        IF N_ELEMENTS(ctab) EQ 0 THEN BEGIN
                set_colour_table,/default
        ENDIF ELSE BEGIN
                set_colour_table,ctab
        ENDELSE
        update_colour_info

; Define plotting region. Default is one panel on screen
        IF N_PARAMS() NE 4 THEN BEGIN & xmaps=1 & ymaps=1 & xmap=0 & ymap=0 & ENDIF
        IF NOT KEYWORD_SET(position) THEN BEGIN
                define_panel,xmaps,ymaps,xmap,ymap,/BAR,/no_charset
                position=!P.POSITION
        ENDIF

; Calculate real data area along the cross-section
	IF line_id EQ 0 OR line_id EQ 1 THEN BEGIN
		real_wid=FIX((siz_wat/2)*COS(rog_wat*!DTOR))  ;!DTOR = PI/180
		img_rng=[siz_wat/2-real_wid,siz_wat/2+real_wid]
	ENDIF
	IF line_id EQ 2 OR line_id EQ 3 THEN BEGIN
		real_wid=FIX((siz_wat/2)*COS((rog_wat+rom_wat)*!DTOR))
		img_rng=[siz_wat/2-real_wid,siz_wat/2+real_wid]
	ENDIF
	IF line_id EQ 4 THEN BEGIN
		real_wid=FIX((siz_wat/2)*COS((rog_wat+32.0-3.24*0.5)*!DTOR))
		img_rng=[siz_wat/2-real_wid,siz_wat/2+real_wid]
             ENDIF

; Plot image with TV
	IF NOT KEYWORD_SET(no_data) THEN BEGIN
                plot_img=WHERE(hms_wat GE stm_wat AND hms_wat LE etm_wat,num_plot_img)
               ; plot_img=WHERE(stm_wat AND etm_wat,num_plot_img)
                IF num_plot_img EQ 0 THEN RETURN
        	img1=REFORM(keo_wat(line_id,plot_img,img_rng(0):img_rng(1)))
        	img2=CONGRID(img1,FIX((position(2)-position(0))*!D.X_Size),FIX((position(3)-position(1))*!D.Y_Size))
        	img3=BYTSCL(img2,MAX=max_wat,MIN=min_wat,TOP=252)
        	TV,img3,FIX(position(0)*!D.X_Size),FIX(position(1)*!D.Y_Size), $
                	XSIZE=FIX((position(2)-position(0))*!D.X_Size),YSIZE=FIX((position(3)-position(1))*!D.Y_Size)
	ENDIF

; Plot frame
	IF line_id EQ 0 THEN ytitle='Zenith Angle (deg)!C!Dfrom South to North!N'
	IF line_id EQ 1 THEN ytitle='Zenith Angle (deg)!C!Dfrom West to East!N'
	IF line_id EQ 2 THEN ytitle='Zenith Angle (deg)!C!Dfrom Mag South to Mag North!N'
	IF line_id EQ 3 THEN ytitle='Zenith Angle (deg)!C!Dfrom Mag West to Mag East!N'
	IF line_id EQ 4 THEN ytitle='Zenith Angle (deg)!C!Daligned with Beam 7!N'
        plot_time_frame,position=position,yrange=[-75,75],no_x=no_x,ytitle=ytitle,yticks=6,yminor=5

; Plot zenith and CP1 fov
        OPLOT,[0,25],[0,0],COL=background,LINE=0
        OPLOT,[0,25],[0,0],COL=foreground,LINE=2

        BLACK='FFFFFF'
        IF KEYWORD_SET(zan) THEN BEGIN
           IF azm EQ  0  THEN BEGIN
           zan_p=zan
           OPLOT,[0,25],[zan_p,zan_p],COL=BLACK,LINE=0,THICK=4
           ENDIF
           IF azm EQ 180 THEN BEGIN
           zan_m=-zan
           OPLOT,[0,25],[zan_m,zan_m],COL=BLACK,LINE=0,THICK=4
           ENDIF
        ENDIF
        

     END


;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       PLOT_WAT_KEO_PANEL2
;
; EXAMPLE:
;
;	Keogram from S to N cut
;	Go > plot_wat_keo_panel,1,3,0,0,line_id=0
;

        PRO plot_wat_keo_panel2,xmaps,ymaps,xmap,ymap,position=position, $
                zan=zan,azm=azm,line_id=line_id,ctab=ctab,no_x=no_x,no_y=no_y,no_data=no_data

        COMMON tjornes

        COMMON colour_info
        update_colour_info

; Line id (line_id=0; N-S cut is default)
        IF NOT KEYWORD_SET(line_id) THEN line_id=0
	IF line_id GT 4 THEN BEGIN
		PRINT,'line_id must be'
		PRINT,'0: GS to GN'
		PRINT,'1: GE to GW'
		PRINT,'2: MS to MN'
		PRINT,'3: ME to MW'
		PRINT,'4: Beam 7 aligned'
		RETURN
	ENDIF

; Tick mark setting
        IF NOT KEYWORD_SET(no_x) THEN no_x=0

; Update colour info
        IF N_ELEMENTS(ctab) EQ 0 THEN BEGIN
                set_colour_table,/default
        ENDIF ELSE BEGIN
                set_colour_table,ctab
        ENDELSE
        update_colour_info

; Define plotting region. Default is one panel on screen
        IF N_PARAMS() NE 4 THEN BEGIN & xmaps=1 & ymaps=1 & xmap=0 & ymap=0 & ENDIF
        IF NOT KEYWORD_SET(position) THEN BEGIN
                define_panel,xmaps,ymaps,xmap,ymap,/BAR,/no_charset
                position=!P.POSITION
        ENDIF

; Calculate real data area along the cross-section
	IF line_id EQ 0 OR line_id EQ 1 THEN BEGIN
		real_wid=FIX((siz_wat2/2)*COS(rog_wat2*!DTOR))
		img_rng=[siz_wat2/2-real_wid,siz_wat2/2+real_wid]
	ENDIF
	IF line_id EQ 2 OR line_id EQ 3 THEN BEGIN
		real_wid=FIX((siz_wat2/2)*COS((rog_wat2+rom_wat2)*!DTOR))
		img_rng=[siz_wat2/2-real_wid,siz_wat2/2+real_wid]
	ENDIF
	IF line_id EQ 4 THEN BEGIN
		real_wid=FIX((siz_wat2/2)*COS((rog_wat2+32.0-3.24*0.5)*!DTOR))
		img_rng=[siz_wat2/2-real_wid,siz_wat2/2+real_wid]
             ENDIF

; Plot image with TV
	IF NOT KEYWORD_SET(no_data) THEN BEGIN
        	plot_img=WHERE(hms_wat2 GE stm_wat2 AND hms_wat2 LE etm_wat2,num_plot_img)
        	IF num_plot_img EQ 0 THEN RETURN
        	img1=REFORM(keo_wat2(line_id,plot_img,img_rng(0):img_rng(1)))
        	img2=CONGRID(img1,FIX((position(2)-position(0))*!D.X_Size),FIX((position(3)-position(1))*!D.Y_Size))
        	img3=BYTSCL(img2,MAX=max_wat2,MIN=min_wat2,TOP=252)
        	TV,img3,FIX(position(0)*!D.X_Size),FIX(position(1)*!D.Y_Size), $
                	XSIZE=FIX((position(2)-position(0))*!D.X_Size),YSIZE=FIX((position(3)-position(1))*!D.Y_Size)
	ENDIF

; Plot frame
	IF line_id EQ 0 THEN ytitle='Zenith Angle (deg)!C!Dfrom South to North!N'
	IF line_id EQ 1 THEN ytitle='Zenith Angle (deg)!C!Dfrom West to East!N'
	IF line_id EQ 2 THEN ytitle='Zenith Angle (deg)!C!Dfrom Mag South to Mag North!N'
	IF line_id EQ 3 THEN ytitle='Zenith Angle (deg)!C!Dfrom Mag West to Mag East!N'
	IF line_id EQ 4 THEN ytitle='Zenith Angle (deg)!C!Daligned with Beam 7!N'
        plot_time_frame,position=position+[0,0,0,0],yrange=[-75,75],no_x=no_x,ytitle=ytitle,yticks=6,yminor=5

; Plot zenith and CP1 fov
        OPLOT,[0,25],[0,0],COL=background,LINE=0
        OPLOT,[0,25],[0,0],COL=foreground,LINE=2

        BLACK='FFFFFF'
        IF KEYWORD_SET(zan) THEN BEGIN
           IF azm EQ  0  THEN BEGIN
           zan_p=zan
           OPLOT,[0,25],[zan_p,zan_p],COL=BLACK,LINE=0,THICK=4
           ENDIF
           IF azm EQ 180 THEN BEGIN
           zan_m=-zan
           OPLOT,[0,25],[zan_m,zan_m],COL=BLACK,LINE=0,THICK=4
           ENDIF
        ENDIF

        
	END


        
;-------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       PLOT_WAT_RTI_PANEL
;

        PRO plot_wat_rti_panel,yrange=yrange,no_x=no_x,no_frame=no_frame,beam_plot=beam_plot,cont=cont

        COMMON syowa

; Variables
        radar_lat=63.77
        radar_lon=-20.54
        radar_boresite=30.0
        gate_length=15.0
        first_range=180.0
        altitude=110.0

; Angle of Tjornes Watec alignment
        rot_ang=rog_wat

; Size of the image
        wid=siz_wat/2

; Find radar cell locations
        x_wat=FLTARR(16,75)
        y_wat=FLTARR(16,75)
        FOR beam=0,15 DO BEGIN
                FOR gate=0,74 DO BEGIN
                        fov_pos=find_cell(radar_lat,radar_lon,16,radar_boresite, $
                                3.24,gate_length,first_range,beam,gate+0.5)
                        lat=fov_pos(0)
                        lon=fov_pos(1)
                        zen=get_zenith_angle_from_tjornes(lat,lon,altitude)
                        ;azm=get_azimuth_from_tjornes(lat,lon)+90.0-rot_ang
                        ;azm=get_azimuth_from_tjornes(lat,lon)+90.0
                        azm=get_azimuth_from_tjornes(lat,lon)
                        IF azm GE 360.0 THEN azm=azm-360.0
                        tmp_pix=get_pixel(zen,azm)
                        x_wat(beam,gate)=tmp_pix(0)
                        y_wat(beam,gate)=tmp_pix(1)
                        ;x_wat(beam,gate)=wid+zen*COS(!PI*azm/180.0)*wid/75.0
                        ;y_wat(beam,gate)=wid+zen*SIN(!PI*azm/180.0)*wid/75.0
                ENDFOR
        ENDFOR

; Extract optical data from 2D array
        rti_wat=INTARR(num_wat,50)
        IF NOT KEYWORD_SET(beam_plot) THEN beam_plot=7
        FOR i=0,num_wat-1 DO BEGIN

                FOR gate=0,49 DO BEGIN
                        IF FIX(x_wat(beam_plot,gate)) LT 367 AND FIX(y_wat(beam_plot,gate)) LT 367 THEN BEGIN

                                rti_wat(i,gate)=img_wat(i,FIX(x_wat(beam_plot,gate)),FIX(y_wat(beam_plot,gate)))
                        ENDIF
                ENDFOR

        ENDFOR

; Plot panel
        IF NOT KEYWORD_SET(yrange) THEN yrange=[10,40]
        IF NOT KEYWORD_SET(no_frame) THEN BEGIN
                plot_time_frame,position=!P.POSITION,yrange=yrange,ytitle='Range Gate',no_x=no_x
        ENDIF

; Plot the data
	IF KEYWORD_SET(cont) THEN BEGIN
		x_cont=tsc_wat(0:num_wat-1)/3600.0
		y_cont=FINDGEN(50)+11.0
        	;contour_lines=5
        	;levels_set=(max_wat-min_wat)*FINDGEN(contour_lines)/FLOAT(contour_lines)+min_wat
        	;colour_set=!D.TABLE_SIZE*INDGEN(contour_lines)/contour_lines
		levels_set=[100,200]
                CONTOUR,rti_wat,x_cont,y_cont,/OVERPLOT,LEVELS=levels_set,NOCLIP=0,COL=254,THICK=6
                CONTOUR,rti_wat,x_cont,y_cont,/OVERPLOT,LEVELS=levels_set,NOCLIP=0,COL=0,THICK=3
                ;CONTOUR,rti_wat,x_cont,y_cont,/OVERPLOT,/FILL,LEVELS=levels_set, $
                ;       	C_COLORS=colour_set,NOCLIP=0
	ENDIF ELSE BEGIN
        	FOR i=0,num_wat-2 DO BEGIN
                	FOR j=0,49 DO BEGIN

                	x_poly=[tsc_wat(i),tsc_wat(i+1),tsc_wat(i+1),tsc_wat(i)]/3600.0
                	y_poly=[j,j,j+1,j+1]+11

                	col=252.0*(rti_wat(i,j)-min_wat)/FLOAT(max_wat-min_wat)
                	IF col LE   0 THEN col=1
                	IF col GE 252 THEN col=252
                	POLYFILL,x_poly,y_poly,COL=col,NOCLIP=0
	
                	ENDFOR
        	ENDFOR
	ENDELSE

; Overplot frame
        IF NOT KEYWORD_SET(no_frame) THEN $
                plot_time_frame,position=!P.POSITION,yrange=yrange,ytitle='Range Gate',no_x=no_x

        END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       PLOT_WAT_POLAR
;
; EXAMPLE:
;
;       Go > plot_wat_polar,2,3
;

        PRO plot_wat_polar,xmaps,ymaps,skip=skip,ctab=ctab,_extra=e

        COMMON syowa
        COMMON ps_info

; Keep current charsize
        charsize_keep=!P.CHARSIZE

; Refresh window
        clear_page

; Setting panel positions
        IF N_PARAMS() NE 2 THEN BEGIN & xmaps=1 & ymaps=1 & ENDIF
        IF NOT KEYWORD_SET(skip) THEN skip=1

; Change charsize
        !P.CHARSIZE=MAX([MIN([1.0/xmaps,1.0/ymaps]),0.24])
        IF psinfo.open EQ 0 THEN !P.CHARSIZE=2.0*!P.CHARSIZE
        IF psinfo.open EQ 1 THEN !P.CHARSIZE=1.5*!P.CHARSIZE

; Plot panels
        FOR ymap=0,ymaps-1 DO BEGIN
                FOR xmap=0,xmaps-1 DO BEGIN

                        plot_wat_polar_panel,xmaps,ymaps,xmap,ymap,ctab=ctab,_extra=e
                        now_wat+=skip
                ENDFOR
        ENDFOR
        IF xmaps EQ 1 AND ymaps EQ 1 THEN now_wat-=skip

; Plot colour bar
        plot_wat_colour_bar,ctab=ctab

; Restore old charsize
        !P.CHARSIZE=charsize_keep

        END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	PLOT_WAT_POLAR_PANEL
;
	
	PRO plot_wat_polar_panel,xmaps,ymaps,xmap,ymap,xrange=xrange,yrange=yrange,clip=clip, $
		position=position,coast=coast,no_data=no_data,ctab=ctab,force_char=force_char, $
		no_title=no_title

	COMMON syowa

        COMMON data
        COMMON prefs

; Lat and lon of Tjornes
	lat_tjo=66.20
	lon_tjo=342.88

; Update colour info
        IF N_ELEMENTS(ctab) EQ 0 THEN BEGIN
                set_colour_table,/default
        ENDIF ELSE BEGIN
                set_colour_table,ctab
        ENDELSE
        update_colour_info

; Charsize
        IF KEYWORD_SET(force_char) THEN BEGIN
                old_charsize=!P.CHARSIZE
                !P.CHARSIZE=force_char
        ENDIF

; Define plotting region. Default is one panel on screen
        IF N_PARAMS() NE 4 THEN BEGIN & xmaps=1 & ymaps=1 & xmap=0 & ymap=0 & ENDIF
        IF NOT KEYWORD_SET(position) THEN BEGIN
                define_panel,xmaps,ymaps,xmap,ymap,/bar,/square,/no_charset
                position=!P.POSITION
        ENDIF
        !P.POSITION=position

; Define plot area
        IF NOT KEYWORD_SET(xrange) OR NOT KEYWORD_SET(yrange) THEN BEGIN
                xrange=[-31,31] & yrange=[-31,31]
                IF KEYWORD_SET(clip) THEN BEGIN
                        pos=cnvcoord(lat_tjo,lon_tjo,1)
                        mla_tjo=pos(0) & mlo_tjo=pos(1)
                        x1= abs(90-mla_tjo)*SIN(mlt(yer_wat(now_wat),yrs_wat(now_wat),mlo_tjo)*!PI/12)
                        y1=-abs(90-mla_tjo)*COS(mlt(yer_wat(now_wat),yrs_wat(now_wat),mlo_tjo)*!Pi/12)
			xyrange=10.0
                        xrange=[x1-xyrange*0.5,x1+xyrange*0.5]
                        yrange=[y1-xyrange*0.5,y1+xyrange*0.5]
                ENDIF
        ENDIF

; Plot frame
        plot_polar_frame,position=!P.POSITION, $
                xrange=xrange,yrange=yrange,black_back=black_back,ticklen=ticklen,no_frame=no_frame

; Plot image as contours
        tmp_bin=WHERE(gla_wat NE 99999.9 AND gla_wat NE 99999.9 AND zan_wat LE 75.0,no_tmp_bin)
        cont=MAKE_ARRAY(4,no_tmp_bin,/FLOAT,VALUE=0.0)
        num_cnt=0L
        FOR xxx=0,siz_wat DO BEGIN
                FOR yyy=0,siz_wat DO BEGIN
                        IF gla_wat(xxx,yyy) NE 99999.9 AND glo_wat(xxx,yyy) NE 99999.9 AND $
                                zan_wat(xxx,yyy) LE 75.0 THEN BEGIN

                                cont(2,num_cnt)=FLOAT(img_wat(now_wat,xxx,yyy))
                                cont(3,num_cnt)=zan_wat(xxx,yyy)
                                IF cont(2,num_cnt) GT max_wat THEN cont(2,num_cnt)=max_wat
                                IF cont(2,num_cnt) LT min_wat THEN cont(2,num_cnt)=min_wat

                                cont(0,num_cnt)= ABS(90-mla_wat(xxx,yyy)) $
                                        *SIN(mlt(yer_wat(now_wat),yrs_wat(now_wat),mlo_wat(xxx,yyy))*!PI/12)
                                cont(1,num_cnt)=-ABS(90-mla_wat(xxx,yyy)) $
                                        *COS(mlt(yer_wat(now_wat),yrs_wat(now_wat),mlo_wat(xxx,yyy))*!PI/12)

                                num_cnt++
                        ENDIF
                ENDFOR
        ENDFOR
        contour_lines=30
        levels_set=(max_wat-min_wat)*FINDGEN(contour_lines)/FLOAT(contour_lines)+min_wat
        levels_lin=(max_wat-min_wat)*FINDGEN(6)/FLOAT(5)+min_wat
        colour_set=!D.TABLE_SIZE*INDGEN(contour_lines)/contour_lines
        IF NOT KEYWORD_SET(no_data) THEN BEGIN
                CONTOUR,cont(2,*),cont(0,*),cont(1,*),/OVERPLOT,/FILL,LEVELS=levels_set, $
                        C_COLORS=colour_set,/IRREGULAR,NOCLIP=0
        ENDIF

; Plot frame
        plot_polar_frame,position=!P.POSITION, $
                xrange=xrange,yrange=yrange,black_back=black_back,ticklen=ticklen,no_frame=no_frame

; Overlay coast
        IF KEYWORD_SET(coast) THEN $
                overlay_polar_coast,force_year=yer_wat(now_wat),force_secs=yrs_wat(now_wat),col=253

; Plot info
        IF NOT KEYWORD_SET(no_title) THEN BEGIN
		pos=!P.POSITION
                XYOUTS,pos(0),pos(3)+0.01,day_wat(now_wat),/NORMAL,COL=foreground
                XYOUTS,pos(2),pos(3)+0.01,tim_wat(now_wat),/NORMAL,COL=foreground,ALI=1
        ENDIF
        PRINT,'Scan'+STRCOMPRESS(now_wat)+' '+day_wat(now_wat)+' '+tim_wat(now_wat)

; Restore old charsize
        IF KEYWORD_SET(force_char) THEN !P.CHARSIZE=old_charsize

	END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	PLOT_WAT_LINE
;

	PRO plot_wat_line,time,time2,zan=zan,azm=azm,_extra=e,click=click,n_smo=n_smo,b_ps=b_ps;,posts=posts

        COMMON TJORNES
        COMMON SYOWA
; Refresh page if click keyword is not specified
	IF NOT KEYWORD_SET(click) THEN BEGIN
		click=0
		clear_page
	ENDIF
        set_scale_both,0,25
        set_scale_wat2,0,30
        time_both_wat,time,time2
        
; Plot only one panel
        plot_wat_line_panel2,1,2,0,0,_extra=e,click=click,zan=zan,azm=azm,n_smo=n_smo,b_ps=b_ps
        plot_wat_line_panel,1,2,0,1,_extra=e,click=click,zan=zan,azm=azm,n_smo=n_smo,b_ps=b_ps   
	END
        
;------------------------------------------------------------------------------------------------------------------
;
; NAME:
;
;	PLOT_WAT_LINE_PANEL2
;

        PRO plot_wat_line_panel2,xmaps,ymaps,xmap,ymap,position=position, $
                colour=colour,no_x=no_x,no_y=no_y,yrange=yrange,zan=zan,azm=azm, $
                right=right,thick=thick,click=click,ytitle=ytitle,yticks=yticks,n_smo=n_smo,plus=plus,b_ps=b_ps

        COMMON tjornes
        COMMON syowa
        COMMON prefs
        COMMON colour_info

; Setting axis
        IF NOT KEYWORD_SET(no_x) THEN no_x=0
        IF NOT KEYWORD_SET(no_y) THEN no_y=0
        IF NOT KEYWORD_SET(yrange) THEN yrange=[min_wat2,max_wat2]

; Define plotting region. Default is one panel on screen
        IF N_PARAMS() NE 4 THEN BEGIN & xmaps=1 & ymaps=1 & xmap=0 & ymap=0 & ENDIF
        IF NOT KEYWORD_SET(position) THEN BEGIN
                define_panel,xmaps,ymaps,xmap,ymap,/bar,/no_charset
                position=!P.POSITION
        ENDIF
        IF KEYWORD_SET(right) THEN ystyle=5 ELSE ystyle=0

; Some settings
	IF NOT KEYWORD_SET(click) THEN BEGIN
        	IF NOT KEYWORD_SET(colour) THEN colour=foreground
		IF NOT KEYWORD_SET(zan) THEN zan=0
		IF NOT KEYWORD_SET(azm) THEN azm=0
		pnt=get_pixel2(zan,azm)
        	cx=pnt(0)    ; calculated with get_pixel
                cy=pnt(1)    ; calculated with get_pixel
             
	ENDIF ELSE BEGIN
		CURSOR,cx,cy,/wait
		cx=FIX(cx)
		cy=FIX(cy)
                ang=get_zan_azm(cx,cy)
               ; print,ang(0)
                ang(0)=50
               ; print,ang(0)
                zan=FIX(ang(0))
               ; print,zan
		azm=FIX(ang(1))
		clear_page
             ENDELSE
        
        print,zan
	PRINT,'Zenith  angle: '+STRING(zan,FORMAT='(F4.1)')
	PRINT,'Azimuth angle: '+STRING(azm,FORMAT='(F5.1)')

; Plot panel
	IF NOT KEYWORD_SET(ytitle) THEN ytitle='Auroral Intensity (Arbitrary Unit)'
        plot_time_frame,position=position+[0,0.05,0,0.05],ytitle=ytitle, $
                no_x=no_x,yrange=yrange,ystyle=ystyle,no_y=no_y,yticks=yticks

; Data for plot
       
        plot1_wat=img_wat2(*,cx,cy)
        plot2_wat=SMOOTH(plot1_wat,3)
        plot_wat2=plot2_wat-SMOOTH(plot2_wat,30)+10
       
       ; plot_wat2=SMOOTH(plot2_wat,30)+12
; 

; Plot the data


        
        
        IF NOT KEYWORD_SET(thick) THEN thick=!P.THICK
	x_plot2=tsc_wat2/3600.0
	IF stm_wat2 GE 240000 THEN x_plot2=x_plot2-24.0

        
        IF NOT KEYWORD_SET(b_ps) THEN b_ps=0     
        IF b_ps EQ 0 THEN BEGIN
           color=255
           IF n_smo EQ 1 THEN BEGIN
              OPLOT,x_plot2,plot1_wat,COLOR=color,THICK=thick
           ENDIF
           IF NOT KEYWORD_SET(n_smo) THEN n_smo=0
           IF n_smo EQ 0 THEN BEGIN
              OPLOT,x_plot2,plot_wat2,COLOR=color,THICK=thick
           ENDIF
        ENDIF

       
        
        IF b_ps EQ 1 THEN BEGIN
           color=0
           IF n_smo EQ 1 THEN BEGIN
              OPLOT,x_plot2,plot1_wat,COLOR=color,THICK=thick
           ENDIF
           IF NOT KEYWORD_SET(n_smo) THEN n_smo=0
           IF n_smo EQ 0 THEN BEGIN
              OPLOT,x_plot2,plot_wat2,COLOR=color,THICK=thick
           ENDIF
        ENDIF
        
        color=240
        ;cx1_wat=DBLARR(100)
        FOR i=1,100-1 DO  BEGIN
           ;CURSOR,cx1,cy1
           ;OPLOT,[cx1_wat(i),cx1_wat(i)],[0,30],color=color,LINESTYLE=0,THICK=3
           ;cx1_wat(i)=cx1     
        ENDFOR


        

; Make zenith angle axis
        IF KEYWORD_SET(right) THEN BEGIN
                AXIS,yaxis=1,yrange=yrange,YSTYLE=1,TICKLEN=-!P.TICKLEN, $
                        ytitle='Auroral Intensity (Arbitrary Unit)',COLOR=colour, $
			CHARTHICK=!P.CHARTHICK
        ENDIF

; Check time info
        

        time_info=STRING(stm_wat2,FORMAT='(I6.6)')+' to '+ $
		STRING(etm_wat2,FORMAT='(I6.6)')+' UT'
	XYOUTS,!P.POSITION(2),!P.POSITION(3)+0.01,time_info,/NORMAL,ALI=1.0
	pixel_info='TJORNES Zenith:'+STRING(zan,FORMAT='(F4.1)')+' deg, '+ $
		'Azm:'+STRING(azm,FORMAT='(F5.1)')+' deg'
	XYOUTS,!P.POSITION(0),!P.POSITION(3)+0.01,pixel_info,/NORMAL,ALI=0.0

       

       
    END


        
;------------------------------------------------------------------------------------------------------------------
;
; NAME:
;
;	PLOT_WAT_LINE_PANEL
;

        PRO plot_wat_line_panel,xmaps,ymaps,xmap,ymap,position=position, $
                  colour=colour,no_x=no_x,no_y=no_y,yrange=yrange,zan=zan,azm=azm, $
                  inf=inf,right=right,thick=thick,click=click,$
                  ytitle=ytitle,yticks=yticks,n_smo=n_smo,b_ps=b_ps

        COMMON syowa
        COMMON tjornes  
        COMMON fft
        COMMON prefs
        COMMON colour_info

; Setting axis
        IF NOT KEYWORD_SET(no_x) THEN no_x=0
        IF NOT KEYWORD_SET(no_y) THEN no_y=0
        IF NOT KEYWORD_SET(yrange) THEN yrange=[min_wat,max_wat]

; Define plotting region. Default is one panel on screen
        IF N_PARAMS() NE 4 THEN BEGIN & xmaps=1 & ymaps=1 & xmap=0 & ymap=0 & ENDIF
        IF NOT KEYWORD_SET(position) THEN BEGIN
                define_panel,xmaps,ymaps,xmap,ymap,/bar,/no_charset
                position=!P.POSITION
        ENDIF
        IF KEYWORD_SET(right) THEN ystyle=5 ELSE ystyle=0

; Some settings
	IF NOT KEYWORD_SET(click) THEN BEGIN
        	IF NOT KEYWORD_SET(colour) THEN colour=foreground
		IF NOT KEYWORD_SET(zan) THEN zan=0
		IF NOT KEYWORD_SET(azm) THEN azm=0

                pnt=get_pixel(26,180) ;syowa
                
        	cx=pnt(0)    ; calculated with get_pixel
                cy=pnt(1)    ; calculated with get_pixel 

                pnt2=get_pixel2(zan,azm) ;tjornes
                
        	cx2=pnt2(0)    ; calculated with get_pixel
                cy2=pnt2(1)    ; calculated with get_pixel


             ENDIF ELSE BEGIN
		CURSOR,cx,cy,/wait
                print,cx
                print,cy
                cx=FIX(cx)
		cy=FIX(cy)
                ang=get_zan_azm(cx,cy)
               ; print,ang(0)
                ang(0)=50
               ; print,ang(0)
                zan=FIX(ang(0))
               ; print,zan
		azm=FIX(ang(1))
		clear_page
             ENDELSE
        
	PRINT,'Zenith  angle: '+STRING(zan,FORMAT='(F4.1)')
	PRINT,'Azimuth angle: '+STRING(azm,FORMAT='(F5.1)')

; Plot panel
	IF NOT KEYWORD_SET(ytitle) THEN ytitle='Auroral Intensity (Arbitrary Unit)'
        plot_time_frame,position=position-[0,0.05,0,0.05],ytitle=ytitle, $
                no_x=no_x,yrange=yrange,ystyle=ystyle,no_y=no_y,yticks=yticks

; Data for plot in TJORNES

        plot1_wat2=img_wat2(*,cx2,cy2)
        plot2_wat2=SMOOTH(plot1_wat2,3)
        plot_wat2=plot2_wat2-SMOOTH(plot2_wat2,30)+10
        
; Data for plot in SYOWA
        IF NOT KEYWORD_SET(plus) THEN plus=0
        ;time_wat,235603,240303
	plot1_wat=img_wat(*,cx,cy)
        plot2_wat=SMOOTH(plot1_wat,3)
        plot_wat=plot2_wat-SMOOTH(plot2_wat,30)+10
     
        

     
        
; Plot the data in BOTH
        IF NOT KEYWORD_SET(thick) THEN thick=!P.THICK
	x_plot=tsc_wat/3600.0
        x_plot2=tsc_wat2/3600.0
        RED='00FFF0'XL
        GREEN='00FFFF'XL
        IF stm_wat GE 240000 THEN x_plot=x_plot-24.0
        IF stm_wat2 GE 240000 THEN x_plot2=x_plot-24.0


        IF NOT KEYWORD_SET(b_ps) THEN b_ps=0     
        IF b_ps EQ 0 THEN BEGIN
          color=255
          IF n_smo EQ 1 THEN BEGIN
            OPLOT,x_plot,plot1_wat,COLOR=color,THICK=thick
          ENDIF
          IF NOT KEYWORD_SET(n_smo) THEN n_smo=0
          IF n_smo EQ 0 THEN BEGIN
            OPLOT,x_plot,plot_wat,COLOR=color,THICK=thick
          ENDIF
        ENDIF

        print,b_ps
        
        IF b_ps EQ 1 THEN BEGIN
           color=0
           IF n_smo EQ 1 THEN BEGIN
              OPLOT,x_plot,plot1_wat,COLOR=color,THICK=thick
           ENDIF
           IF NOT KEYWORD_SET(n_smo) THEN n_smo=0
           IF n_smo EQ 0 THEN BEGIN
              OPLOT,x_plot,plot_wat,COLOR=color,THICK=thick
           ENDIF
        ENDIF



        
       
        color=240
        ;FOR i=1,100-1 DO  BEGIN
        ;   OPLOT,[cx1_wat(i),cx1_wat(i)],[0,30],color=color,LINESTYLE=0,THICK=3
        ;ENDFOR

           
; Make zenith angle axis
        IF KEYWORD_SET(right) THEN BEGIN
                AXIS,yaxis=1,yrange=yrange,YSTYLE=1,TICKLEN=-!P.TICKLEN, $
                        ytitle='Auroral Intensity (Arbitrary Unit)',COLOR=colour, $
			CHARTHICK=!P.CHARTHICK
        ENDIF

; Check time info

        time_info=STRING(stm_wat,FORMAT='(I6.6)')+' to '+ $
		STRING(etm_wat,FORMAT='(I6.6)')+' UT'
	XYOUTS,!P.POSITION(2),!P.POSITION(3)+0.01,time_info,/NORMAL,ALI=1.0
	pixel_info='SYOWA Zenith:'+STRING(zan,FORMAT='(F4.1)')+' deg, '+ $
		'Azm:'+STRING(azm,FORMAT='(F5.1)')+' deg'
	XYOUTS,!P.POSITION(0),!P.POSITION(3)+0.01,pixel_info,/NORMAL,ALI=0.0

      

        print,zan
        
        print,correlate(plot_wat,plot_wat2)
    END



;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	KEO_AND_LINE
;

PRO keo_and_line,time,time2,sc1=sc1,sc2=sc2,zan=zan,azm=azm,click=click,which=which,inf=inf,_extra=e,n_smo=n_smo,b_ps=b_ps

  COMMON syowa
  COMMON tjornes  
	COMMON ps_info

; Refresh window
	clear_page

; Decide Time
  time_both_wat,time,time2

        
; Check time info
  time_info=STRING(stm_wat,FORMAT='(I6.6)')+' to '+STRING(etm_wat,FORMAT='(I6.6)')+' UT'
        

; Keep current charsize
	charsize_keep=!P.CHARSIZE

; Change charsize
	IF psinfo.open EQ 0 THEN !P.CHARSIZE=1.0*!P.CHARSIZE   ;タイトルとかの文字サイズ 
	IF psinfo.open EQ 1 THEN !P.CHARSIZE=1.0*!P.CHARSIZE

; Plot keo & line

;-------SYOWA---------

  IF which EQ 1 THEN BEGIN
    ;-----keo-------
    set_scale_wat,sc1,sc2
    plot_wat_keo_panel,1,2,0,0,line_id=0,_extra=e,zan=zan,azm=azm
  	XYOUTS,!P.POSITION(0),!P.POSITION(3)+0.01,'Keogram from MN to MS in SYOWA', /NORMAL
  	XYOUTS,!P.POSITION(2),!P.POSITION(3)+0.01,time_info,/NORMAL,ALI=1.0
    ;------line-------
    set_scale_wat,50,70
    plot_wat_line_panel,1,2,0,1,_extra=e,click=click,zan=zan,azm=azm,n_smo=n_smo,b_ps=b_ps
      
  ENDIF
        
;------TJORNES--------

  IF which EQ 0 THEN BEGIN
     
    ;-----keo--------
    set_scale_wat2,sc1,sc2
    plot_wat_keo_panel2,1,2,0,0,line_id=2,_extra=e,zan=zan,azm=azm
    XYOUTS,!P.POSITION(0),!P.POSITION(3)+0.01,'Keogram from MS to MN in TJORNES', $
    /NORMAL
    XYOUTS,!P.POSITION(2),!P.POSITION(3)+0.01,time_info,/NORMAL,ALI=1.0     
    ;-----line------
    set_scale_wat2,45,95
    plot_wat_line_panel2,1,2,0,1,_extra=e,click=click,zan=zan,azm=azm,n_smo=n_smo,b_ps=b_ps


  ENDIF

        
; Plot colour bar
	plot_wat_colour_bar,_extra=e

; Restore old charsize
	!P.CHARSIZE=charsize_keep

END

;--------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       SAVE_PS
;
PRO save_ps

  COMMON SYOWA
  COMMON TJORNES


  
  devsav=!d.name
  set_plot,'ps'
  device,filename='pic_grd/syo_keo_line_2352_2355.eps'
  device,/color
  device,xsize=19,ysize=16
  keo_and_line,235200,235500,sc1=20,sc2=90,zan=35,azm=0,which=1,n_smo=1,b_ps=1
  device,/close
  set_plot,devsav
END
        
;------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;      PLOT_FFT
;
;
       PRO plot_fft,time,time2,over=over,zan=zan,azm=azm

         COMMON tjornes
         COMMON syowa
         COMMON fft

         
;------------時間を表示------------

        time_both_wat,time,time2

         
;---------ピクセルを決める----------
        IF NOT KEYWORD_SET(zan) THEN zan=0
        IF NOT KEYWORD_SET(azm) THEN azm=0
        pnt=get_pixel(zan,azm)
        cx=pnt(0)    ; calculated with get_pixel
        cy=pnt(1)    ; calculated with get_pixel

        pnt2=get_pixel(zan,azm)
        cx2=pnt2(0)
        cy2=pnt2(1)
        
;------SMOOTH操作 of SYOWA--------

	plot1_wat=img_wat(*,cx,cy)
        plot2_wat=SMOOTH(plot1_wat,3)
        plot_wat=plot1_wat-SMOOTH(plot2_wat,30)+10
       
        
;------SMOOTH操作 of TJORNES------
        
        plot_wat_tjo1=img_wat2(*,cx2,cy2)
        plot_wat_tjo2=SMOOTH(plot_wat_tjo1,3)
        plot_wat2=plot_wat_tjo2-SMOOTH(plot_wat_tjo2,30)+10
       

;------set timerange for FFT-------

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
        
        IF sm LT 1 THEN BEGIN start_time=ss
        ENDIF ELSE BEGIN  
        start_time=60*sm+ss                ;min&sec change to sec 
        ENDELSE
        
        IF d_hour GT 0 THEN BEGIN em=em+60 ;時間をまたぐ仕組み
        ENDIF
        
        IF em eq 0 and es eq 0 THEN BEGIN end_time=3600
        ENDIF ELSE BEGIN
        end_time=60*em+es
        ENDELSE

        print,start_time
        print,end_time
        
        trim=end_time-start_time
       ; print,trim
        
        
;-------------------[plot FFT]------------------

      over=over
      fft_xrange=[0,80.0]

      define_panel,1,2,0,0,/no_charset
      save_position=!P.POSITION
      fft_position=!P.POSITION
      ytit='Amplitude'

      
;------Prepare array------

      
      num_data=trim                       ;総時間数
      t_sec=1.0                             ;時間分解能
      total_time=t_sec*(num_data-1)

      fft_results_syo=DBLARR(num_data/2)
      fft_results_tjo=DBLARR(num_data/2)
      
      freq=INDGEN(num_data/2)*t_sec      ;総時間の半分
      freq_num=INDGEN(num_data/2)

;-------Make test SIN WAVE--------

     ;test_sin=sin(findgen(1000)/10)
     ; print,SIZE(test_sin)
      
;-------TRIMING for TIME-------
      trim_data_syo=plot_wat(start_time :end_time)    
      trim_data_tjo=plot_wat2(start_time :end_time)
     ;trim_data_sin=test_sin(start_time :end_time)

;--------------FFT-------------     
      fft_results_syo=ABS(FFT(trim_data_syo))          ;ABSは絶対値をとる
      fft_results_tjo=ABS(FFT(trim_data_tjo))
     ;fft_results_sin=ABS(FFT(test_sin))

;--------エイリアシング処理--------     
      fft_results_syo=fft_results_syo[0:num_data/2-1]  
      fft_results_tjo=fft_results_tjo[0:num_data/2-1]
     ;fft_results_sin=fft_results_sin[0:num_data/2-1]
      ;print,fft_results_syo
      ;print,fft_results_tjo
     

      fft_data_syo=REVERSE(fft_results_syo)
      fft_data_tjo=REVERSE(fft_results_tjo)
     ;fft_data_sin=REVERSE(fft_results_sin)
     ;print,fft_data_syo
     ;print,fft_data_tjo


      FOR i=1,(num_data/2)-1 DO BEGIN
         freq(i)=(total_time)/freq_num(i)
      ENDFOR

;振動回数０=無限に発散する
      freq[0]=freq(1)*1.2
      freq=REVERSE(freq)

      
      fft_yrange=[0,1.0]

      fft_data_tjo=SMOOTH(fft_data_tjo,2)
      fft_data_syo=SMOOTH(fft_data_syo,2)

;------Plot FFT of TJORNES
      
         PLOT,[0],[0],title=title,xrange=fft_xrange,yrange=fft_yrange, $
          xtitle=xtit,ytitle=ytit,position=fft_position+[0,0.04,0,0], $
          XTICKFORMAT=xtickformat,YTICKNAME=ytickname,CHARSIZE=!P.CHARSIZE*1.5

         OPLOT,freq,fft_data_tjo,color=color,LINESTYLE=0,SYMSIZE=10,THICK=2
         color=240
         OPLOT,[24,24],[0,2],color=color,LINESTYLE=0,THICK=3
         OPLOT,[39.8,39.8],[0,2],color=color,LINESTYLE=0,THICK=3
         OPLOT,[20,20],[0,2],color=color,LINESTYLE=0,THICK=3
      
;------Check time info------

        time_info=STRING(stm_wat,FORMAT='(I6.6)')+' to '+ $
		STRING(etm_wat,FORMAT='(I6.6)')+' UT'
	XYOUTS,!P.POSITION(2),!P.POSITION(3)+0.01,time_info,/NORMAL,ALI=1.0
	pixel_info='FFT data: '+STRING(zan,FORMAT='(F4.1)')+' deg,  '+ $
		'Azm: '+STRING(azm,FORMAT='(F5.1)')+' deg in TJORNES'
	XYOUTS,!P.POSITION(0),!P.POSITION(3)+0.01,pixel_info,/NORMAL,ALI=0.0


  END
;------------------------------------------------------------------------------------------------------------------------
; NAME:
;
;      PLOT_FFT2
;
;
       PRO plot_fft2,time,time2,over=over,zan=zan,azm=azm

         COMMON tjornes
         COMMON syowa
         COMMON fft

         
;------------時間を表示------------

        time_both_wat,time,time2

         
;---------ピクセルを決める----------
        IF NOT KEYWORD_SET(zan) THEN zan=0
        IF NOT KEYWORD_SET(azm) THEN azm=0
        pnt=get_pixel(zan,azm)
        cx=pnt(0)    ; calculated with get_pixel
        cy=pnt(1)    ; calculated with get_pixel

        pnt2=get_pixel(zan,azm)
        cx2=pnt2(0)
        cy2=pnt2(1)
        
;------SMOOTH操作 of SYOWA--------

	plot1_wat=img_wat(*,cx,cy)
        plot2_wat=SMOOTH(plot1_wat,3)
        plot_wat=plot1_wat-SMOOTH(plot2_wat,30)+10
       
        
;------SMOOTH操作 of TJORNES------
        
        plot_wat_tjo1=img_wat2(*,cx2,cy2)
        plot_wat_tjo2=SMOOTH(plot_wat_tjo1,3)
        plot_wat2=plot_wat_tjo2-SMOOTH(plot_wat_tjo2,30)+10
       

;------set timerange for FFT-------

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
        
        IF sm LT 1 THEN BEGIN start_time=ss
        ENDIF ELSE BEGIN  
        start_time=60*sm+ss                ;min&sec change to sec 
        ENDELSE
        
        IF d_hour GT 0 THEN BEGIN em=em+60 ;時間をまたぐ仕組み
        ENDIF
        
        IF em eq 0 and es eq 0 THEN BEGIN end_time=3600
        ENDIF ELSE BEGIN
        end_time=60*em+es
        ENDELSE

        print,start_time
        print,end_time
        
        trim=end_time-start_time
       ; print,trim
        
        
;-------------------[plot FFT]------------------

      over=over
      fft_xrange=[0,80.0]

      define_panel,1,2,0,1,/no_charset
      fft_position=!P.POSITION
      ytit='Amplitude'

      
;------Prepare array------

      
      num_data=trim                       ;総時間数
      t_sec=1.0                             ;時間分解能
      total_time=t_sec*(num_data-1)

      fft_results_syo=DBLARR(num_data/2)
      fft_results_tjo=DBLARR(num_data/2)
      
      freq=INDGEN(num_data/2)*t_sec      ;総時間の半分
      freq_num=INDGEN(num_data/2)

;-------Make test SIN WAVE--------

     ;test_sin=sin(findgen(1000)/10)
     ; print,SIZE(test_sin)
      
;-------TRIMING for TIME-------
      trim_data_syo=plot_wat(start_time :end_time)    
      trim_data_tjo=plot_wat2(start_time :end_time)
     ;trim_data_sin=test_sin(start_time :end_time)

;--------------FFT-------------     
      fft_results_syo=ABS(FFT(trim_data_syo))          ;ABSは絶対値をとる
      fft_results_tjo=ABS(FFT(trim_data_tjo))
     ;fft_results_sin=ABS(FFT(test_sin))

;--------エイリアシング処理--------     
      fft_results_syo=fft_results_syo[0:num_data/2-1]  
      fft_results_tjo=fft_results_tjo[0:num_data/2-1]
     ;fft_results_sin=fft_results_sin[0:num_data/2-1]
      ;print,fft_results_syo
      ;print,fft_results_tjo
     

      fft_data_syo=REVERSE(fft_results_syo)
      fft_data_tjo=REVERSE(fft_results_tjo)
     ;fft_data_sin=REVERSE(fft_results_sin)
     ;print,fft_data_syo
     ;print,fft_data_tjo


      FOR i=1,(num_data/2)-1 DO BEGIN
         freq(i)=(total_time)/freq_num(i)
      ENDFOR

;振動回数０=無限に発散する
      freq[0]=freq(1)*1.2
      freq=REVERSE(freq)

      
      fft_yrange=[0,1.0]

      fft_data_tjo=SMOOTH(fft_data_tjo,2)
      fft_data_syo=SMOOTH(fft_data_syo,2)

    

;------plot syowa fft-------
      RED='0F00F0'XL
      PLOT,[0],[0],title=title,xrange=fft_xrange,yrange=fft_yrange, $
           xtitle=xtit,ytitle=ytit,position=fft_positio, $
           XTICKFORMAT=xtickformat,YTICKNAME=ytickname,CHARSIZE=!P.CHARSIZE*1.5
       
      OPLOT,freq,fft_data_syo,color=color,LINESTYLE=0,THICK=2
      color=240
      OPLOT,[24,24],[0,1],color=color,LINESTYLE=0,THICK=3
      OPLOT,[30,30],[0,1],color=color,LINESTYLE=0,THICK=3
      OPLOT,[40,40],[0,1],color=color,LINESTYLE=0,THICK=3
      
;------Check time info------

          time_info=STRING(stm_wat,FORMAT='(I6.6)')+' to '+ $
		STRING(etm_wat,FORMAT='(I6.6)')+' UT'
	XYOUTS,!P.POSITION(2),!P.POSITION(3)+0.01,time_info,/NORMAL,ALI=1.0
	pixel_info='FFT data: '+STRING(zan,FORMAT='(F4.1)')+' deg,  '+ $
		'Azm: '+STRING(azm,FORMAT='(F5.1)')+' deg in SYOWA'
	XYOUTS,!P.POSITION(0),!P.POSITION(3)+0.01,pixel_info,/NORMAL,ALI=0.0

  END

;--------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       FFT_BOTH
;
        PRO fft_both,time,time2,zan=zan,azm=azm

          COMMON SYOWA
          COMMON TJORNES

          !p.multi=[0,0,2,0,0]
          plot_fft2,time,time2,zan=zan,azm=azm
          plot_fft,time,time2,zan=zan,azm=azm
          


       END
       
;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	ZAN_CORRELATE
;
        PRO zan_correlate,time,time2,zan=zan,azm=azm

          COMMON TJORNES
          COMMON SYOWA




          memory_cor=DBLARR(151)
          zan_angle=FINDGEN(151)

          FOR i=0,151-1 DO BEGIN
             zan_angle(i)=zan_angle(i)-75
          ENDFOR

          

          
          stime=time
          etime=time2

          
          sh=FIX(stime / 10000)       ;start hour
          sm=FIX((stime mod 10000)/100) ;start min
          ss=stime mod 100              ;start sec
          
          eh=FIX(etime / 10000)       ;end hour
          em=FIX((etime mod 10000)/100) ;end min
          es=etime mod 100              ;end sec

          d_hour=eh-sh

          IF sh EQ 24 AND eh EQ 24 THEN BEGIN
             sm=sm+60
             em=em+60
          ENDIF
          
          IF sm LT 1 THEN BEGIN start_time=ss
          ENDIF ELSE BEGIN  
             start_time=60*sm+ss ;min&sec change to sec 
          ENDELSE
          
          IF d_hour GT 0 THEN BEGIN em=em+60 ;時間をまたぐ仕組み
          ENDIF
          
          IF em eq 0 and es eq 0 THEN BEGIN end_time=3600
          ENDIF ELSE BEGIN
             end_time=60*em+es
          ENDELSE

          print,start_time
          print,end_time
          

          
          IF NOT KEYWORD_SET(zan) THEN zan=0
          IF NOT KEYWORD_SET(azm) THEN azm=0 
;---TJORNES------
          pnt2=get_pixel2(zan,azm)
          cx2=pnt2(0)
          cy2=pnt2(1)
          plot_wat_tjo1=img_wat2(*,cx2,cy2)
          plot_wat_tjo2=SMOOTH(plot_wat_tjo1,3)
          plot_wat2=plot_wat_tjo2-SMOOTH(plot_wat_tjo2,30)+10
          plot_wat2=plot_wat2(start_time :end_time)
;---SYOWA------
          pnt=get_pixel(0,0)
          cx=pnt(0)
          cy=pnt(1)
          plot1_wat=img_wat(*,cx,cy)
          plot2_wat=SMOOTH(plot1_wat,3)
          plot_wat=plot2_wat-SMOOTH(plot2_wat,30)+10
          plot_wat=plot_wat(start_time :end_time)
         
          cor_pos_old=correlate(plot_wat2,plot_wat)
          memory_cor(75)=cor_pos_old

          
;----correlate------       
          FOR i=1,75-1 DO BEGIN
             pnt=get_pixel(i,0)
             cx=pnt(0)
             cy=pnt(1)

             plot1_wat=img_wat(*,cx,cy)
             plot2_wat=SMOOTH(plot1_wat,3)
             plot_prov=plot2_wat-SMOOTH(plot2_wat,30)+10
             plot_prov=plot_prov(start_time :end_time)
             
             cor_pos_new=correlate(plot_wat2,plot_prov)
             memory_cor(i+75)=cor_pos_new
             
             IF cor_pos_old GE cor_pos_new THEN BEGIN
                plot_wat=plot_wat
                cor_pos_old=cor_pos_old
             ENDIF ELSE BEGIN
                plot_wat=plot_prov
                cor_pos_old=cor_pos_new
                answer_zan=i
             ENDELSE
             
          ENDFOR

         

          
;----correlate azm=180------
          FOR i=1,75-1 DO BEGIN
             pnt=get_pixel(i,180)
             cx=pnt(0)
             cy=pnt(1)

             plot1_wat=img_wat(*,cx,cy)
             plot2_wat=SMOOTH(plot1_wat,3)
             plot_prov=plot2_wat-SMOOTH(plot2_wat,30)+10
             plot_prov=plot_prov(start_time :end_time)
             
             cor_pos_new=correlate(plot_wat2,plot_prov)
             memory_cor(75-i)=cor_pos_new
             
             IF cor_pos_new LT cor_pos_old THEN BEGIN
                plot_wat=plot_wat
                cor_pos_old=cor_pos_old
             ENDIF ELSE BEGIN
                plot_wat=plot_prov
                cor_pos_old=cor_pos_new
                answer_zan=-i
             ENDELSE
          ENDFOR
          

         

          
          
;------Plot SOUKAN_KUN--------



          define_panel,1,2,0,0,/no_charset
          position=!P.POSITION
          xtit='zenith angle'
          ytit='Correlation'
          title='Shifting Zenith Angle'
          xrange=[-75,75]
          yrange=[-0.2,0.4]
          PLOT,[0],[0],title=title,xrange=xrange,yrange=yrange, $
          xtitle=xtit,ytitle=ytit,position=position, $
          XTICKFORMAT=xtickformat,YTICKNAME=ytickname,CHARSIZE=!P.CHARSIZE*1.5
        
          OPLOT,zan_angle,memory_cor,color=color,LINESTYLE=0,THICK=2


          
           print,answer_zan
           print,cor_pos_old


        END


;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	TIME_CORRELATE
;
        PRO time_correlate,time,time2,zan=zan,azm=azm

          COMMON TJORNES
          COMMON SYOWA




        memory_cor=DBLARR(61)
        shift_time=FINDGEN(61)

        
        FOR i=0,61-1 DO BEGIN
           shift_time(i)=shift_time(i)-30
        ENDFOR


        
        stime=time
        etime=time2

          
        sh=FIX(stime / 10000)              ;start hour
        sm=FIX((stime mod 10000)/100)    ;start min
        ss=stime mod 100                 ;start sec
        
        eh=FIX(etime / 10000)            ;end hour
        em=FIX((etime mod 10000)/100)    ;end min
        es=etime mod 100                 ;end sec

        d_hour=eh-sh

        IF sh EQ 24 AND eh EQ 24 THEN BEGIN
        sm=sm+60
        em=em+60
        ENDIF
        
        IF sm LT 1 THEN BEGIN start_time=ss
        ENDIF ELSE BEGIN  
        start_time=60*sm+ss                ;min&sec change to sec 
        ENDELSE
        
        IF d_hour GT 0 THEN BEGIN em=em+60 ;時間をまたぐ仕組み
        ENDIF
   
        IF em eq 0 and es eq 0 THEN BEGIN end_time=3600
        ENDIF ELSE BEGIN
        end_time=60*em+es
        ENDELSE

        print,start_time
        print,end_time
        

          
          IF NOT KEYWORD_SET(zan) THEN zan=0
          IF NOT KEYWORD_SET(azm) THEN azm=0 

;---TJORNES------
          pnt2=get_pixel2(zan,azm)
          cx2=pnt2(0)
          cy2=pnt2(1)

          plot_wat_tjo1=img_wat2(*,cx2,cy2)
          plot_wat_tjo2=SMOOTH(plot_wat_tjo1,3)
          plot_wat2=plot_wat_tjo2-SMOOTH(plot_wat_tjo2,30)+10
          plot_wat2=plot_wat2(start_time :end_time)

;---SYOWA------
          pnt=get_pixel(zan,azm)
          cx=pnt(0)
          cy=pnt(1)

          plot1_wat=img_wat(*,cx,cy)
          plot2_wat=SMOOTH(plot1_wat,3)
          plot_wat=plot2_wat-SMOOTH(plot2_wat,30)+10
          plot_wat=plot_wat(start_time :end_time)

;----First correlation----
          cor_pos_old=correlate(plot_wat2,plot_wat)
          memory_cor(30)=cor_pos_old
          print,cor_pos_old
          
;----correlate------       
          FOR i=1,31-1 DO BEGIN
             pnt=get_pixel(zan,azm)
             cx=pnt(0)
             cy=pnt(1)

             plot1_wat=img_wat(*,cx,cy)
             plot2_wat=SMOOTH(plot1_wat,3)
             plot_prov=plot2_wat-SMOOTH(plot2_wat,30)+10
             plot_prov=plot_prov(start_time+i :end_time+i)
             
             cor_pos_new=correlate(plot_wat2,plot_prov)
             memory_cor(30+i)=cor_pos_new
             print,cor_pos_new

             
             IF cor_pos_old GE cor_pos_new THEN BEGIN
                plot_wat=plot_wat
                cor_pos_old=cor_pos_old
             ENDIF ELSE BEGIN
                plot_wat=plot_prov            ;memorize 
                cor_pos_old=cor_pos_new
                answer_time=i
             ENDELSE
             
          ENDFOR


          
;----correlate azm=180------
          FOR i=1,31-1 DO BEGIN
             pnt=get_pixel(zan,azm)
             cx=pnt(0)
             cy=pnt(1)

             plot1_wat=img_wat(*,cx,cy)
             plot2_wat=SMOOTH(plot1_wat,3)
             plot_prov=plot2_wat-SMOOTH(plot2_wat,3)+10
             plot_prov=plot_prov(start_time-i :end_time-i)
             
             cor_pos_new=correlate(plot_wat2,plot_prov)
             memory_cor(30-i)=cor_pos_new
             
             IF cor_pos_new LT cor_pos_old THEN BEGIN
                plot_wat=plot_wat
                cor_pos_old=cor_pos_old
             ENDIF ELSE BEGIN
                plot_wat=plot_prov
                cor_pos_old=cor_pos_new
                answer_time=-i
             ENDELSE
          ENDFOR
          

         

          
          
;------Plot SOUKAN_KUN--------



          define_panel,1,2,0,0,/no_charset
          position=!P.POSITION
          xtit='Shift Time'
          ytit='Correlation'
          title='Shifting time & Correlation'
          xrange=[-30,30]
          yrange=[-0.2,0.4]
          PLOT,[0],[0],title=title,xrange=xrange,yrange=yrange, $
          xtitle=xtit,ytitle=ytit,position=position, $
          XTICKFORMAT=xtickformat,YTICKNAME=ytickname,CHARSIZE=!P.CHARSIZE*1.5
        
          OPLOT,shift_time,memory_cor,color=color,LINESTYLE=0,THICK=2





          
          memory_cor_max=max(memory_cor)
          print,memory_cor_max
          print,answer_time
       END
        

;------------------------------------------------------------------------------------------------------------------

; NAME:
;
;	SOUKAN_SCATTER
;
        PRO soukan_scatter,time,time2,zan=zan,azm=azm,plus=plus

        COMMON TJORNES
        COMMON SYOWA




          stime=time
          etime=time2

          
          sh=FIX(stime / 10000)         ;start hour
          sm=FIX((stime mod 10000)/100) ;start min
          ss=stime mod 100              ;start sec
          
          eh=FIX(etime / 10000)         ;end hour
          em=FIX((etime mod 10000)/100) ;end min
          es=etime mod 100              ;end sec

          d_hour=eh-sh

          IF sh EQ 24 AND eh EQ 24 THEN BEGIN
             sm=sm+60
             em=em+60
          ENDIF
          
          IF sm LT 1 THEN BEGIN start_time=ss
          ENDIF ELSE BEGIN  
             start_time=60*sm+ss        ;min&sec change to sec 
          ENDELSE
          
          IF d_hour GT 0 THEN BEGIN em=em+60 ;時間をまたぐ仕組み
          ENDIF
          
          IF em eq 0 and es eq 0 THEN BEGIN end_time=3600
          ENDIF ELSE BEGIN
             end_time=60*em+es
          ENDELSE

          print,start_time
          print,end_time
          

          
          IF NOT KEYWORD_SET(zan) THEN zan=0
          IF NOT KEYWORD_SET(azm) THEN azm=0 
          IF NOT KEYWORD_SET(plus) THEN plus=0 

          
;---TJORNES------
          pnt2=get_pixel2(zan,azm)
          cx2=pnt2(0)
          cy2=pnt2(1)

          plot_wat_tjo1=img_wat2(*,cx2,cy2)
          plot_wat_tjo2=SMOOTH(plot_wat_tjo1,3)
          plot_wat2=plot_wat_tjo2-SMOOTH(plot_wat_tjo2,30)+10
          plot_wat2=plot_wat2(start_time :end_time)
          
;---SYOWA------
          pnt=get_pixel(zan,azm)
          cx=pnt(0)
          cy=pnt(1)

          plot1_wat=img_wat(*,cx,cy)
          plot2_wat=SMOOTH(plot1_wat,3)
          plot_wat=plot2_wat-SMOOTH(plot2_wat,30)+10
          plot_wat=plot_wat(start_time+plus :end_time+plus)

;-----Correlation------

          cor=correlate(plot_wat2,plot_wat)
          ;print,cor
          measure_errors=SQRT(ABS(plot_wat))
          line_fit=linfit(plot_Wat2,plot_wat,MEASURE_errors=measure_errors)
          ;print,line_fit,measure_errors
;-------Plot Soukan Scatter------

          define_panel,1,2,0,0,/no_charset
          position=!P.POSITION
          xtit='TJORNES'
          ytit='SYOWA'
          ;title='Intensity scatter'
          xrange=[0,25]
          yrange=[0,30]
          PLOT,[0],[0],title=title,xrange=xrange,yrange=yrange, $
          xtitle=xtit,ytitle=ytit,position=position, $
          XTICKFORMAT=xtickformat,YTICKNAME=ytickname,CHARSIZE=!P.CHARSIZE*1.5
        
          OPLOT,plot_wat2,plot_wat,color=color,LINESTYLE=0,THICK=2,psym=1
          color=240
          ;OPLOT,MEASURE_errors,color=color,LINESTYLE=0,THICK=2
          
;------time & pixel & correlation------

        time_info=STRING(stime,FORMAT='(I6.6)')+' to '+ $
		STRING(etime,FORMAT='(I6.6)')+' UT'
	XYOUTS,!P.POSITION(2),!P.POSITION(3)+0.01,time_info,/NORMAL,ALI=1.0
	pixel_info='Intensity scatter : '+STRING(zan,FORMAT='(F4.1)')+' deg,  '+ $
		'Azm: '+STRING(azm,FORMAT='(F5.1)')+' deg'
	XYOUTS,!P.POSITION(0),!P.POSITION(3)+0.01,pixel_info,/NORMAL,ALI=0.0
        cor_info='Correlation : '+STRING(cor,FORMAT='(f0.9)')
     	XYOUTS,!P.POSITION(0),!P.POSITION(3)+0.05,cor_info,/NORMAL,ALI=0.0

     END
        


        
;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	PLOT_WAT_COLOUR_BAR
;

	PRO plot_wat_colour_bar,ymaps,ymap,position=position,ctab=ctab,_extra=e,leg_pos=leg_pos

	COMMON syowa

	COMMON prefs

; Update colour info
        IF N_ELEMENTS(ctab) EQ 0 THEN BEGIN
                set_colour_table,/default
        ENDIF ELSE BEGIN
                set_colour_table,ctab
        ENDELSE
        update_colour_info

; Setup scale
	old_minval=minval
	old_maxval=maxval
	set_scale,min_wat,max_wat,/quiet

; Reverse issue
	old_param=0
	IF parameter EQ 'vel' THEN BEGIN
		old_param=1
		pwr_l
		set_scale,min_wat,max_wat,/quiet
	ENDIF

; Plot colour bar
	IF KEYWORD_SET(position) THEN BEGIN
		plot_colour_bar,position=position,_extra=e, $
			legend='Optical Intensity!C!DArbitrary Unit!N',leg_pos=leg_pos,/no_gnd
			;legend='Optical Intensity (Arbitrary Unit)',leg_pos=leg_pos,/no_gnd
	ENDIF ELSE BEGIN
		IF N_PARAMS() NE 2 THEN BEGIN
			ymaps=1 & ymap=0
		ENDIF
		plot_colour_bar,ymaps,ymap,_extra=e, $
			legend='Optical Intensity!C!DArbitrary Unit!N',leg_pos=leg_pos,/no_gnd
			;legend='Optical Intensity (Arbitrary Unit)',leg_pos=leg_pos,/no_gnd
	ENDELSE

; Restore parameter
	IF old_param EQ 1 THEN vel

; Restore scale
	set_scale,old_minval,old_maxval,/quiet

	END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       OVERLAY_MYOPIC_DATA
;

        PRO overlay_myopic_data,_extra=e

        overlay_myopic_fov,_extra=e,/myopic_data

        END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       OVERLAY_MYOPIC_FOV
;

        PRO overlay_myopic_fov,beams_plot=beams_plot,gates_plot=gates_plot,colour=colour,thick=thick, $
		beams_sep=beams_sep,gates_sep=gates_sep,irev=irev,myopic_data=myopic_data, $
		no_time_set=no_time_set,altitude=altitude

        COMMON syowa

        COMMON data
        COMMON prefs
        COMMON colour_info

; Variables
        radar_lat=63.77
        radar_lon=-20.54
        radar_boresite=30.0
        gate_length=15.0
        first_range=15.0
	IF NOT KEYWORD_SET(altitude) THEN altitude=110.0

; Half width of ATV image
        wid=siz_wat/2

; Find radar cell locations
        x_wat=FLTARR(17,76)
        y_wat=FLTARR(17,76)
        FOR beam=0,16 DO BEGIN
                FOR gate=0,75 DO BEGIN
                        fov_pos=find_cell(radar_lat,radar_lon,16,radar_boresite, $
                                3.24,gate_length,first_range,beam-0.5,gate)
                        lat=fov_pos(0)
                        lon=fov_pos(1)
                        zen=get_zenith_angle_from_tjornes(lat,lon,altitude)
                        IF NOT KEYWORD_SET(irev) THEN BEGIN
                                ;azm=get_azimuth_from_tjornes(lat,lon)+90.0+rog_wat
                                azm=get_azimuth_from_tjornes(lat,lon)
                        ENDIF ELSE BEGIN
                                ;azm=-get_azimuth_from_tjornes(lat,lon)+90.0;+rog_wat
                                azm=-get_azimuth_from_tjornes(lat,lon)
                        ENDELSE
                        IF azm GE 360.0 THEN azm=azm-360.0
                        tmp_pix=get_pixel(zen,azm)
                        x_wat(beam,gate)=tmp_pix(0)
                        y_wat(beam,gate)=tmp_pix(1)
                        ;x_wat(beam,gate)=wid+zen*COS(!PI*azm/180.0)*wid/75.0
                        ;y_wat(beam,gate)=wid+zen*SIN(!PI*azm/180.0)*wid/75.0
                ENDFOR
        ENDFOR

        IF NOT KEYWORD_SET(colour) THEN colour=254
        IF NOT KEYWORD_SET(thick) THEN thick=1

; If you want to plot the data as well ...
        IF KEYWORD_SET(myopic_data) THEN BEGIN

; Get data arrays
		IF NOT KEYWORD_SET(no_time_set) THEN BEGIN
                	go_time,FIX(STRMID(STRING(hms_wat(now_wat),FORMAT='(I6.6)'),0,4)), $
                        	FIX(STRMID(STRING(hms_wat(now_wat),FORMAT='(I6.6)'),4,2))
		ENDIF
                varr=get_scan()
                grnd=get_scan(/gnd_flag)

; Set color bar and levels, and switch colour map for velocity plot
                cin=FIX(FINDGEN(no_colours)*(ncol-4)/(no_colours-1))+1
                lvl=minval+FINDGEN(no_colours)*(maxval-minval)/no_colours
                IF parameter EQ 'vel' THEN cin=ROTATE(cin,2)

; Plot data
                FOR m=0,radar_beams-1 DO BEGIN
                        FOR n=11,75-1 DO BEGIN
                                IF varr(m,n) NE 10000 THEN BEGIN
                                        IF NOT((grnd(m,n) EQ 0 AND scatter EQ 1) OR (grnd(m,n) NE 0 AND scatter EQ 2)) THEN BEGIN
                                                ind=(MAX(WHERE(lvl LE varr(m,n))) > 0)
                                                IF parameter EQ 'vel' AND scatter EQ 3 AND grnd(m,n) EQ 1 THEN col=grey ELSE col=cin(ind)
                                                x_poly=[x_wat(m,n),x_wat(m+1,n),x_wat(m+1,n+1), $
                                                        x_wat(m,n+1),x_wat(m,n)]
                                                y_poly=[y_wat(m,n),y_wat(m+1,n),y_wat(m+1,n+1), $
                                                        y_wat(m,n+1),y_wat(m,n)]
                                                set_colour_table,/default
                                                POLYFILL,x_poly,y_poly,COL=col,NOCLIP=0
                                        ENDIF
                                ENDIF
                        ENDFOR
                ENDFOR
                set_scale_wat,min_wat,max_wat,/quiet
        ENDIF

; Plot the fov on ATV image
        IF KEYWORD_SET(beams_plot) AND KEYWORD_SET(gates_plot) THEN BEGIN
                FOR i=beams_plot,beams_plot DO BEGIN
                        FOR j=gates_plot(0),gates_plot(1) DO BEGIN
                                x_poly=[x_wat(i,j),x_wat(i+1,j),x_wat(i+1,j+1), $
                                        x_wat(i,j+1),x_wat(i,j)]
                                y_poly=[y_wat(i,j),y_wat(i+1,j),y_wat(i+1,j+1), $
                                        y_wat(i,j+1),y_wat(i,j)]
                                ;OPLOT,x_poly,y_poly,COL=254,THICK=thick*2
                                OPLOT,x_poly,y_poly,COL=colour,THICK=thick
                        ENDFOR
                ENDFOR
        ENDIF ELSE BEGIN

; Plot rough grids
                IF NOT KEYWORD_SET(normal) THEN BEGIN
                        FOR i=10,75,10 DO OPLOT,[x_wat(5:11,i)],[y_wat(5:11,i)],COL=colour,THICK=thick
                        ;FOR i=10,75,10 DO OPLOT,[x_wat(*,i)],[y_wat(*,i)],COL=colour,THICK=thick
                        FOR i=5,11,1 DO OPLOT,[x_wat(i,10:75)],[y_wat(i,10:75)],COL=colour,THICK=thick
                        ;FOR i=0,16, 8 DO OPLOT,[x_wat(i,10:75)],[y_wat(i,10:75)],COL=colour,THICK=thick
                ENDIF ELSE BEGIN
                        FOR i=0,75,15 DO PLOTS,[x_wat(*,i)],[y_wat(*,i)],COL=colour, $
                                THICK=thick,NORMAL=normal
                        FOR i=0,16, 8 DO PLOTS,[x_wat(i,*)],[y_wat(i,*)],COL=colour, $
                                THICK=thick,NORMAL=normal
                ENDELSE
        ENDELSE

        END
;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       OVERLAY_MYOPIC_ON_KEOGRAM
;

        PRO overlay_myopic_on_keogram

        COMMON syowa

        COMMON data
        COMMON prefs
        COMMON colour_info

; Variables
        radar_lat=63.77
        radar_lon=-20.54
        radar_boresite=30.0
        gate_length=15.0
        first_range=15.0
        altitude=110.0
        lat_tjo=66.20

; Find radar cell locations
        zan_cell=FLTARR(16,76)
        FOR beam=0,15 DO BEGIN
                FOR gate=11,75 DO BEGIN
                        fov_pos=find_cell(radar_lat,radar_lon,16,radar_boresite, $
                                3.24,gate_length,first_range,beam,gate)
                        lat=fov_pos(0)
                        lon=fov_pos(1)
                        zan_cell(beam,gate)=get_zenith_angle_from_tjornes(lat,lon,altitude)
                        IF lat LT lat_tjo THEN zan_cell(beam,gate)=-zan_cell(beam,gate)
                ENDFOR
        ENDFOR

; Find beam data blocks which have the correct look direction and frequency
        beamno=beam_no
        plot_beams=WHERE(beam_dir EQ beamno AND beam_time GE rti_start-120 AND beam_time LE rti_end+120, $
                no_plot_beams)

; Set color bar and levels, and switch colour map for velocity plot
        cin=FIX(FINDGEN(no_colours)*(ncol-4)/(no_colours-1))+1
        lvl=minval+FINDGEN(no_colours)*(maxval-minval)/no_colours
        IF parameter EQ 'vel' THEN cin=ROTATE(cin,2)

; Cycle through beams to plot
        day_start=FIX(rti_start/86400L)*86400L
        FOR m=0,no_plot_beams-2 DO BEGIN
                FOR n=11,74 DO BEGIN

                        IF a(n,plot_beams(m)) NE 10000 THEN BEGIN

                                start_time=beam_time(plot_beams(m))-day_start
                                end_time=beam_time(plot_beams(m+1))-day_start

                                IF NOT((gscat(n,plot_beams(m)) EQ 0 AND scatter EQ 1) OR $
                                        (gscat(n,plot_beams(m)) NE 0 AND scatter EQ 2)) THEN BEGIN
                                
                                                ind=(MAX(WHERE(lvl LE ((a(n,plot_beams(m)) > minval) < maxval))) > 0)
                                                IF parameter EQ 'vel' AND scatter EQ 3 AND gscat(n,plot_beams(m)) EQ 1 THEN col=grey ELSE col=cin(ind)
                                                POLYFILL,[start_time/3600.0,start_time/3600.0,end_time/3600.0,end_time/3600.0], $
                                                         [zan_cell(beamno,n),zan_cell(beamno,n+1),zan_cell(beamno,n+1),zan_cell(beamno,n)],COL=col,NOCLIP=0

                                ENDIF
                        ENDIF
                ENDFOR
        ENDFOR

        END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       OVERLAY_ZENITH
;

        PRO overlay_zenith

        COMMON syowa

        COMMON colour_info
        update_colour_info

; Zenith angle limit
	zenith_limit=75.0

; Plot circles (from 10 to 70, div is 10 deg)
        xx=FLTARR(361) & yy=FLTARR(361)
        c_img=[siz_wat/2,siz_wat/2+0.5]
        FOR j=15,75,15 DO BEGIN
                d_img=FLOAT(j)*(siz_wat/2)/zenith_limit
                FOR i=0,360 DO BEGIN
                        xx(i)=c_img(0)+1*d_img*COS(i*!PI/180.0)
                        yy(i)=c_img(1)+1*d_img*SIN(i*!PI/180.0)
                ENDFOR
                set_colour_table,/default
                OPLOT,xx,yy,COL=white,NOCLIP=0,LINE=1
        ENDFOR

	END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       OVERLAY_ZENITH_SPECIAL
;

        PRO overlay_zenith_special

        COMMON syowa

        COMMON colour_info
        update_colour_info

; Zenith angle limit
	zenith_limit=75.0

; Plot circles (from 10 to 70, div is 10 deg)
        xx=FLTARR(361) & yy=FLTARR(361)
        c_img=[siz_wat/2,siz_wat/2+0.5]
	ele_ang=[30,60,75]
        FOR k=0,2 DO BEGIN
		j=ele_ang(k)
                d_img=FLOAT(j)*(siz_wat/2)/zenith_limit
                FOR i=0,360 DO BEGIN
                        xx(i)=c_img(0)+1*d_img*COS(i*!PI/180.0)
                        yy(i)=c_img(1)+1*d_img*SIN(i*!PI/180.0)
                ENDFOR
                set_colour_table,/default
                OPLOT,xx,yy,COL=254,NOCLIP=0,LINE=0,THICK=6
                OPLOT,xx,yy,COL=1,NOCLIP=0,LINE=0,THICK=4
        ENDFOR

	END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       MASK_FOV
;

        PRO mask_fov

        COMMON syowa

        xx=FLTARR(361) & yy=FLTARR(361) & xe=FLTARR(361) & ye=FLTARR(361)
        c_img=[siz_wat/2,siz_wat/2]
        d_img=FIX(a_v_wat*!PI/2.0)
        FOR i=0,360 DO BEGIN
                xx(i)=c_img(0)+1*d_img*COS(i*!PI/180.0)
                yy(i)=c_img(1)+1*d_img*SIN(i*!PI/180.0)
                xe(i)=c_img(0)+2*d_img*COS(i*!PI/180.0)
                ye(i)=c_img(1)+2*d_img*SIN(i*!PI/180.0)
        ENDFOR
        set_colour_table,/default
        FOR i=0,359 DO BEGIN
                POLYFILL,[xx(i),xe(i),xe(i+1),xx(i+1)],[yy(i),ye(i),ye(i+1),yy(i+1)], $
                        COL=253,NOCLIP=0
        ENDFOR

	END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	DEFINE_ZENITH_TJORNES
;

	PRO define_zenith_tjornes,preview=preview

; Image used for defining zenith
	fname='./wat_tjo_11-09-30_22-00-04-35.jpg'

; Time info for input to zenith2000
	date='20110930'
	hh='22'
	mmss='0004'

; Read image and identify size
	READ_JPEG,fname,image
	image_size=SIZE(image,/DIM)
	xsize=image_size(0)
	ysize=image_size(1)
	PRINT,xsize,ysize
	WINDOW,XSIZE=xsize,YSIZE=ysize

; Plot image
	LOADCT,0
	LOADCT,39
	TVSCL,image
	!P.POSITION=[0,0,1,1]
	PLOT,[0],[0],XRANGE=[0,xsize],YRANGE=[0,ysize],/XSTYLE,/YSTYLE

	IF KEYWORD_SET(preview) THEN RETURN

	psym8

; Get position of Capera (StarID=4)
	PRINT,'Click Capera (in Auriga)'
	CURSOR,x_tmp,y_tmp,/DEVICE
	SPAWN,'sleep 1'
	P=MAX(image(x_tmp-5:x_tmp+5,y_tmp-5:y_tmp+5),i)
	ind=ARRAY_INDICES(image(x_tmp-5:x_tmp+5,y_tmp-5:y_tmp+5),i)
	x1=x_tmp-5+ind(0)
	y1=y_tmp-5+ind(1)
	x1=x_tmp
	y1=y_tmp
	OPLOT,[x1],[y1],PSYM=8,SYMSIZE=1,COL=254
	PRINT,'Capera  =[',x1,y1,']'

; Get position of Vega (StarID=11)
	PRINT,'Click Vega (in Lyra)'
	CURSOR,x_tmp,y_tmp,/DEVICE
	SPAWN,'sleep 1'
	P=MAX(image(x_tmp-5:x_tmp+5,y_tmp-5:y_tmp+5),i)
	ind=ARRAY_INDICES(image(x_tmp-5:x_tmp+5,y_tmp-5:y_tmp+5),i)
	x2=x_tmp-5+ind(0)
	y2=y_tmp-5+ind(1)
	x2=x_tmp
	y2=y_tmp
	OPLOT,[x2],[y2],PSYM=8,SYMSIZE=1,COL=254
	PRINT,'Vega    =[',x2,y2,']'

; Get position of Polaris (StarID=1)
	PRINT,'Click Polaris (in Ursa Minor)'
	CURSOR,x_tmp,y_tmp,/DEVICE
	SPAWN,'sleep 1'
	P=MAX(image(x_tmp-5:x_tmp+5,y_tmp-5:y_tmp+5),i)
	ind=ARRAY_INDICES(image(x_tmp-5:x_tmp+5,y_tmp-5:y_tmp+5),i)
	x3=x_tmp-5+ind(0)
	y3=y_tmp-5+ind(1)
	x3=x_tmp
	y3=y_tmp
	OPLOT,[x3],[y3],PSYM=8,SYMSIZE=1,COL=254
	PRINT,'Polaris    =[',x3,y3,']'

; Get position of Deneb (StarID=13)
	PRINT,'Click Deneb (in Cygnus)'
	CURSOR,x_tmp,y_tmp,/DEVICE
	SPAWN,'sleep 1'
	P=MAX(image(x_tmp-5:x_tmp+5,y_tmp-5:y_tmp+5),i)
	ind=ARRAY_INDICES(image(x_tmp-5:x_tmp+5,y_tmp-5:y_tmp+5),i)
	x4=x_tmp-5+ind(0)
	y4=y_tmp-5+ind(1)
	x4=x_tmp
	y4=y_tmp
	OPLOT,[x4],[y4],PSYM=8,SYMSIZE=1,COL=254
	PRINT,'Deneb    =[',x4,y4,']'

; Center of the image (identified with CURSOR,x0,y0,/DEVICE)
	x_cnt=xsize/2
	y_cnt=ysize/2

; Try Capera + Vega
	PRINT,' '
	PRINT,'Capera + Vega'
	params=date+hh+'.'+mmss+' '+ $
		STRCOMPRESS(STRING(x_cnt))+' '+STRCOMPRESS(STRING(y_cnt))+' '+ $
		'4 '+STRCOMPRESS(STRING(x1))+' '+STRCOMPRESS(STRING(y1))+' '+ $
		'11 '+STRCOMPRESS(STRING(x2))+' '+STRCOMPRESS(STRING(y2))
	SPAWN,'/radar01/work/Aurora/Starmap/zenith2000_tjornes.sh '+params

; Try Capera + Deneb
	PRINT,' '
	PRINT,'Capera + Deneb'
	params=date+hh+'.'+mmss+' '+ $
		STRCOMPRESS(STRING(x_cnt))+' '+STRCOMPRESS(STRING(y_cnt))+' '+ $
		'4 '+STRCOMPRESS(STRING(x1))+' '+STRCOMPRESS(STRING(y1))+' '+ $
		'13 '+STRCOMPRESS(STRING(x4))+' '+STRCOMPRESS(STRING(y4))
	SPAWN,'/radar01/work/Aurora/Starmap/zenith2000_tjornes.sh '+params

; Try Polaris + Deneb
	PRINT,' '
	PRINT,'Polaris + Deneb'
	params=date+hh+'.'+mmss+' '+ $
		STRCOMPRESS(STRING(x_cnt))+' '+STRCOMPRESS(STRING(y_cnt))+' '+ $
		'1 '+STRCOMPRESS(STRING(x3))+' '+STRCOMPRESS(STRING(y3))+' '+ $
		'13 '+STRCOMPRESS(STRING(x4))+' '+STRCOMPRESS(STRING(y4))
	SPAWN,'/radar01/work/Aurora/Starmap/zenith2000_tjornes.sh '+params

	END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	DEFINE_ZENITH_SYOWA
;

	PRO define_zenith_syowa,preview=preview,confirm=confirm

; Image used for defining zenith
	fname='./SYO_Watec_11-09-30_22-00-08-27.jpg'

; Time info for input to zenith2000
	date='20110930'
	hh='22'
	mmss='0008'

; Read image and identify size
	READ_JPEG,fname,image
	image_size=SIZE(image,/DIM)
	xsize=image_size(0)
	ysize=image_size(1)
	PRINT,xsize,ysize
	WINDOW,XSIZE=xsize,YSIZE=ysize

; Plot image
	LOADCT,0
	;LOADCT,39
	;cut_col_tab
	image=BYTSCL(image,MAX=50.0,MIN=0.0,TOP=252.0)
	TVSCL,image
	!P.POSITION=[0,0,1,1]
	PLOT,[0],[0],XRANGE=[0,xsize],YRANGE=[0,ysize],/XSTYLE,/YSTYLE
	psym8

	IF KEYWORD_SET(preview) THEN RETURN

	IF KEYWORD_SET(confirm) THEN BEGIN

		X0=318
		Y0=248
		A=123.0
		ang=-106.0

		OPLOT,[X0],[Y0],PSYM=8,COL=254,SYMSIZE=2
		OPLOT,X0+[-500,500],[Y0,Y0],LINE=1
		OPLOT,[X0,X0],Y0+[-500,500],LINE=1
		dim=A*!PI*0.5
		x_cir=FLTARR(361)
		y_cir=FLTARR(361)
		FOR i=0,360 DO BEGIN
			x_cir(i)=X0+dim*COS(i*!DTOR)
			y_cir(i)=Y0+dim*SIN(i*!DTOR)
		ENDFOR	
		OPLOT,x_cir,y_cir

		RETURN
	ENDIF


; Get position of Sirius (StarID=6)
	PRINT,'Click Sirius'
	CURSOR,x_tmp,y_tmp,/DEVICE
	SPAWN,'sleep 1'
	P=MAX(image(x_tmp-5:x_tmp+5,y_tmp-5:y_tmp+5),i)
	ind=ARRAY_INDICES(image(x_tmp-5:x_tmp+5,y_tmp-5:y_tmp+5),i)
	x1=x_tmp-5+ind(0)
	y1=y_tmp-5+ind(1)
	x1=x_tmp
	y1=y_tmp
	OPLOT,[x1],[y1],PSYM=8,SYMSIZE=1,COL=254
	PRINT,'Sirius  =[',x1,y1,']'

; Get position of Fomalhaut (StarID=14)
	PRINT,'Click Fomalhaut'
	CURSOR,x_tmp,y_tmp,/DEVICE
	SPAWN,'sleep 1'
	P=MAX(image(x_tmp-5:x_tmp+5,y_tmp-5:y_tmp+5),i)
	ind=ARRAY_INDICES(image(x_tmp-5:x_tmp+5,y_tmp-5:y_tmp+5),i)
	x2=x_tmp-5+ind(0)
	y2=y_tmp-5+ind(1)
	x2=x_tmp
	y2=y_tmp
	OPLOT,[x2],[y2],PSYM=8,SYMSIZE=1,COL=254
	PRINT,'Fomalhaut    =[',x2,y2,']'

; Get position of Canopus (StarID=16)
	PRINT,'Click Canopus'
	CURSOR,x_tmp,y_tmp,/DEVICE
	SPAWN,'sleep 1'
	P=MAX(image(x_tmp-5:x_tmp+5,y_tmp-5:y_tmp+5),i)
	ind=ARRAY_INDICES(image(x_tmp-5:x_tmp+5,y_tmp-5:y_tmp+5),i)
	x3=x_tmp-5+ind(0)
	y3=y_tmp-5+ind(1)
	x3=x_tmp
	y3=y_tmp
	OPLOT,[x3],[y3],PSYM=8,SYMSIZE=1,COL=254
	PRINT,'Canopus    =[',x3,y3,']'

; Get position of Acrux (StarID=17)
	PRINT,'Click Acrux' 
	CURSOR,x_tmp,y_tmp,/DEVICE
	SPAWN,'sleep 1'
	P=MAX(image(x_tmp-5:x_tmp+5,y_tmp-5:y_tmp+5),i)
	ind=ARRAY_INDICES(image(x_tmp-5:x_tmp+5,y_tmp-5:y_tmp+5),i)
	x4=x_tmp-5+ind(0)
	y4=y_tmp-5+ind(1)
	x4=x_tmp
	y4=y_tmp
	OPLOT,[x4],[y4],PSYM=8,SYMSIZE=1,COL=254
	PRINT,'Acrux    =[',x4,y4,']'

; Center of the image (identified with CURSOR,x0,y0,/DEVICE)
	x_cnt=xsize/2
	y_cnt=ysize/2

; Try Sirius + Fomalhaut
	PRINT,' '
	PRINT,'Sirius + Fomalhaut'
	params=date+hh+'.'+mmss+' '+ $
		STRCOMPRESS(STRING(x_cnt))+' '+STRCOMPRESS(STRING(y_cnt))+' '+ $
		'6 '+STRCOMPRESS(STRING(x1))+' '+STRCOMPRESS(STRING(y1))+' '+ $
		'14 '+STRCOMPRESS(STRING(x2))+' '+STRCOMPRESS(STRING(y2))
	SPAWN,'/radar01/work/Aurora/Starmap/zenith2000_syowa.sh '+params

; Try Sirius + Acrux
	PRINT,' '
	PRINT,'Sirius + Acrux'
	params=date+hh+'.'+mmss+' '+ $
		STRCOMPRESS(STRING(x_cnt))+' '+STRCOMPRESS(STRING(y_cnt))+' '+ $
		'6 '+STRCOMPRESS(STRING(x1))+' '+STRCOMPRESS(STRING(y1))+' '+ $
		'17 '+STRCOMPRESS(STRING(x4))+' '+STRCOMPRESS(STRING(y4))
	SPAWN,'/radar01/work/Aurora/Starmap/zenith2000_syowa.sh '+params

; Try Canopus + Fomalhaut
	PRINT,' '
	PRINT,'Canopus + Fomalhaut'
	params=date+hh+'.'+mmss+' '+ $
		STRCOMPRESS(STRING(x_cnt))+' '+STRCOMPRESS(STRING(y_cnt))+' '+ $
		'16 '+STRCOMPRESS(STRING(x3))+' '+STRCOMPRESS(STRING(y3))+' '+ $
		'14 '+STRCOMPRESS(STRING(x2))+' '+STRCOMPRESS(STRING(y2))
	SPAWN,'/radar01/work/Aurora/Starmap/zenith2000_syowa.sh '+params

; Try Acrux + Fomalhaut
	PRINT,' '
	PRINT,'Acrux + Fomalhaut'
	params=date+hh+'.'+mmss+' '+ $
		STRCOMPRESS(STRING(x_cnt))+' '+STRCOMPRESS(STRING(y_cnt))+' '+ $
		'14 '+STRCOMPRESS(STRING(x2))+' '+STRCOMPRESS(STRING(y2))+' '+ $
		'17 '+STRCOMPRESS(STRING(x4))+' '+STRCOMPRESS(STRING(y4))
	SPAWN,'/radar01/work/Aurora/Starmap/zenith2000_syowa.sh '+params

	END


























































        



















        ;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       RADAR_TO_GROUND_RANGE2
;
; NOTE:
;
;       Function for converting the radar range to the ground range
;

        FUNCTION radar_to_ground_range2,radar_range,altitude

; Earth radius
        Re=6370.0

; Just using "YOGEN-TEIRI"
        sin_el=((altitude+Re)^2-Re^2-radar_range^2)/(2*Re*radar_range)
        ground_range=Re*ATAN(SQRT(1-sin_el^2)/(Re/radar_range+sin_el))
        
        RETURN,ground_range

        END

;-------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       FIND_CELL2
;
; NOTE:
;
;       Function for finding radar range gates without rbpos
;

        FUNCTION find_cell2,radarlat,radarlon,nobeams,boresite,beamsep,gatelen,firstrange,beamno,gateno

; Earth radius
        Re=6370.0

        lon=radarlon*!DTOR
        colat=(90-radarlat)*!DTOR
        azimuth=((((boresite+beamsep*(beamno-0.5*(nobeams-1)))+900) MOD 360)-180)*!DTOR
        range=firstrange+gateno*gatelen
        
; Uncomment this line to take account of ground range being less than radar range
        range=radar_to_ground_range(range,110.0)

; Just using "KYUMEN-SANKAKU"
        ccolat=ACOS(COS(range/Re)*COS(colat)+SIN(range/Re)*SIN(colat)*COS(azimuth))
        clon=ACOS(MIN([(COS(range/Re)-COS(ccolat)*COS(colat))/(SIN(ccolat)*SIN(colat)), $
                1.0]))
        clat=90-ccolat/!DTOR
        IF azimuth GT 0 THEN clon=(lon+clon)/!DTOR ELSE clon=(lon-clon)/!DTOR
        
        RETURN,[clat,clon]

        END

;-------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       RANGE_TO_COORDS2
;
; NOTES:
;
;       Given a start "lat" and "lon" find the position of the point "range" in km away
;       at an azimuth of "azimuth" from north.
;

        FUNCTION range_to_coords2,lat,lon,azimuth,range

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
;       ANG_TO_LOC2
;
; NOTE:
;
;       Given an ground position (geog coords) and azimuth and the zenith of
;       observations and the assumed emission altitude return the geographic
;       lat and lon (i.e. [ele,azm,alt] -> [lat,lon])
;

        FUNCTION ang_to_loc2,lat,lon,azimuth,zenith,altitude

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

;-------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       GET_MAPPED_POINT2
;
; PURPOSE:
;
;       Calculate ionospheric point in lat and lon for given ele and azm
;       assuming altitude of the ionosphere.
;

        FUNCTION get_mapped_point2,ele,azm,alt

        lat_tjo=66.20
        lon_tjo=342.88

        lat_and_lon=ang_to_loc(lat_tjo,lon_tjo,azm,90-ele,alt)

        RETURN,lat_and_lon

        END

;-------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	GET_GROUND_DISTANCE_FROM_TJORNES2
; 
;       Get distance (on the surface of the Earth) in km from Tjornes to given point.
;

	FUNCTION get_ground_distance_from_tjornes2,lat_giv,lon_giv

; Location of Tjornes
	lat_tjo= 66.20
	lon_tjo=342.88

; Radius of the Earth
        Re=6370.

; Co-latitude from pole to Tjornes
        colata=(90.0-lat_tjo)/180.*!PI

; Co-latitude from pole to point B
        colatb=(90.0-lat_giv)/180.*!PI

; Angle between the lines from pole to A and from pole to B
        inpro=ABS(lon_giv-lon_tjo)/180.0*!PI

; Distance between point A and point B
        ground_dist=Re*ACOS(COS(colata)*COS(colatb)+SIN(colata)*SIN(colatb)*COS(inpro))

; Return the results
        RETURN,ground_dist

	END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	GET_DISTANCE_FROM_TJORNES2
;

	FUNCTION get_distance_from_tjornes2,lat_giv,lon_giv,height

; Radius of the Earth
        Re=6370.

; Get ground distance from Tjornes
	ground_distance=get_ground_distance_from_tjornes(lat_giv,lon_giv)

; Get offset angle on the big circle including Tjornes and given point
	ang_offset=360.0*ground_distance/(2*!PI*Re)

; Get straigh line distance
	distance=SQRT(Re^2+(Re+height)^2-2*Re*(Re+height)*COS(!PI*ang_offset/180.0))

	RETURN,distance

	END

;-------------------------------------------------------------------------------------------------------------------
; NAME:
;
; 	GET_ZENITH_ANGLE_FROM_TJORNES2
;

	FUNCTION get_zenith_angle_from_tjornes2,lat_giv,lon_giv,height

; Radius of the Earth
        Re=6370.

; Get straight line distance from Tjornes to the point specified
	distance=get_distance_from_tjornes(lat_giv,lon_giv,height)

; Get elevation angle
	angle1=ACOS((distance^2-height^2-2*Re*height)/(2*distance*Re))*180.0/!PI
	zenith_angle=180.0-angle1

	RETURN,zenith_angle

	END

;-------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	GET_AZIMUTH_FROM_TJORNES2
;

	FUNCTION get_azimuth_from_tjornes2,lat_giv,lon_giv

; Location of Tjornes
	lat_tjo= 66.20
	lon_tjo=-17.12

; Radius of the Earth
        Re=6370.

; Angle between the lines from pole to A and from pole to B
;        small_c=!PI*(lon_giv-lon_tjo)/180.0
        large_c=!PI*(lon_giv-lon_tjo)/180.0

; Get colatitude of Tjornes and given point
	colat_tjo=90.0-lat_tjo
	colat_giv=90.0-lat_giv
	small_b=!PI*colat_tjo/180.0
	small_a=!PI*colat_giv/180.0

; Get ground distance from Tjornes
	cos_small_c=COS(small_a)*COS(small_b)+SIN(small_a)*SIN(small_b)*COS(large_c)
	small_c=ACOS(cos_small_c)

; Get azimuth angle
	cos_large_a=(COS(small_a)-COS(small_b)*COS(small_c))/(SIN(small_b)*SIN(small_c))
	azimuth=180.0*ACOS(cos_large_a)/!PI
	IF lon_giv LT lon_tjo THEN azimuth=-azimuth
	IF azimuth LT 0 THEN azimuth=azimuth+360.0

	RETURN,azimuth

	END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	GET_PIXEL2
;
; PURPOSE:
;
;	Get x and y of pixel of the specified zenith and azimuth angle.
;

	FUNCTION get_pixel2,zan,azm

	COMMON tjornes

; Values of fish-eye lens define with define_zenith procedure
	X0=siz_wat2/2
	Y0=siz_wat2/2

	azm2=360-azm
	azm3=azm2 MOD 360
	azm4=azm3-rog_wat2
	x=FIX(X0+(0.5*siz_wat2*zan/75.0)*SIN(azm4*!PI/180.0))
	y=FIX(Y0+(0.5*siz_wat2*zan/75.0)*COS(azm4*!PI/180.0))
	
	RETURN,[x,y]

	END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	GET_ZAN_AZM2
;
; PURPOSE:
;
;	Get zenith and azimuth angle of specified x and y pixel on the image.
;

	FUNCTION get_zan_azm2,x,y

	COMMON tjornes

; Values of fish-eye lens define with define_zenith procedure
	X0=siz_wat2/2
	Y0=siz_wat2/2

; Displacement from zenith
	dx=FLOAT(x-X0)
	dy=FLOAT(y-Y0)

; Convert (dx,dy) into zenith and azimuth
	zan=75.0*SQRT(dx^2+dy^2)/FLOAT(siz_wat2*0.5)
	azm0=90.0-180.0*ATAN(dy,dx)/!PI
	azm1=azm0+rog_wat2
	azm2=360.0-azm1
	IF azm2 LT   0.0 THEN azm2=azm2+360.0
	IF azm2 GT 360.0 THEN azm2=azm2-360.0

	RETURN,[zan,azm2]
	;RETURN,[ROUND(zan),ROUND(azm2)]

	END

;-------------------------------------------------------------------------------------------------------------------
; NAME
;
;	FILE_WAT2
;

	PRO file_wat2,path=path,retry=retry

	COMMON tjornes

; Imomusi data directory
	IF NOT KEYWORD_SET(path) THEN path='/raid/data/Aurora/Tjorwat/20110930_htr_tks/tjo/'
	;IF NOT KEYWORD_SET(path) THEN path='/raid/data/Aurora/Tjorwat/20110930_htr_all/'

; Search the data files
        found_file=FILE_SEARCH(path+'*.jpg',COUNT=num_found_file)
        num_wat2=num_found_file
       

; Values of fish-eye lens define with define_zenith procedure (see alignment.txt in detail)
        X0=316         ; X of central pixel 320
        Y0=246         ; Y of central pixel 240
	a_v_wat2=156    ; A-value (diameter/3.1415)
	rog_wat2=15.18   ; GN is rotated to the East by this angle in degree
	rom_wat2=-14.56 ; Declination at Tjornes (taken from WDC) in 2009

; Grid adjustment
	cnt_wat=[X0,Y0]
	siz_wat2=FIX(a_v_wat2*!PI*75.0/90.0)            ; zenith angle limit of 75 deg
	IF siz_wat2 MOD 2 EQ 1 THEN siz_wat2=siz_wat2+1  ; siz_wat2 must be an even number

; Initialize arrays

        ;num_wat2=1000
	img_wat2=INTARR(num_wat2,siz_wat2+1,siz_wat2+1)
	keo_wat2=INTARR(5,num_wat2,siz_wat2+1)
	day_wat2=STRARR(num_wat2)
	yer_wat2=STRARR(num_wat2)
	yrs_wat2=STRARR(num_wat2)
	jdy_wat2=STRARR(num_wat2)
	tsc_wat2=LONARR(num_wat2)
	tim_wat2=STRARR(num_wat2)
	hms_wat2=STRARR(num_wat2)

; Read the data
	FOR i=0,num_wat2-1 DO BEGIN

                hh=STRMID(found_file(i),STRLEN(path)+17,2)
                IF hh EQ 00 THEN hh=STRING(24)
                mm=STRMID(found_file(i),STRLEN(path)+20,2)
                ss=STRMID(found_file(i),STRLEN(path)+23,2)
                tsc_wat2(i)=3600L*FIX(hh)+60L*FIX(mm)+FIX(ss)
                hms_wat2(i)=hh+mm+ss
                READ_JPEG,found_file(i),img_tmp
		PRINT,SIZE(img_tmp,/DIM)

; Get time info
                yr=STRMID(found_file(i),STRLEN(path)+08,2)
                mo=STRMID(found_file(i),STRLEN(path)+11,2)
                dy=STRMID(found_file(i),STRLEN(path)+14,2)
                yer_wat2(i)=2000+FIX(yr)
        	day_wat2(i)='20'+yr+mo+dy
                tim_wat2(i)=hh+mm+' '+ss+'s UT'
		jdy_wat2(i)=JULDAY(FIX(mo),FIX(dy),yer_wat2(i))-JULDAY(1,0,yer_wat2(i))
		yrs_wat2(i)=90000L*(jdy_wat2(i)-1)+tsc_wat2(i)
        	PRINT,'Scan'+STRING(i,FORMAT='(I5)')+': '+day_wat2(i)+' '+tim_wat2(i)+ $
			' ('+STRMID(found_file(i),STRLEN(path),32)+')'

; Get 2D image
                img_wat2(i,*,*)=img_tmp(cnt_wat(0)-siz_wat2/2: $
			cnt_wat(0)+siz_wat2/2,cnt_wat(1)-siz_wat2/2:cnt_wat(1)+siz_wat2/2)
                print,SIZE(img_wat2,/DIM)
                ;loadct,[0]
                ;TV,img_wat2

; Get keogram from GS to GN (Line ID 0)
		ang_for_keo=rog_wat2
		offset=FIX(TAN(ang_for_keo*!PI/180.0)*(siz_wat2/2))
		FOR k=0,siz_wat2 DO BEGIN
			x_keo=FIX(siz_wat2/2+offset-k*TAN(ang_for_keo*!PI/180.0))
			y_keo=k
			keo_wat2(0,i,k)=img_wat2(i,x_keo,y_keo)
		ENDFOR

; Get keogram from GE to GW (Line ID 1)
		ang_for_keo=rog_wat2
		offset=FIX(TAN(ang_for_keo*!PI/180.0)*(siz_wat2/2))
		FOR k=0,siz_wat2 DO BEGIN
			x_keo=k
			y_keo=FIX(siz_wat2/2-offset+k*TAN(ang_for_keo*!PI/180.0))
			keo_wat2(1,i,k)=img_wat2(i,x_keo,y_keo)
		ENDFOR

; Get keogram from MS to MN (Line ID 2)
		ang_for_keo=rog_wat2+rom_wat2
		offset=FIX(TAN(ang_for_keo*!PI/180.0)*(siz_wat2/2))
		FOR k=0,siz_wat2 DO BEGIN
			x_keo=FIX(siz_wat2/2+offset-k*TAN(ang_for_keo*!PI/180.0))
			y_keo=k
			keo_wat2(2,i,k)=img_wat2(i,x_keo,y_keo)
		ENDFOR

; Get keogram from ME to MW (Line ID 3)
		ang_for_keo=rog_wat2+rom_wat2
		offset=FIX(TAN(ang_for_keo*!PI/180.0)*(siz_wat2/2))
		FOR k=0,siz_wat2 DO BEGIN
			x_keo=k
			y_keo=FIX(siz_wat2/2-offset+k*TAN(ang_for_keo*!PI/180.0))
			keo_wat2(3,i,k)=img_wat2(i,x_keo,y_keo)
		ENDFOR

; Get keogram (2: beam 7 align)
                boresite_from_tjornes=32.0 ; NOTE: This is not the radar boresite! Azimuth of the beam 7 from TJO.
               ; ang_for_keo=rog_wat2+boresite_from_tjornes-3.24*0.5
		;offset=FIX(TAN(ang_for_keo*!PI/180.0)*(siz_wat2/2))
		;FOR k=0,siz_wat2 DO BEGIN
		;	x_keo=FIX(siz_wat2/2+offset-k*TAN(ang_for_keo*!PI/180.0))
		;	y_keo=k
		;	keo_wat2(4,i,k)=img_wat2(i,x_keo,y_keo)
	;	ENDFOR

	ENDFOR

; Summary of settings
	PRINT,''
	PRINT,'Summary of current settings'
	PRINT,'------------------------------------------'

; Current image
	go_wat2,0L

; Start and end time
	time_wat2,hms_wat2(0),hms_wat2(num_wat2-1)
	IF STRMID(hms_wat2(num_wat2-1),3,1) EQ '9' THEN time_wat2,hms_wat2(0),hms_wat2(num_wat2-1)

; Set scale
	set_scale_wat2,0,200
	PRINT,'------------------------------------------'
	PRINT,''

; Set default charsize
	!P.CHARSIZE=1.0

	END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	MAKE_FOV_WAT2
;
	
	PRO make_fov_wat2

	COMMON tjornes

	height=110.0

; Initialize array
	gla_wat2=FLTARR(siz_wat2+1,siz_wat2+1)
	glo_wat2=FLTARR(siz_wat2+1,siz_wat2+1)
	mla_wat2=FLTARR(siz_wat2+1,siz_wat2+1)
	mlo_wat2=FLTARR(siz_wat2+1,siz_wat2+1)
	zan_wat2=FLTARR(siz_wat2+1,siz_wat2+1)

	;PLOT,[0],[0],XRANGE=[330,360],YRANGE=[60,75]
	FOR xx=0,siz_wat2 DO BEGIN
		FOR yy=0,siz_wat2 DO BEGIN

			tmp=get_zan_azm(xx,yy)
			zan=tmp(0)
			azm=tmp(1)
			IF zan LT 0.0 OR zan GT 75.0 THEN BEGIN
				zan=99999.9
				azm=99999.9
			ENDIF
			IF FINITE(zan,/NAN) THEN BEGIN
				zan=99999.9
				azm=99999.9
			ENDIF
			zan_wat2(xx,yy)=zan
			IF zan NE 99999.9 AND azm NE 99999.9 THEN BEGIN
        			tmp=get_mapped_point(90-zan,azm,height)
				;!P.SYMSIZE=1
				;OPLOT,[tmp(1)],[tmp(0)],COL=252,SYMSIZE=0.2*!P.SYMSIZE,PSYM=1
				gla_wat2(xx,yy)=tmp(0)
				glo_wat2(xx,yy)=tmp(1)
				tmp=cnvcoord(gla_wat2(xx,yy),glo_wat2(xx,yy),height)
				mla_wat2(xx,yy)=tmp(0)
				mlo_wat2(xx,yy)=tmp(1)
			ENDIF ELSE BEGIN
				gla_wat2(xx,yy)=99999.9
				glo_wat2(xx,yy)=99999.9
				mla_wat2(xx,yy)=99999.9
				mlo_wat2(xx,yy)=99999.9
			ENDELSE

		ENDFOR
	ENDFOR

	END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	GO_WAT2
;
; EXAMPLE:
;
;       Go > go_wat2,10
;

	PRO go_wat2,img_to_go

	COMMON tjornes

; Set new current image
	now_wat2=0L
	now_wat2=img_to_go

; Verbose
        current_hms=STRING(hms_wat2(now_wat2),FORMAT='(I6.6)')
        PRINT,'Time of current image: '+current_hms+' UT'

	END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	NEXT_IMAGE2
;

	PRO next_image2

	COMMON tjornes

	go_wat,now_wat2+1

	END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       GO_TIME_WAT2
;
; EXAMPLE:
;
;       Go > go_time_wat,025800
;

        PRO go_time_wat2,image_time

        COMMON tjornes

; Finding image
        found_image=WHERE(hms_wat2 EQ STRING(image_time,FORMAT='(I6.6)'),num_found_image)
        print,found_image
; Jump image
        IF num_found_image EQ 0 THEN BEGIN
                PRINT,'scan out of range!'
        ENDIF ELSE BEGIN
                go_wat2,found_image
        ENDELSE

	END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	TIME_WAT2
;
; EXAMPLE:
;
;	Go > time_wat,030000,031000
;
; NOTE:
;
;	Time should be specified in 6-digits (hhmmss).
;

        PRO time_wat2,stime,etime

        COMMON tjornes

; Store the data
        stm_wat2=stime
        etm_wat2=etime

; Verbose
	PRINT,'Start time: '+STRMID(STRING(stime,FORMAT='(I6.6)'),0,4)+' '+ $
		STRMID(STRING(stime,FORMAT='(I6.6)'),4,2)+'s UT'
	PRINT,'End   time: '+STRMID(STRING(etime,FORMAT='(I6.6)'),0,4)+' '+ $
		STRMID(STRING(etime,FORMAT='(I6.6)'),4,2)+'s UT'

; Change time info on Go
        time,stime/100L,etime/100L

        END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       SET_SCALE_WAT2
;
; EXAMPLE:
;
;	Go > set_scale_wat,0,300
;

        PRO set_scale_wat2,min_val,max_val,quiet=quiet

        COMMON tjornes

; Store the data in common block
        min_wat2=min_val
        max_wat2=max_val

; Verbose
	IF NOT KEYWORD_SET(quiet) THEN BEGIN
		PRINT,'Min intensity: '+STRING(min_val,FORMAT='(I4)')
		PRINT,'Max intensity: '+STRING(max_val,FORMAT='(I4)')
	ENDIF

	END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       PLOT_WAT2
;
; EXAMPLE:
;
;	Go > plot_wat,2,3
;

        PRO plot_wat2,xmaps,ymaps,skip=skip,ctab=ctab,_extra=e

        COMMON tjornes
	COMMON ps_info

; Keep current charsize
	charsize_keep=!P.CHARSIZE

; Refresh window
	clear_page

; Setting panel positions
        IF N_PARAMS() NE 2 THEN BEGIN & xmaps=1 & ymaps=1 & ENDIF
        IF NOT KEYWORD_SET(skip) THEN skip=1

; Change charsize
	!P.CHARSIZE=MAX([MIN([1.0/xmaps,1.0/ymaps]),0.24])
	IF psinfo.open EQ 0 THEN !P.CHARSIZE=2.0*!P.CHARSIZE
	IF psinfo.open EQ 1 THEN !P.CHARSIZE=1.5*!P.CHARSIZE

; Plot panels
        FOR ymap=0,ymaps-1 DO BEGIN
                FOR xmap=0,xmaps-1 DO BEGIN

                        plot_wat_panel2,xmaps,ymaps,xmap,ymap,ctab=ctab,_extra=e
                        now_wat2+=skip
                ENDFOR
        ENDFOR
	IF xmaps EQ 1 AND ymaps EQ 1 THEN now_wat2-=skip

; Plot colour bar
        plot_wat_colour_bar,ctab=ctab

; Restore old charsize
	!P.CHARSIZE=charsize_keep

        END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       PLOT_WAT_PANEL2
;
; NOTE:
;
;	1. if ctab is not specified, cut_col_tab is used for plot.
;	2. if cont is specified, the data will be plotted with CONTOUR. Otherwise, plotted with TV.
;	3. if irev is specified, the image is reversed in the east-west direction.
;

        PRO plot_wat_panel2,xmaps,ymaps,xmap,ymap,position=position,ctab=ctab,cont=cont,irev=irev, $
		myopic_data=myopic_data,myopic_fov=myopic_fov,force_char=force_char,white_grid=white_grid, $
		mnorth=mnorth,no_title=no_title,no_grid=no_grid,special_grid=special_grid,no_ew=no_ew

        COMMON colour_info
        COMMON tjornes

	COMMON data
	COMMON prefs

; Update colour info
        IF N_ELEMENTS(ctab) EQ 0 THEN BEGIN
                set_colour_table,/default
        ENDIF ELSE BEGIN
                set_colour_table,ctab
        ENDELSE
        update_colour_info

        
; Charsize
	IF KEYWORD_SET(force_char) THEN BEGIN
		old_charsize=!P.CHARSIZE
		!P.CHARSIZE=force_char
	ENDIF

; Define plotting region. Default is one panel on screen
        IF N_PARAMS() NE 4 THEN BEGIN & xmaps=1 & ymaps=1 & xmap=0 & ymap=0 & ENDIF
        IF NOT KEYWORD_SET(position) THEN BEGIN
                define_panel,xmaps,ymaps,xmap,ymap,/bar,/square,/no_charset
        	position=!P.POSITION
        ENDIF
        !P.POSITION=position

; Plot frame
        PLOT,[0],[0],XRANGE=[0,siz_wat2],YRANGE=[0,siz_wat2],/XSTYLE,/YSTYLE, $
                XTICKFORMAT='no_ticks',YTICKFORMAT='no_ticks',XTICKS=1,YTICKS=1

; Plot image with contour
        IF KEYWORD_SET(cont) THEN BEGIN
                c_lines=30
                levels_set=(max_wat2-min_wat2)*FINDGEN(c_lines)/FLOAT(c_lines)+min_wat2
                colour_set=252.0*FINDGEN(c_lines)/FLOAT(c_lines)
                colour_set(0)=0
                x_cont=INDGEN(siz_wat2+1)
                y_cont=INDGEN(siz_wat2+1)
		IF KEYWORD_SET(irev) THEN BEGIN
              		z_cont=REVERSE(REFORM(img_wat2(now_wat2,*,*)))
		ENDIF ELSE BEGIN
                	z_cont=REFORM(img_wat2(now_wat2,*,*))
		ENDELSE
                CONTOUR,z_cont,x_cont,y_cont,/FILL,LEVELS=levels_set, $
			C_COLORS=colour_set,/OVERPLOT
        ENDIF ELSE BEGIN

; Plot the data with TV
		IF KEYWORD_SET(irev) THEN BEGIN
                	img1=REVERSE(REFORM(img_wat2(now_wat2,*,*)))
		ENDIF ELSE BEGIN
                	img1=REFORM(img_wat2(now_wat2,*,*))
		ENDELSE
                img2=CONGRID(img1,FIX((position(2)-position(0))*!D.X_Size), $
			FIX((position(3)-position(1))*!D.Y_Size))
                img3=BYTSCL(img2,MAX=max_wat2,MIN=min_wat2,TOP=252.0)
                TV,img3,position(0)*!D.X_Size,position(1)*!D.Y_Size, $
                        XSIZE=(position(2)-position(0))*!D.X_Size, $
			YSIZE=(position(3)-position(1))*!D.Y_Size
        ENDELSE

; Plot info
	IF NOT KEYWORD_SET(no_title) THEN BEGIN
        ;	XYOUTS,0,siz_wat2+5,day_wat2(now_wat2),/DATA,COL=foreground,CHARSIZE=1.1    ;charsize kaetayo
        ;	XYOUTS,siz_wat2-1,siz_wat2+5,tim_wat2(now_wat2),/DATA,COL=foreground,ALI=1,CHARSIZE=1.3  
	ENDIF
        PRINT,'Scan'+STRCOMPRESS(now_wat2)+' '+day_wat2(now_wat2)+' '+tim_wat2(now_wat2)

; Plot the geographic grid (axis is inclined by rog_wat2)
	IF NOT KEYWORD_SET(no_grid) THEN BEGIN
		ofs_wid=FLOAT(siz_wat2/2)*TAN(rog_wat2*!PI/180.0)
		IF KEYWORD_SET(irev) THEN sign=-1 ELSE sign=+1
        	OPLOT,0.5*(siz_wat2)*[1,1]+sign*[ofs_wid,-ofs_wid],[0,siz_wat2],COL=white,THICK=1,LINE=1
        	OPLOT,[0,siz_wat2],0.5*(siz_wat2)*[1,1]-sign*[ofs_wid,-ofs_wid],COL=white,THICK=1,LINE=1
	ENDIF

; Plot the geomagnetic S-N (axis is inclined by rog_wat2+rom_wat2)
	ofs_wid=FLOAT(siz_wat2/2)*TAN((rog_wat2+rom_wat2)*!PI/180.0)
	IF KEYWORD_SET(irev) THEN sign=-1 ELSE sign=+1
        OPLOT,0.5*(siz_wat2)*[1,1]+sign*[ofs_wid,-ofs_wid],[0,siz_wat2],COL=white,THICK=1,LINE=2

; Plot direction
	IF NOT KEYWORD_SET(no_grid) THEN BEGIN
		IF NOT KEYWORD_SET(irev) THEN BEGIN
			XYOUTS,siz_wat2+5,siz_wat2/2+(siz_wat2/2)*TAN(rog_wat2*!DTOR)-5,'W',ALI=0.0,/DATA
			XYOUTS,-5,siz_wat2/2-(siz_wat2/2)*TAN(rog_wat2*!DTOR)-5,'E',ALI=1.0,/DATA
			XYOUTS,siz_wat2/2-(siz_wat2/2)*TAN(rog_wat2*!DTOR),siz_wat2+5,'N',ALI=0.5,/DATA
			XYOUTS,siz_wat2/2+(siz_wat2/2)*TAN(rog_wat2*!DTOR),-15,'S',ALI=0.5,/DATA
			XYOUTS,siz_wat2/2-(siz_wat2/2)*TAN((rog_wat2+rom_wat2)*!DTOR),siz_wat2+5,'MN',ALI=0.5,/DATA
		ENDIF ELSE BEGIN
			XYOUTS,-5,siz_wat2/2+(siz_wat2/2)*TAN(rog_wat2*!DTOR)-5,'W',ALI=1.0,/DATA
			XYOUTS,siz_wat2+5,siz_wat2/2-(siz_wat2/2)*TAN(rog_wat2*!DTOR)-5,'E',ALI=0.0,/DATA
			XYOUTS,siz_wat2/2+(siz_wat2/2)*TAN(rog_wat2*!DTOR),siz_wat2+5,'N',ALI=0.5,/DATA
			XYOUTS,siz_wat2/2-(siz_wat2/2)*TAN(rog_wat2*!DTOR),-15,'S',ALI=0.5,/DATA
			XYOUTS,siz_wat2/2+(siz_wat2/2)*TAN((rog_wat2+rom_wat2)*!DTOR),siz_wat2+5,'MN',ALI=0.5,/DATA
		ENDELSE
	ENDIF

; Mask edge and overlay zenith
	mask_fov2
	IF NOT KEYWORD_SET(no_grid) THEN overlay_zenith2
	IF KEYWORD_SET(special_grid) THEN overlay_zenith_special2

; Overplot SuperDARN stuff?
	IF KEYWORD_SET(myopic_fov) THEN overlay_myopic_fov2,irev=irev
	IF KEYWORD_SET(myopic_data) THEN overlay_myopic_data2,irev=irev

; Overplot frame again
        PLOT,[0],[0],XRANGE=[0,siz_wat2],YRANGE=[0,siz_wat2],/XSTYLE,/YSTYLE, $
                XTICKFORMAT='no_ticks',YTICKFORMAT='no_ticks',XTICKS=1,YTICKS=1
	IF NOT KEYWORD_SET(no_title) THEN BEGIN
		;XYOUTS,10,10,'Outermost circle is 75!9'+STRING("260B")+'!3 zenith angle',CHARSIZE=0.6*!P.CHARSIZE,COL=white
		IF KEYWORD_SET(irev) THEN $
			XYOUTS,siz_wat2+3,siz_wat2,'East-West Reversed',ORI=270,CHARSIZE=0.5*!P.CHARSIZE
	ENDIF

; Plot SuperDARN scan time
	IF KEYWORD_SET(myopic_data) THEN BEGIN
		scan_beams=WHERE(beam_scan EQ scan_no,no_scan_beams)
		pos=!P.POSITION
		XYOUTS,pos(0)+0.01,pos(1)+0.02,str_time(beam_time(scan_beams(0)),1)+' UT',/NORMAL,CHARSIZE=1.3
	ENDIF

; Restore old charsize
	IF KEYWORD_SET(force_char) THEN !P.CHARSIZE=old_charsize

	END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       OVERLAY_WAT_CONT2
;

        PRO overlay_wat_cont2,irev=irev

        COMMON colour_info
        COMMON tjornes

	COMMON data
	COMMON prefs

; Update colour info
        IF N_ELEMENTS(ctab) EQ 0 THEN BEGIN
                set_colour_table,/default
        ENDIF ELSE BEGIN
                set_colour_table,ctab
        ENDELSE
        update_colour_info

; Plot image with contour
        ;c_lines=30
        ;levels_set=(max_wat2-min_wat2)*FINDGEN(c_lines)/FLOAT(c_lines)+min_wat2
        ;colour_set=252.0*FINDGEN(c_lines)/FLOAT(c_lines)
        ;colour_set(0)=0
	levels_set=[80,160]
        x_cont=INDGEN(siz_wat2+1)
        y_cont=INDGEN(siz_wat2+1)
	IF KEYWORD_SET(irev) THEN BEGIN
      		z_cont=REVERSE(REFORM(img_wat2(now_wat2,*,*)))
	ENDIF ELSE BEGIN
        	z_cont=REFORM(img_wat2(now_wat2,*,*))
	ENDELSE
        CONTOUR,z_cont,x_cont,y_cont,LEVELS=levels_set,/OVERPLOT,COL=0,THICK=6
        CONTOUR,z_cont,x_cont,y_cont,LEVELS=levels_set,/OVERPLOT,COL=254,THICK=2

	END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	PLOT_WAT_KEO2
;

	PRO plot_wat_keo2,_extra=e

	COMMON tjornes
	COMMON ps_info

; Refresh window
	clear_page

; Check time info
	time_info=STRING(stm_wat2,FORMAT='(I6.6)')+' to '+STRING(etm_wat2,FORMAT='(I6.6)')+' UT'

; Keep current charsize
	charsize_keep=!P.CHARSIZE

; Change charsize
	IF psinfo.open EQ 0 THEN !P.CHARSIZE=1.4*!P.CHARSIZE
	IF psinfo.open EQ 1 THEN !P.CHARSIZE=0.7*!P.CHARSIZE

; Plot two keograms
	plot_wat_keo_panel2,1,3,0,0,line_id=2,_extra=e
	XYOUTS,!P.POSITION(0),!P.POSITION(3)+0.01,'Keogram from Mag South to Mag North', $
		/NORMAL
	XYOUTS,!P.POSITION(2),!P.POSITION(3)+0.01,time_info,/NORMAL,ALI=1.0
	plot_wat_keo_panel2,1,3,0,1,line_id=3,_extra=e
	XYOUTS,!P.POSITION(0),!P.POSITION(3)+0.01,'Keogram from Mag East to Mag West', $
		/NORMAL
	XYOUTS,!P.POSITION(2),!P.POSITION(3)+0.01,time_info,/NORMAL,ALI=1.0
	plot_wat_keo_panel2,1,3,0,2,line_id=4,_extra=e
	XYOUTS,!P.POSITION(0),!P.POSITION(3)+0.01,'Beam 7 aligned Keogram', $
		/NORMAL
	XYOUTS,!P.POSITION(2),!P.POSITION(3)+0.01,time_info,/NORMAL,ALI=1.0

; Plot colour bar
	plot_wat_colour_bar2,_extra=e

; Restore old charsize
	!P.CHARSIZE=charsize_keep

	END


;-------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       PLOT_WAT_RTI_PANEL2
;

        PRO plot_wat_rti_panel2,yrange=yrange,no_x=no_x,no_frame=no_frame,beam_plot=beam_plot,cont=cont

        COMMON tjornes

; Variables
        radar_lat=63.77
        radar_lon=-20.54
        radar_boresite=30.0
        gate_length=15.0
        first_range=180.0
        altitude=110.0

; Angle of Tjornes Watec alignment
        rot_ang=rog_wat2

; Size of the image
        wid=siz_wat2/2

; Find radar cell locations
        x_wat=FLTARR(16,75)
        y_wat=FLTARR(16,75)
        FOR beam=0,15 DO BEGIN
                FOR gate=0,74 DO BEGIN
                        fov_pos=find_cell2(radar_lat,radar_lon,16,radar_boresite, $
                                3.24,gate_length,first_range,beam,gate+0.5)
                        lat=fov_pos(0)
                        lon=fov_pos(1)
                        zen=get_zenith_angle_from_tjornes(lat,lon,altitude)
                        ;azm=get_azimuth_from_tjornes(lat,lon)+90.0-rot_ang
                        ;azm=get_azimuth_from_tjornes(lat,lon)+90.0
                        azm=get_azimuth_from_tjornes(lat,lon)
                        IF azm GE 360.0 THEN azm=azm-360.0
                        tmp_pix=get_pixel(zen,azm)
                        x_wat(beam,gate)=tmp_pix(0)
                        y_wat(beam,gate)=tmp_pix(1)
                        ;x_wat(beam,gate)=wid+zen*COS(!PI*azm/180.0)*wid/75.0
                        ;y_wat(beam,gate)=wid+zen*SIN(!PI*azm/180.0)*wid/75.0
                ENDFOR
        ENDFOR

; Extract optical data from 2D array
        rti_wat2=INTARR(num_wat2,50)
        IF NOT KEYWORD_SET(beam_plot) THEN beam_plot=7
        FOR i=0,num_wat2-1 DO BEGIN

                FOR gate=0,49 DO BEGIN
                        IF FIX(x_wat(beam_plot,gate)) LT 367 AND FIX(y_wat(beam_plot,gate)) LT 367 THEN BEGIN

                                rti_wat2(i,gate)=img_wat2(i,FIX(x_wat(beam_plot,gate)),FIX(y_wat(beam_plot,gate)))
                        ENDIF
                ENDFOR

        ENDFOR

; Plot panel
        IF NOT KEYWORD_SET(yrange) THEN yrange=[10,40]
        IF NOT KEYWORD_SET(no_frame) THEN BEGIN
                plot_time_frame,position=!P.POSITION,yrange=yrange,ytitle='Range Gate',no_x=no_x
        ENDIF

; Plot the data
	IF KEYWORD_SET(cont) THEN BEGIN
		x_cont=tsc_wat2(0:num_wat2-1)/3600.0
		y_cont=FINDGEN(50)+11.0
        	;contour_lines=5
        	;levels_set=(max_wat2-min_wat2)*FINDGEN(contour_lines)/FLOAT(contour_lines)+min_wat2
        	;colour_set=!D.TABLE_SIZE*INDGEN(contour_lines)/contour_lines
		levels_set=[100,200]
                CONTOUR,rti_wat2,x_cont,y_cont,/OVERPLOT,LEVELS=levels_set,NOCLIP=0,COL=254,THICK=6
                CONTOUR,rti_wat2,x_cont,y_cont,/OVERPLOT,LEVELS=levels_set,NOCLIP=0,COL=0,THICK=3
                ;CONTOUR,rti_wat2,x_cont,y_cont,/OVERPLOT,/FILL,LEVELS=levels_set, $
                ;       	C_COLORS=colour_set,NOCLIP=0
	ENDIF ELSE BEGIN
        	FOR i=0,num_wat2-2 DO BEGIN
                	FOR j=0,49 DO BEGIN

                	x_poly=[tsc_wat2(i),tsc_wat2(i+1),tsc_wat2(i+1),tsc_wat2(i)]/3600.0
                	y_poly=[j,j,j+1,j+1]+11

                	col=252.0*(rti_wat2(i,j)-min_wat2)/FLOAT(max_wat2-min_wat2)
                	IF col LE   0 THEN col=1
                	IF col GE 252 THEN col=252
                	POLYFILL,x_poly,y_poly,COL=col,NOCLIP=0
	
                	ENDFOR
        	ENDFOR
	ENDELSE

; Overplot frame
        IF NOT KEYWORD_SET(no_frame) THEN $
                plot_time_frame,position=!P.POSITION,yrange=yrange,ytitle='Range Gate',no_x=no_x

        END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       PLOT_WAT_POLAR2
;
; EXAMPLE:
;
;       Go > plot_wat_polar,2,3
;

        PRO plot_wat_polar,xmaps,ymaps,skip=skip,ctab=ctab,_extra=e

        COMMON tjornes
        COMMON ps_info

; Keep current charsize
        charsize_keep=!P.CHARSIZE

; Refresh window
        clear_page

; Setting panel positions
        IF N_PARAMS() NE 2 THEN BEGIN & xmaps=1 & ymaps=1 & ENDIF
        IF NOT KEYWORD_SET(skip) THEN skip=1

; Change charsize
        !P.CHARSIZE=MAX([MIN([1.0/xmaps,1.0/ymaps]),0.24])
        IF psinfo.open EQ 0 THEN !P.CHARSIZE=2.0*!P.CHARSIZE
        IF psinfo.open EQ 1 THEN !P.CHARSIZE=1.5*!P.CHARSIZE

; Plot panels
        FOR ymap=0,ymaps-1 DO BEGIN
                FOR xmap=0,xmaps-1 DO BEGIN

                        plot_wat_polar_panel,xmaps,ymaps,xmap,ymap,ctab=ctab,_extra=e
                        now_wat2+=skip
                ENDFOR
        ENDFOR
        IF xmaps EQ 1 AND ymaps EQ 1 THEN now_wat2-=skip

; Plot colour bar
        plot_wat_colour_bar,ctab=ctab

; Restore old charsize
        !P.CHARSIZE=charsize_keep

        END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	PLOT_WAT_POLAR_PANEL2
;
	
	PRO plot_wat_polar_panel2,xmaps,ymaps,xmap,ymap,xrange=xrange,yrange=yrange,clip=clip, $
		position=position,coast=coast,no_data=no_data,ctab=ctab,force_char=force_char, $
		no_title=no_title

	COMMON tjornes

        COMMON data
        COMMON prefs

; Lat and lon of Tjornes
	lat_tjo=66.20
	lon_tjo=342.88

; Update colour info
        IF N_ELEMENTS(ctab) EQ 0 THEN BEGIN
                set_colour_table,/default
        ENDIF ELSE BEGIN
                set_colour_table,ctab
        ENDELSE
        update_colour_info

; Charsize
        IF KEYWORD_SET(force_char) THEN BEGIN
                old_charsize=!P.CHARSIZE
                !P.CHARSIZE=force_char
        ENDIF

; Define plotting region. Default is one panel on screen
        IF N_PARAMS() NE 4 THEN BEGIN & xmaps=1 & ymaps=1 & xmap=0 & ymap=0 & ENDIF
        IF NOT KEYWORD_SET(position) THEN BEGIN
                define_panel,xmaps,ymaps,xmap,ymap,/bar,/square,/no_charset
                position=!P.POSITION
        ENDIF
        !P.POSITION=position

; Define plot area
        IF NOT KEYWORD_SET(xrange) OR NOT KEYWORD_SET(yrange) THEN BEGIN
                xrange=[-31,31] & yrange=[-31,31]
                IF KEYWORD_SET(clip) THEN BEGIN
                        pos=cnvcoord(lat_tjo,lon_tjo,1)
                        mla_tjo=pos(0) & mlo_tjo=pos(1)
                        x1= abs(90-mla_tjo)*SIN(mlt(yer_wat2(now_wat2),yrs_wat2(now_wat2),mlo_tjo)*!PI/12)
                        y1=-abs(90-mla_tjo)*COS(mlt(yer_wat2(now_wat2),yrs_wat2(now_wat2),mlo_tjo)*!Pi/12)
			xyrange=10.0
                        xrange=[x1-xyrange*0.5,x1+xyrange*0.5]
                        yrange=[y1-xyrange*0.5,y1+xyrange*0.5]
                ENDIF
        ENDIF

; Plot frame
        plot_polar_frame,position=!P.POSITION, $
                xrange=xrange,yrange=yrange,black_back=black_back,ticklen=ticklen,no_frame=no_frame

; Plot image as contours
        tmp_bin=WHERE(gla_wat2 NE 99999.9 AND gla_wat2 NE 99999.9 AND zan_wat2 LE 75.0,no_tmp_bin)
        cont=MAKE_ARRAY(4,no_tmp_bin,/FLOAT,VALUE=0.0)
        num_cnt=0L
        FOR xxx=0,siz_wat2 DO BEGIN
                FOR yyy=0,siz_wat2 DO BEGIN
                        IF gla_wat2(xxx,yyy) NE 99999.9 AND glo_wat2(xxx,yyy) NE 99999.9 AND $
                                zan_wat2(xxx,yyy) LE 75.0 THEN BEGIN

                                cont(2,num_cnt)=FLOAT(img_wat2(now_wat2,xxx,yyy))
                                cont(3,num_cnt)=zan_wat2(xxx,yyy)
                                IF cont(2,num_cnt) GT max_wat2 THEN cont(2,num_cnt)=max_wat2
                                IF cont(2,num_cnt) LT min_wat2 THEN cont(2,num_cnt)=min_wat2

                                cont(0,num_cnt)= ABS(90-mla_wat2(xxx,yyy)) $
                                        *SIN(mlt(yer_wat2(now_wat2),yrs_wat2(now_wat2),mlo_wat2(xxx,yyy))*!PI/12)
                                cont(1,num_cnt)=-ABS(90-mla_wat2(xxx,yyy)) $
                                        *COS(mlt(yer_wat2(now_wat2),yrs_wat2(now_wat2),mlo_wat2(xxx,yyy))*!PI/12)

                                num_cnt++
                        ENDIF
                ENDFOR
        ENDFOR
        contour_lines=30
        levels_set=(max_wat2-min_wat2)*FINDGEN(contour_lines)/FLOAT(contour_lines)+min_wat2
        levels_lin=(max_wat2-min_wat2)*FINDGEN(6)/FLOAT(5)+min_wat2
        colour_set=!D.TABLE_SIZE*INDGEN(contour_lines)/contour_lines
        IF NOT KEYWORD_SET(no_data) THEN BEGIN
                CONTOUR,cont(2,*),cont(0,*),cont(1,*),/OVERPLOT,/FILL,LEVELS=levels_set, $
                        C_COLORS=colour_set,/IRREGULAR,NOCLIP=0
        ENDIF

; Plot frame
        plot_polar_frame,position=!P.POSITION, $
                xrange=xrange,yrange=yrange,black_back=black_back,ticklen=ticklen,no_frame=no_frame

; Overlay coast
        IF KEYWORD_SET(coast) THEN $
                overlay_polar_coast,force_year=yer_wat2(now_wat2),force_secs=yrs_wat2(now_wat2),col=253

; Plot info
        IF NOT KEYWORD_SET(no_title) THEN BEGIN
		pos=!P.POSITION
                XYOUTS,pos(0),pos(3)+0.01,day_wat2(now_wat2),/NORMAL,COL=foreground
                XYOUTS,pos(2),pos(3)+0.01,tim_wat2(now_wat2),/NORMAL,COL=foreground,ALI=1
        ENDIF
        PRINT,'Scan'+STRCOMPRESS(now_wat2)+' '+day_wat2(now_wat2)+' '+tim_wat2(now_wat2)

; Restore old charsize
        IF KEYWORD_SET(force_char) THEN !P.CHARSIZE=old_charsize

	END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	PLOT_WAT_LINE2
;

	PRO plot_wat_line2,_extra=e,click=click

; Refresh page if click keyword is not specified
	IF NOT KEYWORD_SET(click) THEN BEGIN
		click=0
		clear_page
	ENDIF

; Plot only one panel
	plot_wat_line_panel2,_extra=e,click=click

	END



;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	PLOT_WAT_COLOUR_BAR2
;

	PRO plot_wat_colour_bar2,ymaps,ymap,position=position,ctab=ctab,_extra=e,leg_pos=leg_pos

	COMMON tjornes

	COMMON prefs

; Update colour info
        IF N_ELEMENTS(ctab) EQ 0 THEN BEGIN
                set_colour_table,/default
        ENDIF ELSE BEGIN
                set_colour_table,ctab
        ENDELSE
        update_colour_info

; Setup scale
	old_minval=minval
	old_maxval=maxval
	set_scale,min_wat2,max_wat2,/quiet

; Reverse issue
	old_param=0
	IF parameter EQ 'vel' THEN BEGIN
		old_param=1
		pwr_l
		set_scale,min_wat2,max_wat2,/quiet
	ENDIF

; Plot colour bar
	IF KEYWORD_SET(position) THEN BEGIN
		plot_colour_bar,position=position,_extra=e, $
			legend='Optical Intensity!C!DArbitrary Unit!N',leg_pos=leg_pos,/no_gnd
			;legend='Optical Intensity (Arbitrary Unit)',leg_pos=leg_pos,/no_gnd
	ENDIF ELSE BEGIN
		IF N_PARAMS() NE 2 THEN BEGIN
			ymaps=1 & ymap=0
		ENDIF
		plot_colour_bar,ymaps,ymap,_extra=e, $
			legend='Optical Intensity!C!DArbitrary Unit!N',leg_pos=leg_pos,/no_gnd
			;legend='Optical Intensity (Arbitrary Unit)',leg_pos=leg_pos,/no_gnd
	ENDELSE

; Restore parameter
	IF old_param EQ 1 THEN vel

; Restore scale
	set_scale,old_minval,old_maxval,/quiet

	END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       OVERLAY_MYOPIC_DATA2
;

        PRO overlay_myopic_data2,_extra=e

        overlay_myopic_fov,_extra=e,/myopic_data

        END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       OVERLAY_MYOPIC_FOV2
;

        PRO overlay_myopic_fov2,beams_plot=beams_plot,gates_plot=gates_plot,colour=colour,thick=thick, $
		beams_sep=beams_sep,gates_sep=gates_sep,irev=irev,myopic_data=myopic_data, $
		no_time_set=no_time_set,altitude=altitude

        COMMON tjornes

        COMMON data
        COMMON prefs
        COMMON colour_info

; Variables
        radar_lat=63.77
        radar_lon=-20.54
        radar_boresite=30.0
        gate_length=15.0
        first_range=15.0
	IF NOT KEYWORD_SET(altitude) THEN altitude=110.0

; Half width of ATV image
        wid=siz_wat2/2

; Find radar cell locations
        x_wat=FLTARR(17,76)
        y_wat=FLTARR(17,76)
        FOR beam=0,16 DO BEGIN
                FOR gate=0,75 DO BEGIN
                        fov_pos=find_cell2(radar_lat,radar_lon,16,radar_boresite, $
                                3.24,gate_length,first_range,beam-0.5,gate)
                        lat=fov_pos(0)
                        lon=fov_pos(1)
                        zen=get_zenith_angle_from_tjornes(lat,lon,altitude)
                        IF NOT KEYWORD_SET(irev) THEN BEGIN
                                ;azm=get_azimuth_from_tjornes(lat,lon)+90.0+rog_wat2
                                azm=get_azimuth_from_tjornes(lat,lon)
                        ENDIF ELSE BEGIN
                                ;azm=-get_azimuth_from_tjornes(lat,lon)+90.0;+rog_wat2
                                azm=-get_azimuth_from_tjornes(lat,lon)
                        ENDELSE
                        IF azm GE 360.0 THEN azm=azm-360.0
                        tmp_pix=get_pixel(zen,azm)
                        x_wat(beam,gate)=tmp_pix(0)
                        y_wat(beam,gate)=tmp_pix(1)
                        ;x_wat(beam,gate)=wid+zen*COS(!PI*azm/180.0)*wid/75.0
                        ;y_wat(beam,gate)=wid+zen*SIN(!PI*azm/180.0)*wid/75.0
                ENDFOR
        ENDFOR

        IF NOT KEYWORD_SET(colour) THEN colour=254
        IF NOT KEYWORD_SET(thick) THEN thick=1

; If you want to plot the data as well ...
        IF KEYWORD_SET(myopic_data) THEN BEGIN

; Get data arrays
		IF NOT KEYWORD_SET(no_time_set) THEN BEGIN
                	go_time,FIX(STRMID(STRING(hms_wat2(now_wat2),FORMAT='(I6.6)'),0,4)), $
                        	FIX(STRMID(STRING(hms_wat2(now_wat2),FORMAT='(I6.6)'),4,2))
		ENDIF
                varr=get_scan()
                grnd=get_scan(/gnd_flag)

; Set color bar and levels, and switch colour map for velocity plot
                cin=FIX(FINDGEN(no_colours)*(ncol-4)/(no_colours-1))+1
                lvl=minval+FINDGEN(no_colours)*(maxval-minval)/no_colours
                IF parameter EQ 'vel' THEN cin=ROTATE(cin,2)

; Plot data
                FOR m=0,radar_beams-1 DO BEGIN
                        FOR n=11,75-1 DO BEGIN
                                IF varr(m,n) NE 10000 THEN BEGIN
                                        IF NOT((grnd(m,n) EQ 0 AND scatter EQ 1) OR (grnd(m,n) NE 0 AND scatter EQ 2)) THEN BEGIN
                                                ind=(MAX(WHERE(lvl LE varr(m,n))) > 0)
                                                IF parameter EQ 'vel' AND scatter EQ 3 AND grnd(m,n) EQ 1 THEN col=grey ELSE col=cin(ind)
                                                x_poly=[x_wat(m,n),x_wat(m+1,n),x_wat(m+1,n+1), $
                                                        x_wat(m,n+1),x_wat(m,n)]
                                                y_poly=[y_wat(m,n),y_wat(m+1,n),y_wat(m+1,n+1), $
                                                        y_wat(m,n+1),y_wat(m,n)]
                                                set_colour_table,/default
                                                POLYFILL,x_poly,y_poly,COL=col,NOCLIP=0
                                        ENDIF
                                ENDIF
                        ENDFOR
                ENDFOR
                set_scale_wat,min_wat2,max_wat2,/quiet
        ENDIF

; Plot the fov on ATV image
        IF KEYWORD_SET(beams_plot) AND KEYWORD_SET(gates_plot) THEN BEGIN
                FOR i=beams_plot,beams_plot DO BEGIN
                        FOR j=gates_plot(0),gates_plot(1) DO BEGIN
                                x_poly=[x_wat(i,j),x_wat(i+1,j),x_wat(i+1,j+1), $
                                        x_wat(i,j+1),x_wat(i,j)]
                                y_poly=[y_wat(i,j),y_wat(i+1,j),y_wat(i+1,j+1), $
                                        y_wat(i,j+1),y_wat(i,j)]
                                ;OPLOT,x_poly,y_poly,COL=254,THICK=thick*2
                                OPLOT,x_poly,y_poly,COL=colour,THICK=thick
                        ENDFOR
                ENDFOR
        ENDIF ELSE BEGIN

; Plot rough grids
                IF NOT KEYWORD_SET(normal) THEN BEGIN
                        FOR i=10,75,10 DO OPLOT,[x_wat(5:11,i)],[y_wat(5:11,i)],COL=colour,THICK=thick
                        ;FOR i=10,75,10 DO OPLOT,[x_wat(*,i)],[y_wat(*,i)],COL=colour,THICK=thick
                        FOR i=5,11,1 DO OPLOT,[x_wat(i,10:75)],[y_wat(i,10:75)],COL=colour,THICK=thick
                        ;FOR i=0,16, 8 DO OPLOT,[x_wat(i,10:75)],[y_wat(i,10:75)],COL=colour,THICK=thick
                ENDIF ELSE BEGIN
                        FOR i=0,75,15 DO PLOTS,[x_wat(*,i)],[y_wat(*,i)],COL=colour, $
                                THICK=thick,NORMAL=normal
                        FOR i=0,16, 8 DO PLOTS,[x_wat(i,*)],[y_wat(i,*)],COL=colour, $
                                THICK=thick,NORMAL=normal
                ENDELSE
        ENDELSE

        END
;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       OVERLAY_MYOPIC_ON_KEOGRAM2
;

        PRO overlay_myopic_on_keogram2

        COMMON tjornes

        COMMON data
        COMMON prefs
        COMMON colour_info

; Variables
        radar_lat=63.77
        radar_lon=-20.54
        radar_boresite=30.0
        gate_length=15.0
        first_range=15.0
        altitude=110.0
        lat_tjo=66.20

; Find radar cell locations
        zan_cell=FLTARR(16,76)
        FOR beam=0,15 DO BEGIN
                FOR gate=11,75 DO BEGIN
                        fov_pos=find_cell2(radar_lat,radar_lon,16,radar_boresite, $
                                3.24,gate_length,first_range,beam,gate)
                        lat=fov_pos(0)
                        lon=fov_pos(1)
                        zan_cell(beam,gate)=get_zenith_angle_from_tjornes(lat,lon,altitude)
                        IF lat LT lat_tjo THEN zan_cell(beam,gate)=-zan_cell(beam,gate)
                ENDFOR
        ENDFOR

; Find beam data blocks which have the correct look direction and frequency
        beamno=beam_no
        plot_beams=WHERE(beam_dir EQ beamno AND beam_time GE rti_start-120 AND beam_time LE rti_end+120, $
                no_plot_beams)

; Set color bar and levels, and switch colour map for velocity plot
        cin=FIX(FINDGEN(no_colours)*(ncol-4)/(no_colours-1))+1
        lvl=minval+FINDGEN(no_colours)*(maxval-minval)/no_colours
        IF parameter EQ 'vel' THEN cin=ROTATE(cin,2)

; Cycle through beams to plot
        day_start=FIX(rti_start/86400L)*86400L
        FOR m=0,no_plot_beams-2 DO BEGIN
                FOR n=11,74 DO BEGIN

                        IF a(n,plot_beams(m)) NE 10000 THEN BEGIN

                                start_time=beam_time(plot_beams(m))-day_start
                                end_time=beam_time(plot_beams(m+1))-day_start

                                IF NOT((gscat(n,plot_beams(m)) EQ 0 AND scatter EQ 1) OR $
                                        (gscat(n,plot_beams(m)) NE 0 AND scatter EQ 2)) THEN BEGIN
                                
                                                ind=(MAX(WHERE(lvl LE ((a(n,plot_beams(m)) > minval) < maxval))) > 0)
                                                IF parameter EQ 'vel' AND scatter EQ 3 AND gscat(n,plot_beams(m)) EQ 1 THEN col=grey ELSE col=cin(ind)
                                                POLYFILL,[start_time/3600.0,start_time/3600.0,end_time/3600.0,end_time/3600.0], $
                                                         [zan_cell(beamno,n),zan_cell(beamno,n+1),zan_cell(beamno,n+1),zan_cell(beamno,n)],COL=col,NOCLIP=0

                                ENDIF
                        ENDIF
                ENDFOR
        ENDFOR

        END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       OVERLAY_ZENITH2
;

        PRO overlay_zenith2

        COMMON tjornes

        COMMON colour_info
        update_colour_info

; Zenith angle limit
	zenith_limit=75.0

; Plot circles (from 10 to 70, div is 10 deg)
        xx=FLTARR(361) & yy=FLTARR(361)
        c_img=[siz_wat2/2,siz_wat2/2+0.5]
        FOR j=15,75,15 DO BEGIN
                d_img=FLOAT(j)*(siz_wat2/2)/zenith_limit
                FOR i=0,360 DO BEGIN
                        xx(i)=c_img(0)+1*d_img*COS(i*!PI/180.0)
                        yy(i)=c_img(1)+1*d_img*SIN(i*!PI/180.0)
                ENDFOR
                set_colour_table,/default
                OPLOT,xx,yy,COL=white,NOCLIP=0,LINE=1
        ENDFOR

	END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       OVERLAY_ZENITH_SPECIAL2
;

        PRO overlay_zenith_special2

        COMMON tjornes

        COMMON colour_info
        update_colour_info

; Zenith angle limit
	zenith_limit=75.0

; Plot circles (from 10 to 70, div is 10 deg)
        xx=FLTARR(361) & yy=FLTARR(361)
        c_img=[siz_wat2/2,siz_wat2/2+0.5]
	ele_ang=[30,60,75]
        FOR k=0,2 DO BEGIN
		j=ele_ang(k)
                d_img=FLOAT(j)*(siz_wat2/2)/zenith_limit
                FOR i=0,360 DO BEGIN
                        xx(i)=c_img(0)+1*d_img*COS(i*!PI/180.0)
                        yy(i)=c_img(1)+1*d_img*SIN(i*!PI/180.0)
                ENDFOR
                set_colour_table,/default
                OPLOT,xx,yy,COL=254,NOCLIP=0,LINE=0,THICK=6
                OPLOT,xx,yy,COL=1,NOCLIP=0,LINE=0,THICK=4
        ENDFOR

	END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;       MASK_FOV2
;

        PRO mask_fov2

        COMMON tjornes

        xx=FLTARR(361) & yy=FLTARR(361) & xe=FLTARR(361) & ye=FLTARR(361)
        c_img=[siz_wat2/2,siz_wat2/2]
        d_img=FIX(a_v_wat2*!PI/2.0)
        FOR i=0,360 DO BEGIN
                xx(i)=c_img(0)+1*d_img*COS(i*!PI/180.0)
                yy(i)=c_img(1)+1*d_img*SIN(i*!PI/180.0)
                xe(i)=c_img(0)+2*d_img*COS(i*!PI/180.0)
                ye(i)=c_img(1)+2*d_img*SIN(i*!PI/180.0)
        ENDFOR
        set_colour_table,/default
        FOR i=0,359 DO BEGIN
                POLYFILL,[xx(i),xe(i),xe(i+1),xx(i+1)],[yy(i),ye(i),ye(i+1),yy(i+1)], $
                        COL=253,NOCLIP=0
        ENDFOR

	END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	DEFINE_ZENITH_TJORNES2
;

	PRO define_zenith_tjornes2,preview=preview

; Image used for defining zenith
	fname='./wat_tjo_11-09-30_22-00-04-35.jpg'

; Time info for input to zenith2000
	date='20110930'
	hh='22'
	mmss='0004'

; Read image and identify size
	READ_JPEG,fname,image
	image_size=SIZE(image,/DIM)
	xsize=image_size(0)
	ysize=image_size(1)
	PRINT,xsize,ysize
	WINDOW,XSIZE=xsize,YSIZE=ysize

; Plot image
	LOADCT,0
	LOADCT,39
	TVSCL,image
	!P.POSITION=[0,0,1,1]
	PLOT,[0],[0],XRANGE=[0,xsize],YRANGE=[0,ysize],/XSTYLE,/YSTYLE

	IF KEYWORD_SET(preview) THEN RETURN

	psym8

; Get position of Capera (StarID=4)
	PRINT,'Click Capera (in Auriga)'
	CURSOR,x_tmp,y_tmp,/DEVICE
	SPAWN,'sleep 1'
	P=MAX(image(x_tmp-5:x_tmp+5,y_tmp-5:y_tmp+5),i)
	ind=ARRAY_INDICES(image(x_tmp-5:x_tmp+5,y_tmp-5:y_tmp+5),i)
	x1=x_tmp-5+ind(0)
	y1=y_tmp-5+ind(1)
	x1=x_tmp
	y1=y_tmp
	OPLOT,[x1],[y1],PSYM=8,SYMSIZE=1,COL=254
	PRINT,'Capera  =[',x1,y1,']'

; Get position of Vega (StarID=11)
	PRINT,'Click Vega (in Lyra)'
	CURSOR,x_tmp,y_tmp,/DEVICE
	SPAWN,'sleep 1'
	P=MAX(image(x_tmp-5:x_tmp+5,y_tmp-5:y_tmp+5),i)
	ind=ARRAY_INDICES(image(x_tmp-5:x_tmp+5,y_tmp-5:y_tmp+5),i)
	x2=x_tmp-5+ind(0)
	y2=y_tmp-5+ind(1)
	x2=x_tmp
	y2=y_tmp
	OPLOT,[x2],[y2],PSYM=8,SYMSIZE=1,COL=254
	PRINT,'Vega    =[',x2,y2,']'

; Get position of Polaris (StarID=1)
	PRINT,'Click Polaris (in Ursa Minor)'
	CURSOR,x_tmp,y_tmp,/DEVICE
	SPAWN,'sleep 1'
	P=MAX(image(x_tmp-5:x_tmp+5,y_tmp-5:y_tmp+5),i)
	ind=ARRAY_INDICES(image(x_tmp-5:x_tmp+5,y_tmp-5:y_tmp+5),i)
	x3=x_tmp-5+ind(0)
	y3=y_tmp-5+ind(1)
	x3=x_tmp
	y3=y_tmp
	OPLOT,[x3],[y3],PSYM=8,SYMSIZE=1,COL=254
	PRINT,'Polaris    =[',x3,y3,']'

; Get position of Deneb (StarID=13)
	PRINT,'Click Deneb (in Cygnus)'
	CURSOR,x_tmp,y_tmp,/DEVICE
	SPAWN,'sleep 1'
	P=MAX(image(x_tmp-5:x_tmp+5,y_tmp-5:y_tmp+5),i)
	ind=ARRAY_INDICES(image(x_tmp-5:x_tmp+5,y_tmp-5:y_tmp+5),i)
	x4=x_tmp-5+ind(0)
	y4=y_tmp-5+ind(1)
	x4=x_tmp
	y4=y_tmp
	OPLOT,[x4],[y4],PSYM=8,SYMSIZE=1,COL=254
	PRINT,'Deneb    =[',x4,y4,']'

; Center of the image (identified with CURSOR,x0,y0,/DEVICE)
	x_cnt=xsize/2
	y_cnt=ysize/2

; Try Capera + Vega
	PRINT,' '
	PRINT,'Capera + Vega'
	params=date+hh+'.'+mmss+' '+ $
		STRCOMPRESS(STRING(x_cnt))+' '+STRCOMPRESS(STRING(y_cnt))+' '+ $
		'4 '+STRCOMPRESS(STRING(x1))+' '+STRCOMPRESS(STRING(y1))+' '+ $
		'11 '+STRCOMPRESS(STRING(x2))+' '+STRCOMPRESS(STRING(y2))
	SPAWN,'/radar01/work/Aurora/Starmap/zenith2000_tjornes.sh '+params

; Try Capera + Deneb
	PRINT,' '
	PRINT,'Capera + Deneb'
	params=date+hh+'.'+mmss+' '+ $
		STRCOMPRESS(STRING(x_cnt))+' '+STRCOMPRESS(STRING(y_cnt))+' '+ $
		'4 '+STRCOMPRESS(STRING(x1))+' '+STRCOMPRESS(STRING(y1))+' '+ $
		'13 '+STRCOMPRESS(STRING(x4))+' '+STRCOMPRESS(STRING(y4))
	SPAWN,'/radar01/work/Aurora/Starmap/zenith2000_tjornes.sh '+params

; Try Polaris + Deneb
	PRINT,' '
	PRINT,'Polaris + Deneb'
	params=date+hh+'.'+mmss+' '+ $
		STRCOMPRESS(STRING(x_cnt))+' '+STRCOMPRESS(STRING(y_cnt))+' '+ $
		'1 '+STRCOMPRESS(STRING(x3))+' '+STRCOMPRESS(STRING(y3))+' '+ $
		'13 '+STRCOMPRESS(STRING(x4))+' '+STRCOMPRESS(STRING(y4))
	SPAWN,'/radar01/work/Aurora/Starmap/zenith2000_tjornes.sh '+params

	END

;------------------------------------------------------------------------------------------------------------------
; NAME:
;
;	DEFINE_ZENITH_SYOWA
;

	PRO define_zenith_syowa2,preview=preview,confirm=confirm

; Image used for defining zenith
	fname='./SYO_Watec_11-09-30_22-00-08-27.jpg'

; Time info for input to zenith2000
	date='20110930'
	hh='22'
	mmss='0008'

; Read image and identify size
	READ_JPEG,fname,image
	image_size=SIZE(image,/DIM)
	xsize=image_size(0)
	ysize=image_size(1)
	PRINT,xsize,ysize
	WINDOW,XSIZE=xsize,YSIZE=ysize

; Plot image
	LOADCT,0
	;LOADCT,39
	;cut_col_tab
	image=BYTSCL(image,MAX=50.0,MIN=0.0,TOP=252.0)
	TVSCL,image
	!P.POSITION=[0,0,1,1]
	PLOT,[0],[0],XRANGE=[0,xsize],YRANGE=[0,ysize],/XSTYLE,/YSTYLE
	psym8

	IF KEYWORD_SET(preview) THEN RETURN

	IF KEYWORD_SET(confirm) THEN BEGIN

		X0=318
		Y0=248
		A=123.0
		ang=-106.0

		OPLOT,[X0],[Y0],PSYM=8,COL=254,SYMSIZE=2
		OPLOT,X0+[-500,500],[Y0,Y0],LINE=1
		OPLOT,[X0,X0],Y0+[-500,500],LINE=1
		dim=A*!PI*0.5
		x_cir=FLTARR(361)
		y_cir=FLTARR(361)
		FOR i=0,360 DO BEGIN
			x_cir(i)=X0+dim*COS(i*!DTOR)
			y_cir(i)=Y0+dim*SIN(i*!DTOR)
		ENDFOR	
		OPLOT,x_cir,y_cir

		RETURN
	ENDIF


; Get position of Sirius (StarID=6)
	PRINT,'Click Sirius'
	CURSOR,x_tmp,y_tmp,/DEVICE
	SPAWN,'sleep 1'
	P=MAX(image(x_tmp-5:x_tmp+5,y_tmp-5:y_tmp+5),i)
	ind=ARRAY_INDICES(image(x_tmp-5:x_tmp+5,y_tmp-5:y_tmp+5),i)
	x1=x_tmp-5+ind(0)
	y1=y_tmp-5+ind(1)
	x1=x_tmp
	y1=y_tmp
	OPLOT,[x1],[y1],PSYM=8,SYMSIZE=1,COL=254
	PRINT,'Sirius  =[',x1,y1,']'

; Get position of Fomalhaut (StarID=14)
	PRINT,'Click Fomalhaut'
	CURSOR,x_tmp,y_tmp,/DEVICE
	SPAWN,'sleep 1'
	P=MAX(image(x_tmp-5:x_tmp+5,y_tmp-5:y_tmp+5),i)
	ind=ARRAY_INDICES(image(x_tmp-5:x_tmp+5,y_tmp-5:y_tmp+5),i)
	x2=x_tmp-5+ind(0)
	y2=y_tmp-5+ind(1)
	x2=x_tmp
	y2=y_tmp
	OPLOT,[x2],[y2],PSYM=8,SYMSIZE=1,COL=254
	PRINT,'Fomalhaut    =[',x2,y2,']'

; Get position of Canopus (StarID=16)
	PRINT,'Click Canopus'
	CURSOR,x_tmp,y_tmp,/DEVICE
	SPAWN,'sleep 1'
	P=MAX(image(x_tmp-5:x_tmp+5,y_tmp-5:y_tmp+5),i)
	ind=ARRAY_INDICES(image(x_tmp-5:x_tmp+5,y_tmp-5:y_tmp+5),i)
	x3=x_tmp-5+ind(0)
	y3=y_tmp-5+ind(1)
	x3=x_tmp
	y3=y_tmp
	OPLOT,[x3],[y3],PSYM=8,SYMSIZE=1,COL=254
	PRINT,'Canopus    =[',x3,y3,']'

; Get position of Acrux (StarID=17)
	PRINT,'Click Acrux' 
	CURSOR,x_tmp,y_tmp,/DEVICE
	SPAWN,'sleep 1'
	P=MAX(image(x_tmp-5:x_tmp+5,y_tmp-5:y_tmp+5),i)
	ind=ARRAY_INDICES(image(x_tmp-5:x_tmp+5,y_tmp-5:y_tmp+5),i)
	x4=x_tmp-5+ind(0)
	y4=y_tmp-5+ind(1)
	x4=x_tmp
	y4=y_tmp
	OPLOT,[x4],[y4],PSYM=8,SYMSIZE=1,COL=254
	PRINT,'Acrux    =[',x4,y4,']'

; Center of the image (identified with CURSOR,x0,y0,/DEVICE)
	x_cnt=xsize/2
	y_cnt=ysize/2

; Try Sirius + Fomalhaut
	PRINT,' '
	PRINT,'Sirius + Fomalhaut'
	params=date+hh+'.'+mmss+' '+ $
		STRCOMPRESS(STRING(x_cnt))+' '+STRCOMPRESS(STRING(y_cnt))+' '+ $
		'6 '+STRCOMPRESS(STRING(x1))+' '+STRCOMPRESS(STRING(y1))+' '+ $
		'14 '+STRCOMPRESS(STRING(x2))+' '+STRCOMPRESS(STRING(y2))
	SPAWN,'/radar01/work/Aurora/Starmap/zenith2000_syowa.sh '+params

; Try Sirius + Acrux
	PRINT,' '
	PRINT,'Sirius + Acrux'
	params=date+hh+'.'+mmss+' '+ $
		STRCOMPRESS(STRING(x_cnt))+' '+STRCOMPRESS(STRING(y_cnt))+' '+ $
		'6 '+STRCOMPRESS(STRING(x1))+' '+STRCOMPRESS(STRING(y1))+' '+ $
		'17 '+STRCOMPRESS(STRING(x4))+' '+STRCOMPRESS(STRING(y4))
	SPAWN,'/radar01/work/Aurora/Starmap/zenith2000_syowa.sh '+params

; Try Canopus + Fomalhaut
	PRINT,' '
	PRINT,'Canopus + Fomalhaut'
	params=date+hh+'.'+mmss+' '+ $
		STRCOMPRESS(STRING(x_cnt))+' '+STRCOMPRESS(STRING(y_cnt))+' '+ $
		'16 '+STRCOMPRESS(STRING(x3))+' '+STRCOMPRESS(STRING(y3))+' '+ $
		'14 '+STRCOMPRESS(STRING(x2))+' '+STRCOMPRESS(STRING(y2))
	SPAWN,'/radar01/work/Aurora/Starmap/zenith2000_syowa.sh '+params

; Try Acrux + Fomalhaut
	PRINT,' '
	PRINT,'Acrux + Fomalhaut'
	params=date+hh+'.'+mmss+' '+ $
		STRCOMPRESS(STRING(x_cnt))+' '+STRCOMPRESS(STRING(y_cnt))+' '+ $
		'14 '+STRCOMPRESS(STRING(x2))+' '+STRCOMPRESS(STRING(y2))+' '+ $
		'17 '+STRCOMPRESS(STRING(x4))+' '+STRCOMPRESS(STRING(y4))
	SPAWN,'/radar01/work/Aurora/Starmap/zenith2000_syowa.sh '+params

	END
