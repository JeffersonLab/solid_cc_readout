#!/bin/bash

###########
#  USER   #
###########

source /home/rich/LaserTED/swFE/bin/RUNME.cfg



DISABLEFULLGRID=0 # skip [gain,thr] = [x2,25][x4,25][x4,50] use it for high gain mapmt



############
#  EXPERT  #
############
debug=0 # use 0 for data taking
ENABLECALG=0 # enable online gain calibration
ADDR2="192.168.2.10" # TOP
ADDR1="192.168.1.10" # BOTTOM

ADCNBIT=12 # 12, 10 or  8


##############
#  FUNCTION  #
##############

ADC_RESOLUTION(){

  BITS=$1
  TEN=0
  EIGHT=0

  if [ $BITS -eq 8 ]; then
    EIGHT=1
  fi

  if [ $BITS -eq 10 ]; then
    TEN=1
  fi
}

SET_FPGA(){
  ip=$1 
}


SET_MAROC_ID(){
  mrc0=$1
  mrc1=$2
  mrc2=$3
  printf -v TILEFOLDER '../tile/%s_%s_%s/' $mrc0 $mrc1 $mrc2
  printf -v TILEFOLDERBIS '../tile/%s_%s_%s_' $mrc0 $mrc1 $mrc2
  printf -v TILEFOLDERTRIS '../tile/%s_%s_%s' $mrc0 $mrc1 $mrc2
  printf -v TILENOPATH '%s_%s_%s' $mrc0 $mrc1 $mrc2
  printf -v DB '../db/%s_%s_%s' $mrc0 $mrc1 $mrc2
  printf -v OUTFOLDER ../out/

  rm -rf $TILEFOLDER 
  mkdir $TILEFOLDER

}

SET_MAPMT(){
  mapmt0=$1
  mapmt1=$2
  mapmt2=$3
}

SET_LASER_POSITION(){
  x=$1
  y=$2
}


SET_TDC(){
  SOURCE=1; # 0 CTEST, 1 ANODES
}


SET_TDC_EVBUILDER(){
  TRIG_DELAY=$1 
  LOOKBACK=$2
  WINDOW=$3 
}


SET_PULSER(){
#  PULSE_FREQ=10000 # TDC +ADC
#  PULSE_FREQ=1500 # TDC +ADC
  PULSE_FREQ=1000 # TDC +ADC
  #PULSE_FREQ=500000 # TDC only
  PULSE_DUTYCYCLE=0.9
  PULSE_REPETITION=100000000
}


SET_STATISTICS(){
  TIME=360000  # time in s
  EVENTS=1000000   # 10^6
 # EVENTS=100000   # 10^5
 # EVENTS=10000   # 10^4
 #  EVENTS=1000 # 10^3
}


SET_GAIN_MODE(){
  GAIN_MODE=$1 # 0 FOR NO EQUALIZAITON, 1 FOR MAPS
}


SET_GAIN_DEFAULT(){
  GAIN=$1 
}

SET_THRESHOLD(){
  # absolute threshold 
  #THR_MODE=0
  #THR=250  

  # relative to the pedestal file in ../db/<tile>/chip0.txt
  THR=$1
  THR_MODE=1
}

SET_ADC(){
 ENABLE_ADC=1
}



SET_CTEST(){
  # CHANNEL TEST MASK (CTEST PIN, AUXILIARY INPUT)
  # 0 = disabled 
  # 1 = single channel (ch_sel), 
  # 2 = all channels
  # 3 = all but one (ch_sel)
  CTEST_MODE=0 
  CTEST_CHSEL=23 
}



SET_CHARGE(){
  CTEST_AMPLITUDE=$1 
}

SET_PROBING(){

  # 0..63 to enable an individual channel (oscilloscope or for DC measurement)
  # -1 to disable  
  PROBE=$1
}

SET_RUN_ID(){
 
   RUNID=$1
}

###########
# STORAGE #
###########


STORE_TILE(){

  NOW=$(date +"%Y_%m_%d_%H_%M")
  echo $NOW

  if [ $debug -eq 1 ]; then 
    mv $TILEFOLDER $TILEFOLDERBIS$NOW_$debug
    echo Data stored in $TILEFOLDERBIS$NOW_$debug
  else
    mv $TILEFOLDER $TILEFOLDERBIS$NOW
    echo Data stored in $TILEFOLDERBIS$NOW
  fi
  cp RUNME.cfg $TILEFOLDERBIS$NOW

  #echo ./PARSING $TILEFOLDERBIS$NOW 1 60 1
  #./PARSING.sh $TILEFOLDERBIS$NOW 1 60 1 &
}


############
# COMMANDS #
############

DO_WHEEL(){
   w=$1 
   wheel $w
}

DO_LASER_POSITION(){
  lts_run $x $y
}

DO_HV(){

  hvon $hv
  echo HV RAMPING, WAIT $RAMP seconds
  sleep $RAMP
  echo HV SETTLING, WAIT $WARMUP seconds
  sleep $WARMUP  
  echo DO_HV $hv DONE

}

DO_SLOWCONTROL(){
  ./daq run.address=$ip run.daq_mode=2 \
  run.maroc.id.[0]=$mrc0  run.maroc.id.[1]=$mrc1  run.maroc.id.[2]=$mrc2 
}

DO_RUN(){

   PREFIX="run"
   DAQMODE=1 # events

   if [ $debug -eq 1 ]; then 
     echo W$w TILE $TILE  RUN $totRun HV $hv G $GG THR $thr 
     EVENTS=1000

  else

    laser 0
    echo  sleep 1
    sleep 1

    laser 1
    echo  sleep 1
    sleep 1


    ./daq run.name=$PREFIX run.address=$ip run.id=$RUNID run.daq_mode=$DAQMODE \
    run.maroc.gain_mode=$GAIN_MODE run.maroc.gain_default=$GAIN \
    run.maroc.thr_default=$THR run.maroc.thr_mode=$THR_MODE \
    run.tdc.trigger_delay=$TRIG_DELAY \
    run.tdc.evtb_lookback=$LOOKBACK run.tdc.evtb_windowwidth=$WINDOW \
    run.maroc.EN_ADC=$ENABLE_ADC run.adc.enable_adc=$ENABLE_ADC run.adc.hold_delay=13 \
     run.maroc.ramp_8bit=$EIGHT run.maroc.ramp_10bit=$TEN \
    run.event_preset=$EVENTS run.time_preset=$TIME \
    run.ctest_amplitude=$CTEST_AMPLITUDE run.maroc.ctest_mode=$CTEST_MODE run.maroc.ch_sel=$CTEST_CHSEL \
    run.maroc.ch_probe=$PROBE \
    run.source_type=$SOURCE \
    run.pulser.frequency=$PULSE_FREQ run.pulser.dutycycle=$PULSE_DUTYCYCLE run.pulser.repetition=$PULSE_REPETITION \
    run.maroc.id.[0]=$mrc0  run.maroc.id.[1]=$mrc1  run.maroc.id.[2]=$mrc2 \
    run.mapmt.id.[0]=$mapmt0 run.mapmt.id.[1]=$mapmt1 run.mapmt.id.[2]=$mapmt2 \
    run.mapmt.hv=$hv \
    run.laser.y=$y run.laser.x=$x run.laser.w=$w 

    laser 0
  fi

  DO_LOG

  totRun=$((totRun+1))
}

SET_SCALER(){

  PREFIX="scaler"
  DAQMODE=0 # scaler
  DURATION=10000 # Scaler window in milliseconds
  REPETITION=1 # Number of iterations
}

SET_SCALER_DARK(){

  PREFIX="scaler_dark"
  PULSE_REPETITION=0
  DAQMODE=0 # scaler
  DURATION=10000 # Scaler window in milliseconds
  REPETITION=1 # Number of iterations
}


DO_SCALER(){
   
  laser $1;
  time ./daq run.name=$PREFIX run.address=$ip run.id=$RUNID run.daq_mode=$DAQMODE \
  run.slowcontrol.repetition=$REPETITION run.slowcontrol.time_interval=$DURATION \
  run.maroc.gain_mode=$GAIN_MODE run.maroc.gain_default=$GAIN \
  run.maroc.thr_default=$THR run.maroc.thr_mode=$THR_MODE \
  run.ctest_amplitude=$CTEST_AMPLITUDE run.maroc.ctest_mode=$CTEST_MODE run.maroc.ch_sel=$CTEST_CHSEL \
  run.maroc.ch_probe=$PROBE \
  run.source_type=$SOURCE \
  run.pulser.frequency=$PULSE_FREQ run.pulser.dutycycle=$PULSE_DUTYCYCLE \
  run.pulser.repetition=$PULSE_REPETITION \
  run.maroc.id.[0]=$mrc0  run.maroc.id.[1]=$mrc1  run.maroc.id.[2]=$mrc2 \
  run.mapmt.id.[0]=$mapmt0 run.mapmt.id.[1]=$mapmt1 run.mapmt.id.[2]=$mapmt2 \
  run.mapmt.hv=$hv \
  run.laser.y=$y run.laser.x=$x run.laser.w=$w 
  DO_LOG
}


DO_LOG(){
 printf -v ENTRY  '%3d W %d X %4d Y %4d HV %5d GMODE %2d GAIN %5s THR %5d %s\n'  $totRun $w $x $y $hv $GAIN_MODE $GAIN $THR $g$PREFIX
 echo $ENTRY >> ./logbook.txt
 #echo $ENTRY appended to ./logbook.txt
}



SET_HV(){

  hv=$1

  RAMP=60
  WARMUP=60
  if [ $debug -eq 1 ]; then 
    RAMP=0 
    WARMUP=0
  fi

  DO_HV
  DO_SLOWCONTROL
}



DO_SLEEP(){

  SECONDS=$1 
  echo sleep $SECONDS
  sleep $SECONDS
}


DO_LVOFF(){
  lvoff
  DO_SLEEP 1
}

DO_LVON(){
  lvon
  DO_SLEEP 3
}




INIT_SW(){  
 
  rm -f  ./fpgaMonitor.txt
  rm -f ./logbook.txt

  if [ "$TILE" = "TOP" ]; then 
    SET_LASER_POSITION 200 25 
    SET_FPGA $ADDR2
    SET_MAROC_ID $TOP_MAROC0 $TOP_MAROC1 $TOP_MAROC2
    SET_MAPMT $TOP_MAPMT0 $TOP_MAPMT1 $TOP_MAPMT2
  else
    SET_LASER_POSITION 200 125 
    SET_FPGA $ADDR1
    SET_MAROC_ID $BOTTOM_MAROC0 $BOTTOM_MAROC1 $BOTTOM_MAROC2
    SET_MAPMT $BOTTOM_MAPMT0 $BOTTOM_MAPMT1 $BOTTOM_MAPMT2
  fi

  SET_ADC 1
  SET_TDC_EVBUILDER 30 30 30 
  SET_PROBING -1
  SET_TDC
  SET_STATISTICS
  SET_GAIN_MODE 0
  SET_GAIN_DEFAULT 64
  SET_THRESHOLD 50
  SET_PULSER
  SET_CTEST
  #SET_CHARGE 400 # [0..4095] that corresponds roughly to [5fC..2.5pC]

  w=0 

  echo ADC_RESOLUTION $ADCNBIT
  ADC_RESOLUTION $ADCNBIT    
  echo sleep 1 
  sleep 1

}

DO_PEDESTAL(){
    laser 0
    hvoff 0


  if [ $debug -eq 1 ]; then 
    echo " DEBUG MODE -------------------> PEDESTAL"
  else

  for GG in 32 64 128 #16 32 64 128 255
  do
    ./PED.sh $mrc0 $mrc1 $mrc2 $GG 0 $ip  
    printf -v GPEDFILEFORANA '%s/pedestal0_%03d.txt' $TILEFOLDER  $GG
    printf -v GPEDFILEFORDAQ '%s/chip0_%03d.txt' $TILEFOLDER  $GG
    echo cp $DB/pedestal0.txt  $GPEDFILEFORANA
    cp $DB/pedestal0.txt  $GPEDFILEFORANA
    echo cp $DB/chip0.txt  $GPEDFILEFORDAQ
    cp $DB/chip0.txt  $GPEDFILEFORDAQ
  done 
  fi   
  rm -f $TILEFOLDER/*.log
  rm -f $TILEFOLDER/*.bin 
  rm -f $TILEFOLDER/*.pdf
  rm -f $TILEFOLDER/*.root 
  rm -f $TILEFOLDER/pedestal0.txt
  #rm -f $TILEFOLDER/chip0.txt
  rm -f $TILEFOLDER/ped*ADC* 

    DO_SLOWCONTROL


}



DO_DUMMY(){

  # DUMMY RUN to skip the problem of no signal in the first run...
  TIME=10
  SET_GAIN_MODE 0
  SET_GAIN_DEFAULT 64
  SET_THRESHOLD 50
  SET_HV 1050
  SET_RUN_ID 1050
  DO_RUN
  rm $TILEFOLDER/*1050*
  sed -i 1d ./logbook.txt  
}

DARK(){

  SET_SCALER_DARK

  rm -f $OUTFOLDER$PREFIX*

  for hv in 1000 1010 1020 1030 1040 1050 1060 1070 1080 1090 1100 #1025 1050 1075 1100  
  do
    SET_HV $hv
    for GG in 64
    do
      SET_GAIN_MODE 0
      SET_GAIN_DEFAULT $GG
      for threshold in 10 25 50 75 100 200  
      do
        SET_THRESHOLD $threshold
        DO_SLOWCONTROL
        SET_RUN_ID $hv
        totRun=$hv
        DO_SCALER 0
        printf -v DARKFILE_OLD '%sscaler_dark_%06d.txt' $OUTFOLDER  $hv
        printf -v DARKLOG_OLD '%sscaler_dark_%06d.log' $OUTFOLDER  $hv
        printf -v DARKFILE_NEW '%sscaler_dark_HV%04d_GAIN%04d_THR%04d.txt' $OUTFOLDER  $hv $GAIN $THR
        printf -v DARKLOG_NEW '%sscaler_dark_HV%04d_GAIN%04d_THR%04d.log' $OUTFOLDER  $hv $GAIN $THR
        echo mv $DARKFILE_OLD $DARKFILE_NEW
        mv $DARKFILE_OLD $DARKFILE_NEW
        mv $DARKLOG_OLD $DARKLOG_NEW
        DO_SLOWCONTROL
      done # loop on THRESHOLD
    done # loop on GAIN
  done # loop on HV

  echo  ls -ltr $OUTFOLDER$PREFIX*
  ls -ltr $OUTFOLDER$PREFIX*
  echo  ls -ltr $TILEFOLDER
  ls -ltr $TILEFOLDER


  DARKFOLDER="SCALERDARK"
  echo mkdir $DARKFOLDER
  mkdir $DARKFOLDER
  echo mv $OUTFOLDER$PREFIX* $DARKFOLDER
  mv $OUTFOLDER$PREFIX* $DARKFOLDER
  echo mv $DARKFOLDER $TILEFOLDER
  mv $DARKFOLDER $TILEFOLDER
#  rm -r $DARKFOLDER
}


COPY_LOGBOOK(){
  cat ./logbook.txt
  mv ./logbook.txt $TILEFOLDER
}

COPY_FPGAMONITOR(){
  cat ./fpgaMonitor.txt
  mv ./fpgaMonitor.txt $TILEFOLDER
}



DO_DATA_ACQUISITION(){

 echo "-----------------------------------------------> RUN "$totRun
  #totRun=1
  TIME=36000 # max run time in s
#  SET_GAIN_DEFAULT 64 
#  SET_THRESHOLD 50
  
  for hv in 1000
  do
   # SET_HV $hv   # turn HV on
    for thr in 600
    do
      SET_THRESHOLD $thr
#      for GG in 128 # use gain 2x
#      for GG in 64 # use gain 1x
#      for GG in 48 # use gain 0.75x
#      for GG in 32  # use gain 0.5x
      for GG in 24  # use gain 0.375x
#      for GG in 16  # use gain 0.25x
#      for GG in 8  # use gain 0.125x
#      for GG in 4  # use gain 0.0625x
      do
       
	sleep 40
        SET_GAIN_MODE 0  #use gain default (8 bit MAROC amp from 0 to 4x, 255 for 4x, 64 for 1x)
#        SET_GAIN_MODE 1  #use gain map from dir cfg/gain, orginial from gain_test
        SET_GAIN_DEFAULT $GG
        SET_RUN_ID $totRun
        DO_RUN

        prevRun=$((totRun-1))
        printf -v FILESEED '%srun_%06d' $TILEFOLDER $prevRun
        printf -v PEDFILE '%s.ped' $FILESEED   
        printf -v GPEDFILEFORANA '%s/pedestal0_%03d.txt' $TILEFOLDER  $GG
        #cp $GPEDFILEFORANA $PEDFILE
        cp $GPEDFILEFORANA $PEDFILE 2>/dev/null # to suppress error output in bash
        
        if [ $ENABLECALG -eq 1 ]; then  
          ### ONLINE CALIB
          #eqmode=BES
           eqmode=TOT
          # eqmode=PMT
          printf -v LOGFILE '%s.log' $FILESEED
          ./reco $LOGFILE  

          ./ana  $LOGFILE $PEDFILE
          printf -v PARSTRING 'W%s_HV%s_THR%s' $w $hv $thr
          printf -v cmd 'plotADC.C++("%s","%s",192)' $FILESEED $PARSTRING    
          root -l -b -q  $cmd
          ./MapsGenerator_clean.sh $hv $thr $prevRun $TILEFOLDER $eqmode $mapmt0 $mapmt1 $mapmt2
          SET_GAIN_MODE 1
          SET_RUN_ID $totRun        
          DO_RUN
          prevRun=$((totRun-1))
          printf -v FILESEED '%srun_%06d' $TILEFOLDER $prevRun
          printf -v PEDFILE '%s.ped' $FILESEED   
          echo cp $TILEFOLDER/pedestal0.txt $PEDFILE
          cp $TILEFOLDER/pedestal0.txt $PEDFILE
        fi
      done
    done
  done 
}



LASER(){

  w=$1
  EVENTS=$2
  DO_WHEEL $w 
#  DO_LASER_POSITION 
  DO_SLOWCONTROL
  DO_DATA_ACQUISITION
  DO_SLOWCONTROL

}

#########
# MAIN  #
#########

for TILE in "TOP" "BOTTOM"
#for TILE in "TOP"
#for TILE in "BOTTOM"
do
  INIT_SW    
#  DO_LVOFF
#  DO_LVON
#  DO_PEDESTAL
#  DO_DUMMY
#  echo sleep 30s and start code go!
#  sleep 30
  totRun=1 # Run Counter
	 LASER 3 100000
#  LASER 1 10000000
#  LASER 2 10000000
#  LASER 3 10000000
#  DARK
  COPY_LOGBOOK
  COPY_FPGAMONITOR
  STORE_TILE
done
#eof
