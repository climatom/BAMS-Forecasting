#!/bin/sh
HOST='climate.appstate.edu'
USER='appalair'
PASSWD='app1.air'
DIR='/home/appalair/Everest_Forecasting-master/Data/ECMWF/ecaccess/Logs/'
FILE='*.txt'
WORKDIR='/home/lunet/gytm3/Everest2019/Forecasting/ECMWF/Logging'
cd $WORKDIR

ftp -n $HOST <<END_SCRIPT
quote USER $USER
quote PASS $PASSWD
binary
passive
prompt
cd $DIR
mget $FILE
quit
END_SCRIPT
exit 0
