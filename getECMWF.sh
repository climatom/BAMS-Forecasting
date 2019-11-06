#!/bin/bash

# Simple script to grab ECMWF forecast logs from Baker's server
# Balcony log
wget --user=appalair --password=app1.air appalair@climate.appstate.edu:Everest_Forecasting-master/Data/ECMWF/ecaccess/Logs/fclog_Balcony.txt -O /home/lunet/gytm3/Everest2019/Forecasting/ECMWF/Logging/fclog_Balcony.txt

# South Col
wget --user=appalair --password=app1.air appalair@climate.appstate.edu:Everest_Forecasting-master/Data/ECMWF/ecaccess/Logs/fclog_SouthCol.txt -O /home/lunet/gytm3/Everest2019/Forecasting/ECMWF/Logging/fclog_SouthCol.txt

# CampII
wget --user=appalair --password=app1.air appalair@climate.appstate.edu:Everest_Forecasting-master/Data/ECMWF/ecaccess/Logs/fclog_CampII.txt -O /home/lunet/gytm3/Everest2019/Forecasting/ECMWF/Logging/fclog_CampII.txt



