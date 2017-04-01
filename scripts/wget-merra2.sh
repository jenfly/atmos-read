#! /bin/bash


cd /net/eady/data1/jwalker/datastore/merra2/wget
filenm=/net/eady/data1/jwalker/dynamics/python/atmos-read/scripts/merra_urls/merra_urls_ubudget.txt
wget  --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookie --no-check-certificate --content-disposition -i $filenm

# ------------------------------
#cd /home/jwalker/datastore/merra2/dailyrad/
#
#for yr in {1980..2015}
#do
#    mkdir $yr
#    cd $yr
#    filenm=/home/jwalker/dynamics/python/atmos-read/scripts/merra_urls/merra2_rad_$yr.txt
#    wget  --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookie --no-check-certificate --content-disposition -i $filenm
#    cd ..
#done
# ------------------------------

# url="http://goldsmr4.gesdisc.eosdis.nasa.gov/daac-bin/OTF/HTTP_services.cgi?FILENAME=%2Fdata%2FMERRA2%2FM2T1NXRAD.5.12.4%2F1980%2F01%2FMERRA2_100.tavg1_2d_rad_Nx.19800101.nc4&FORMAT=bmM0Lw&BBOX=-65%2C40%2C65%2C120&LABEL=MERRA2_100.tavg1_2d_rad_Nx.19800101.SUB.nc4&FLAGS=1&SHORTNAME=M2T1NXRAD&SERVICE=SUBSET_MERRA2&LAYERS=&VERSION=1.02&VARIABLES=lwgnt%2Clwtup%2Cswgnt%2Cswtnt"
#
# echo $url
#
# wget  --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookie --no-check-certificate --content-disposition $url
#
# ------------------------------------------------------------

# urldir=http://goldsmr4.sci.gsfc.nasa.gov/opendap/hyrax/MERRA2_MONTHLY/M2TMNXRAD.5.12.4/
# datadir=/home/jwalker/datastore/merra2/monthly/
#
# for yr in {2011..2015}
# do
#     for m in 01 02 03 04 05 06 07 08 09 10 11 12
#     do
#         filenm=MERRA2_400.tavgM_2d_rad_Nx.$yr$m.nc4
#         url=$urldir$yr/$filenm.nc4
#         echo $url
#         wget  --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookie --no-check-certificate $url
#         mv $filenm.nc4 $datadir/$filenm
#     done
# done

# --------------------------------------------------

# url=http://goldsmr4.sci.gsfc.nasa.gov/opendap/hyrax/MERRA2_MONTHLY/M2TMNXFLX.5.12.4/1980/MERRA2_100.tavgM_2d_flx_Nx.198001.nc4.nc4
#
# wget  --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookie --no-check-certificate http://goldsmr4.sci.gsfc.nasa.gov/opendap/hyrax/MERRA2_MONTHLY/M2TMNXFLX.5.12.4/1980/MERRA2_100.tavgM_2d_flx_Nx.198001.nc4.nc4
#
#
# wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookies $url

#wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --auth-no-challenge=on --keep-session-cookies http://server[:port]/path/file[.format[?subset]] -O output[.format]

# For an entire Directory
#wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --auth-no-challenge=on --keep-session-cookies -npr http://server[:port]/path/file[.format[?subset]] -O output[.format]
