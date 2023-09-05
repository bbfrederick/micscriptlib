#!/bin/bash

MYIPADDRESS=`ifconfig en0 | grep 'inet ' | awk '{print $2}'`

#docker run -it micscriptlib bash

# allow network connections in Xquartz Security settings
xhost +
#xhost +si:localuser:micscriptlib

# Allow your local user access via xhost: xhost +SI:localuser:micscriptlib and create a similar user with docker run option: --user=$(id -u):$(id -g)
docker run \
    --network host\
    --volume=/Users/frederic:/Users/frederic \
    -it \
    -e DISPLAY=${MYIPADDRESS}:0 \
    -v /tmp/.X11-unix:/tmp/.X11-unix \
    -u micscriptlib fredericklab/micscriptlib:latest \
    PICAchooser \
        fix \
        --featdir /Users/frederic/Dropbox_PHC/MR_data/gradertest/079N_resting_visit1.feat \
        --melodicdir /Users/frederic/Dropbox_PHC/MR_data/gradertest/079N_resting_visit1.feat/filtered_func_data.ica \
        --scalemotiontodata
