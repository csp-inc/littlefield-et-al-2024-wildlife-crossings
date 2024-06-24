#!/usr/bin/env bash
sudo docker run -it --rm \
	-v /home/azureuser/datadrive/pew-connectivity:/home/omniscape \
	-w /home/omniscape \
	-e JULIA_NUM_THREADS=8 \
	vlandau/omniscape:0.5.3
