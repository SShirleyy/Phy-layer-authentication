#!/bin/bash

cd /home/${USER}/prefix/default/
source setup_env.sh
cd /home/${USER}/

IFNAME=$(ls /sys/class/net | grep ^e)

if [-z $IFNAME]
then
		echo "no ethernet device found"
		exit 1
fi

if [ -n $(nmcli -terse dev show $IFNAME | grep GENERAL.STATE | grep '100\|70')]
then
	sudo nmcli dev disconnect $IFNAME
fi
sudo ifconfig $IFNAME 192.168.50.1
sudo ifconfig $IFNAME promisc

sudo sysctl -w net.core.rmem_max=50000000
sudo sysctl -w net.core.wmem_max=1048576

### create tap interface
if [[ `ifconfig -a | grep tap0 | wc -l` -eq 0 ]]
then
	sudo ip tuntap add dev tap0 user $USER mode tap
fi

### reconfigure it in any case, just to be sure it's up
sudo ifconfig tap0 down
sudo ifconfig tap0 hw ether 12:34:56:78:90:ab
sudo ifconfig tap0 up
sudo ifconfig tap0 192.168.123.1


FIFO="/tmp/wifi.pcap"

if [ -e $FIFO ]
then
    sudo rm -rf $FIFO
fi
mkfifo $FIFO
