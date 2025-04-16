#!/bin/bash

UBUNTU_DEPENDENCIES="binutils cmake dpkg-dev g++ gcc libssl-dev git libx11-dev libxext-dev libxft-dev libxpm-dev python3 libtbb-dev"
CENTOS_DEPENDENCIES="git make cmake gcc-c++ gcc binutils libX11-devel libXpm-devel libXft-devel libXext-devel python openssl-devel xrootd-client-devel xrootd-libs-devel"

if [ -f /etc/os-release ]; then
    . /etc/os-release
    OS=$ID
else
    echo "Impossible de détecter le système d'exploitation."
    exit 1
fi

if [ "$OS" == "ubuntu" ]; then
    DEPENDENCIES=$UBUNTU_DEPENDENCIES
    INSTALL_CMD="sudo apt-get install -y"
elif [ "$OS" == "centos" ]; then
    DEPENDENCIES=$CENTOS_DEPENDENCIES
    INSTALL_CMD="sudo yum install -y"
else
    echo "Système d'exploitation non supporté."
    exit 1
fi

for PACKAGE in $DEPENDENCIES; do
    if ! dpkg -s $PACKAGE &>/dev/null && ! rpm -q $PACKAGE &>/dev/null; then
        echo "Le package $PACKAGE n'est pas installé."
        echo "Veuillez exécuter la commande suivante pour l'installer :"
        echo "$INSTALL_CMD $PACKAGE"
        exit 1
    fi
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="$SCRIPT_DIR/Common/config.py"

if grep -q "root_path" "$CONFIG_FILE"; then
    sed -i "s|^root_path.*|root_path = '$SCRIPT_DIR'|" "$CONFIG_FILE"
else
    echo "root_path = '$SCRIPT_DIR'" >> "$CONFIG_FILE"
fi

echo "Configuration mise à jour dans $CONFIG_FILE"
