FROM almalinux:9

LABEL maintainer="theo.reymermier@cern.ch"
LABEL description="MadGraph on AlmaLinux 9 with Python 3.9 and dependencies"

RUN dnf update -y && \
    dnf install -y \
    epel-release \
    python39 \
    python-devel \
    gcc \
    gcc-c++ \
    gcc-gfortran \
    make \
    wget \
    tar \
    which \
    libffi-devel \
    zlib-devel \
    bzip2 \
    git \
    perl \
    rsync \
    && dnf clean all

RUN ln -sf /usr/bin/python3.9 /usr/bin/python3 && \
    python3 --version

WORKDIR /External_Integration

RUN wget https://launchpad.net/mg5amcnlo/3.0/3.6.x/+download/MG5_aMC_v3.5.7.tar.gz && \
    tar -xvzf MG5_aMC_v3.5.7.tar.gz && \
    rm MG5_aMC_v3.5.7.tar.gz && \
    mv MG5_aMC_v3_5_7 MG5_aMC

RUN cd MG5_aMC && \
    echo "install lhapdf6" | ./bin/mg5_aMC && \
    echo "install pythia8" | ./bin/mg5_aMC && \
    echo "install mg5amc_py8_interface" | ./bin/mg5_aMC

ENTRYPOINT ["/bin/bash"]