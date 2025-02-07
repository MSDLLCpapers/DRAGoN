#!docker build

# syntax=docker/dockerfile-upstream:master-labs
# hadolint global ignore=SC1091

ARG FROM=amazonlinux:2023.6.20250128.0

FROM $FROM

LABEL       dock.img.name="scottnortonphd/msdllcpapers/dragon" \
            dock.img.description="AmazonLinux2023-based image with DRAGoN pipeline" \
            \
            dock.maintainer.isid="nortonsc" \
            dock.maintainer.name="Scott Norton" \
            dock.maintainer.email="scott.norton@merck.com" \
            dock.maintainer.division="Merck Research Labs IT" \
            dock.maintainer.team="DSSI Bioinformatics" \
            \
            dock.docker.run="docker container run --rm -it  IMAGE" \
            \
            dock.os="linux" \
            org.opencontainers.image.source="https://github.com/MSDLLCPapers/DRAGoN" \
            org.opencontainers.image.description="AmazonLinux2023-based image with DRAGoN pipeline" \
            org.opencontainers.image.licenses="MIT"

SHELL ["/bin/bash", "-c"]

# ADD proxy.conf /etc/yum.repos.d
# COPY yum.conf /etc/yum.repos.d

RUN update-ca-trust enable && update-ca-trust extract \
 && yum update -y && yum install -y --nogpgcheck --allowerasing \
    sudo-1.9.15-1.p5.amzn2023.0.1 \
    curl-8.5.0-1.amzn2023.0.4 \
    wget-1.21.3-1.amzn2023.0.4 \
    git-2.47.1-1.amzn2023.0.2 \
    glibc-langpack-en-2.34-117.amzn2023.0.1 \
    gnupg2-2.3.7-1.amzn2023.0.4 \
    procps-ng-3.3.17-1.amzn2023.0.2 \
    libxcrypt-4.4.33-7.amzn2023 \
    bc-1.07.1-14.amzn2023.0.2 \
    tar-1.34-1.amzn2023.0.4 \
    findutils-1:4.8.0-2.amzn2023.0.2 \
    which-2.21-26.amzn2023.0.2 \
    unzip-6.0-57.amzn2023.0.2 \
    && yum clean -y all
#RUN rm -rf /var/cache/yum

ENV LANG="en_US.UTF-8" LANGUAGE="en_US:en" LC_ALL="en_US.UTF-8"

## Install miniforge3
ENV ANACONDA3_VERSION=24.11.3-0 ANACONDA3_ARCH=Linux-x86_64
ADD https://github.com/conda-forge/miniforge/releases/download/${ANACONDA3_VERSION}/Miniforge3-${ANACONDA3_VERSION}-${ANACONDA3_ARCH}.sh /Miniforge3.sh
RUN bash Miniforge3.sh -b -p "/opt/conda" \
 && bash "/opt/conda/etc/profile.d/conda.sh"

# Create conda environment
ENV CONDA_DEFAULT_ENV="drugseq-env"
ENV PATH=$PATH:/opt/conda/bin
COPY conf/conda /opt/envs
# shellcheck disable=SC1091
RUN conda env create -f /opt/envs/$CONDA_DEFAULT_ENV.yml \
 && source activate ${CONDA_DEFAULT_ENV} \
 && echo "source activate $CONDA_DEFAULT_ENV" >> ~/.bashrc
ENV PATH=/opt/conda/envs/${CONDA_DEFAULT_ENV}/bin:$PATH
#RUN conda clean -y --all

# Compile C++ utilities
COPY src /src
COPY CMakeLists.txt /CMakeLists.txt
COPY version.txt /version.txt
ENV DRUGSEQ_INSTALL_PREFIX="/opt/drugseq"
ENV DRUGSEQ_BIN_DIR="$DRUGSEQ_INSTALL_PREFIX/bin"
ENV PATH="$PATH:$DRUGSEQ_BIN_DIR"
ENV CONDA_PREFIX=/opt/conda/envs/${CONDA_DEFAULT_ENV}
ENV LD_LIBRARY_PATH=${CONDA_PREFIX}/lib
RUN mkdir -p build \
 && cmake . -B build \
        -DCMAKE_INSTALL_PREFIX=${DRUGSEQ_INSTALL_PREFIX} \
        -DCMAKE_PREFIX_PATH=${CONDA_PREFIX} \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_MODULE_PATH=${CONDA_PREFIX}/unpacked_source/cmake \
        -DDRUGSEQ_VERSION_STR="$(cat /version.txt)" \
 && cmake --build build --parallel 4 \
 && cmake --install build --strip

# install awscli v2
RUN curl -s "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "/tmp/awscliv2.zip" \
 && unzip -q /tmp/awscliv2.zip -d /tmp \
 && /tmp/aws/install -b /usr/bin \
 && rm -rf /tmp/aws*

RUN hash -r

# Setup entrypoint
CMD ["/bin/bash"]
