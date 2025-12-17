FROM ubuntu:22.04

# Set environment variables for the build process
ENV MINICONDA_PATH=/opt/conda
ENV PATH=$MINICONDA_PATH/bin:$PATH
ENV CC3D_ENV=cc3d_compile

# Create a dedicated user for security
RUN adduser --disabled-password --gecos "" cc3duser

# Install dependencies and download the latest Miniconda
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    wget \
    bzip2 \
    ca-certificates \
    libxml2-dev \
    libxslt-dev \
    libgl1-mesa-glx \
    libgl1-mesa-dev \
    git \
    build-essential \
    cmake \
    patchelf \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Download Miniconda and install
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py310_25.7.0-2-Linux-x86_64.sh -O miniconda.sh \
    && /bin/bash miniconda.sh -b -p $MINICONDA_PATH \
    && rm miniconda.sh

# Set file permissions for the new user
# The user needs to own the /opt/conda directory to manage packages
RUN chown -R cc3duser:cc3duser $MINICONDA_PATH

# Use the new user for all subsequent commands
USER cc3duser

WORKDIR /src-cc3d

# Accept channels' Terms of Service
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main \
    && conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

RUN conda init

# Install Mamba
RUN conda install --yes mamba -c conda-forge \
    && conda clean --all --yes

COPY --chown=cc3duser:cc3duser environment.yaml environment.yaml

# Create a new conda environment and install packages
RUN mamba env create -f environment.yaml --name $CC3D_ENV

# Clean up cache files to reduce the image size
RUN mamba clean --all --yes

# Clone the sources
RUN git clone --branch 4.6.0.1 --depth 1 https://github.com/CompuCell3D/CompuCell3D.git \
    && git clone https://github.com/CompuCell3D/cc3d-player5.git \
    && git clone https://github.com/CompuCell3D/cc3d-twedit5.git

COPY DeveloperZone /src-cc3d/CompuCell3D/CompuCell3D/DeveloperZone

# Let conda to be launched without having to do `conda activate`
SHELL ["conda", "run", "-n", "cc3d_compile", "/bin/bash", "-c"]

RUN echo "conda activate $CC3D_ENV" >> ~/.bashrc

# Demonstrate the environment is activated
# RUN echo "Make sure pandas is installed:"
# RUN python -c "import pandas"

RUN cmake -S /src-cc3d/CompuCell3D/CompuCell3D \
    -B /src-cc3d/CompuCell3D_build \
    -DPython3_EXECUTABLE=$MINICONDA_PATH/envs/$CC3D_ENV/bin/python \
    -DNO_OPENCL=ON  \
    -DBUILD_STANDALONE=OFF \
    -DOPENGL_gl_LIBRARY=/usr/lib64/libGL.so \
    -DOPENGL_glx_LIBRARY=/usr/lib64/libGLX.so \
    -G "Unix Makefiles" \
    -DCMAKE_INSTALL_PREFIX=/src-cc3d/CompuCell3D_install

WORKDIR /src-cc3d/CompuCell3D_build

RUN make -j $(nproc) && make install

RUN ln -s /src-cc3d/cc3d-twedit5/cc3d/twedit5 $MINICONDA_PATH/envs/$CC3D_ENV/lib/python3.10/site-packages/cc3d/twedit5 \
    && ln -s /src-cc3d/cc3d-player5/cc3d/player5 $MINICONDA_PATH/envs/$CC3D_ENV/lib/python3.10/site-packages/cc3d/player5

WORKDIR /src-cc3d

# Custom Python script to avoid to run twedit5 GUI for configuring DeveloperZone
COPY --chown=cc3duser:cc3duser configure_dev_zone.py configure_dev_zone.py

USER root

# 1. On force bash pour la commande "source"
RUN ln -sf /bin/bash /bin/sh

# 2. On crée un lien symbolique pour le script "activate" manquant dans l'env
RUN ln -s /opt/conda/bin/activate /opt/conda/envs/cc3d_compile/bin/activate

USER cc3duser

WORKDIR /src-cc3d/CompuCell3D_dev_zone_build

# On repasse en utilisateur normal
USER cc3duser

# 1. On crée le dossier de build explicitement
WORKDIR /src-cc3d/CompuCell3D_dev_zone_build

# 2. On lance CMake DIRECTEMENT (sans le script Python fragile)
# On reprend exactement les arguments que le script Python tentait d'utiliser (visibles dans vos logs précédents)
RUN cmake -G "Unix Makefiles" \
    -DCMAKE_BUILD_TYPE:STRING=RelWithDebInfo \
    -DCMAKE_INSTALL_PREFIX:PATH=$MINICONDA_PATH/envs/$CC3D_ENV/lib/python3.10/site-packages \
    -DCOMPUCELL3D_GIT_DIR:PATH=/src-cc3d/CompuCell3D \
    -DCOMPUCELL3D_INSTALL_PATH:PATH=$MINICONDA_PATH/envs/$CC3D_ENV/lib/python3.10/site-packages \
    -DPYTHON_INCLUDE_DIR:PATH=$MINICONDA_PATH/envs/$CC3D_ENV/include/python3.10 \
    -DPYTHON_LIBRARY:PATH=$MINICONDA_PATH/envs/$CC3D_ENV/lib/libpython3.10.so \
    -S /src-cc3d/CompuCell3D/CompuCell3D/DeveloperZone \
    -B .

# 3. Maintenant le Makefile existe, on peut compiler
RUN make -j $(nproc) && make install

# RUN make && make install

RUN patchelf --set-rpath "$MINICONDA_PATH/envs/cc3d_compile/lib:$MINICONDA_PATH/envs/cc3d_compile/lib/python3.10/site-packages/cc3d/cpp/CompuCell3DPlugins:$MINICONDA_PATH/envs/cc3d_compile/lib/python3.10/site-packages/cc3d/cpp/CompuCell3DSteppables" \
    $MINICONDA_PATH/envs/cc3d_compile/lib/python3.10/site-packages/cc3d/cpp/_CompuCellExtraModules.so

# New directory for running custom scripts
WORKDIR /script_cc3d

# COPY --chown=cc3duser:cc3duser script_cc3d /script_cc3d

CMD ["/bin/bash"]
