FROM --platform=linux/amd64 rockylinux:8

# Ref: https://docs.conda.io/projects/conda/en/latest/user-guide/install/rpm-debian.html
# Import our GPG public key
RUN rpm --import https://repo.anaconda.com/pkgs/misc/gpgkeys/anaconda.asc

# Add the Anaconda repository
RUN <<EOF
echo "[conda]" >> /etc/yum.repos.d/conda.repo
echo "name=Conda" >> /etc/yum.repos.d/conda.repo
echo "baseurl=https://repo.anaconda.com/pkgs/misc/rpmrepo/conda" >> /etc/yum.repos.d/conda.repo
echo "enabled=1" >> /etc/yum.repos.d/conda.repo
echo "gpgcheck=1" >> /etc/yum.repos.d/conda.repo
echo "gpgkey=https://repo.anaconda.com/pkgs/misc/gpgkeys/anaconda.asc" >> /etc/yum.repos.d/conda.repo
EOF

COPY prokevo.yml /opt

ENV PATH="/opt/conda/bin:${PATH}"

RUN dnf install -y conda openssh-clients python3 findutils epel-release && \
    dnf install -y --enablerepo=powertools R && \
    cp /opt/conda/etc/profile.d/conda.sh /etc/profile.d && \
    mkdir -p /opt/ProkEvo_dir/prokevo && \
    conda env create -f /opt/prokevo.yml -p /opt/ProkEvo_dir/prokevo
