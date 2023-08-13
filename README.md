# proteinANLYS_streamlit
***
## Set Up by docker
- build docker image
```sh
git clone https://github.com/jie151/proteinANLYS_streamlit.git
cd proteinANLYS_streamlit/docker/
docker build -t protein_anlys_streamlit:latest .
docker run -p 8501:8501 protein_anlys_streamlit:latest
```
- pull docker image
```
docker pull cguwebproteinanls/protein_anlys_stremalit:rscript
```
## Set Up (ubuntu:18.04)
```
ubuntu:18:04
python3.10.4
R 4.3.1, Bioconductor version 3.17 (BiocManager 1.30.21.1)
```
1. R (4.2~4.3都可以)
```
apt-get update && apt-get install -y --no-install-recommends build-essential&& apt-get install -y wget
apt update -qq
wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
apt install -y --no-install-recommends r-base
```

2. python3.10.4 (virtualenv)
```
add-apt-repository ppa:deadsnakes/ppa
apt update && apt upgrade -y
apt-get install -y  python3.10

apt install python3.10 python3-venv python3.10-venv python3-dev -y
python3.10 -m venv myvenv
```


3. 安裝python、R會用到的package
```
# 在虛擬環境中
apt-get install python3-pip -y
python3 -m pip install --upgrade pip
# 在requirements.txt 在docker資料夾中
pip3 install -r requirements.txt

apt-get update && upgrade -y
apt-get install -y r-base-core r-base-dev libcurl4-openssl-dev libssl-dev
apt-get install -y libcairo2-dev libxt-dev libnetcdf-dev libgdal-dev libnlopt-dev
apt install -y cmake
apt-get install -y default-jdk
# 在 docker 資料夾中
Rscript requirements.r
```

