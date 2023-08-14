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
R 4.3.1, Bioconductor version 3.17
```
1. R (4.2~4.3都可以)
```
sudo apt update && apt install -y --no-install-recommends build-essential&& apt install wget
sudo apt update -qq
wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
sudo apt install -y --no-install-recommends r-base
```

2. python3.10.4 (virtualenv)
```
sudo apt update && sudo apt install zlib1g-dev libncurses5-dev libgdbm-dev libnss3-dev libssl-dev libreadline-dev libffi-dev libsqlite3-dev libbz2-dev pkg-config -y
wget https://www.python.org/ftp/python/3.10.4/Python-3.10.4.tgz
tar -xf Python-3.10.*.tgz
cd Python-3.10.4/
./configure --enable-optimizations
make -j $(nproc)
sudo make altinstall

sudo apt install python3-venv python3-dev python3-pip -y
python3.10 -m venv myvenv
source myvenv/bin/activate
```

3. 安裝python、R會用到的package
```
git clone https://github.com/jie151/proteinANLYS_streamlit.git
cd proteinANLYS_streamlit/docker/
# virtualenv
python3 -m pip install --upgrade pip
# 在requirements.txt 在docker資料夾中
pip3 install -r requirements.txt

sudo apt update && upgrade -y
sudo apt install -y r-base-core r-base-dev libcurl4-openssl-dev libssl-dev
sudo apt install -y libcairo2-dev libxt-dev libnetcdf-dev libgdal-dev libnlopt-dev
sudo apt install -y cmake default-jdk
# 在 docker 資料夾中
Rscript requirements.r
```

4. streamlit run (在 proteinANLYS_streamlit/docker/ 執行)
```
streamlit run website.py
```

