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
還沒把最新的傳上去
```
## Set Up (ubuntu:18.04)
```
ubuntu:18:04
python3.10.4
R 4.3.1, Bioconductor version 3.17 (BiocManager 1.30.21.1)
```

1. python3.10.4 (virtualenv)
```
add-apt-repository ppa:deadsnakes/ppa && apt-get update
apt-get install -y  python3.10

apt install python3.10 python3-venv python3.10-venv python3-dev -y
python3.10 -m venv myvenv
activate PATH=/myvenv/bin:$PATH

apt-get install python3-pip -y
python3 -m pip install --upgrade pip

pip3 install -r requirements.txt
```
2. R 4.3.1
```

```

2.
***
<!-- ## 檔案介紹(/docker/)
- **pic1_10_2.r**</br>
畫DEP部分的 venn diagram</br>
(另外放一個檔的原因: 用rpy2執行，要用rJava會有問題...)</br>
- **pic2_10.r**</br>
畫 18. running score and preranked list of GSEA result的圖</br>
(另外放的原因: 用rpy2執行，無法儲存正確的檔案)</br>
- **website.py**</br>
畫其他的所有圖，並顯示</br>
- **requirements.r, requirements.txt**</br>
需要安裝的 R, python package</br>
- **proteinGroups_HsinYuan_Rat.txt**</br>
範例資料</br> -->
