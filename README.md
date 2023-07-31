# proteinANLYS_streamlit
***
## Set Up 
- build docker image
```sh
git clone https://github.com/jie151/proteinANLYS_streamlit.git
cd proteinANLYS_streamlit/docker/
docker build -t protein_anlys_streamlit:latest . 
docker run -p 8501:8501 protein_anlys_streamlit:latest
```
- pull docker image
```sh
還沒把最新的傳上去
```
***
## 檔案介紹(/docker/)
- **pic1_10_2.r**
畫DEP部分的 venn diagram
(另外放一個檔的原因: 用rpy2執行，要用rJava會有問題...)
- **pic2_10.r**
畫 18. running score and preranked list of GSEA result的圖
(另外放的原因: 用rpy2執行，無法儲存正確的檔案)
- **website.py**
畫其他的所有圖，並顯示
- **requirements.r, requirements.txt**
需要安裝的 R, python package
- **proteinGroups_HsinYuan_Rat.txt**
範例資料
