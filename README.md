# proteinANLYS_streamlit
***
# Set Up 
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
# 檔案介紹(/docker/)
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
範例資料</br>
