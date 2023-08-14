import os
# from rpy2 import robjects
import streamlit as st
import pandas as pd
from PIL import Image
from contextlib import contextmanager, redirect_stdout
from io import StringIO
import base64
import uuid
import re
import sys
import numpy as np
import statistics
from streamlit.runtime.scriptrunner.script_run_context import get_script_run_ctx
import matplotlib.pyplot as plt


# 在網頁顯示print的內容
@contextmanager
def st_capture(output_func):
    with StringIO() as stdout, redirect_stdout(stdout):
        old_write = stdout.write
        def new_write(string):
            ret = old_write(string)
            output_func(stdout.getvalue())
            return ret
        stdout.write = new_write
        yield

def download_button(download_filename, button_text):
    # Load file
    with open(download_filename, 'rb') as f:
        object_to_download = f.read()
    b64 = base64.b64encode(object_to_download).decode()

    button_uuid = str(uuid.uuid4()).replace('-', '')
    button_id = re.sub('\d+', '', button_uuid)

    custom_css = f"""
        <style>
            #{button_id} {{
                background-color: rgb(255, 255, 255);
                color: rgb(38, 39, 48);
                padding: 0.25em 0.38em;
                position: relative;
                text-decoration: none;
                border-radius: 4px;
                border-width: 1px;
                border-style: solid;
                border-color: rgb(230, 234, 241);
                border-image: initial;
            }}
            #{button_id}:hover {{
                border-color: rgb(246, 51, 102);
                color: rgb(246, 51, 102);
            }}
            #{button_id}:active {{
                box-shadow: none;
                background-color: rgb(246, 51, 102);
                color: white;
                }}
        </style> """

    dl_link = custom_css + f'<a download="{download_filename}" id="{button_id}" href="data:file/txt;base64,{b64}">{button_text}</a><br></br>'
    st.markdown(dl_link, unsafe_allow_html=True)

def save_uploadedfile(uploadedfile):
    with open(os.path.join("uploadFile",uploadedfile.name),"wb") as f:
        f.write(uploadedfile.getbuffer())
# 使用者上傳檔案，若有會儲存起來、無則使用預設檔案
def upload_file():
    if not os.path.exists("./file/"):
        os.system("mkdir file")
    if not os.path.exists(f"./file/{session_id}/"):
        # st.write(session_id)
        os.system(f"mkdir ./file/{session_id}/")
    # 上傳檔案
    uploaded_file = st.sidebar.file_uploader('Upload a TXT/CSV file', type=['csv', 'txt'], accept_multiple_files=False)
    if uploaded_file is not None:
        # 將上傳的檔案儲存下來 (由R開啟)
        filename = f"./file/{session_id}/uploadFile.txt"
        with open(filename,"wb") as f:
            f.write( uploaded_file.getbuffer())

        sep_word = "," if uploaded_file.type == "text/csv" else "\t"
        data = pd.read_csv(uploaded_file, sep = sep_word, error_bad_lines=False, dtype=object)
    else:
        filename = "./proteinGroups_HsinYuan_Rat.txt"
        data = pd.read_csv(filename, sep="\t", dtype=object)
    data.columns = data.columns.str.replace(' ', '.')
    return filename, data

# 刪除上一次舊的資料與結果
@st.cache_data
def delete_file(data):
    if os.path.exists(f"./file/{session_id}/image"):
        os.system(f"rm -r ./file/{session_id}/image")

    os.system(f"mkdir file ./file/{session_id}/image")


# 從檔案中的Fasta.headers判斷資料的物種，並看有幾種condition, 讓使用者選擇要分析哪一種condition
# return species, sorted(condition_py), numCondition_py
def initialize_DEP_data(data):
    # 確認資料的物種
    species_dict = {'mouse': "OX=10090", "rat":"OX=10116", "human":"OX=9606"}
    species = "human"
    try:
        species = [key for key, value in species_dict.items() if re.search(value, data.iloc[1]["Fasta.headers"])][0]
    except:
        print("確認物種是否為: mouse, rat, human") # warning
    # 確認資料包含哪些condition
    data = data.filter(regex=("Reporter.intensity.corrected.*"))

    condition_py = set([ re.sub("Reporter.intensity.corrected.[0-9]+", "",colname) for colname in data.columns])
    numCondition_py = len(condition_py)
    return species, sorted(condition_py), numCondition_py

# 選定condition後，畫出該condition中每個欄位各自的分布，計算median、outlier，讓使用者選擇要留下的欄位
# return data_exp_design_dict, data
def select_condition_median_outlier(data, option):
    # 選定condition後，畫出該condition中每個欄位各自的分布，計算median、outlier，讓使用者選擇要留下的欄位
    data = data.filter(regex=("Reporter.intensity.corrected.*"))
    if option:
        data = data.filter(regex=(option))
    data = data.apply(pd.to_numeric)
    data = np.log2(data)

    data_exp_design_dict = {}
    data_exp_design_dict['colnames'] = data.columns
    data_exp_design_dict['medians'] = [statistics.median(data[colname]) for colname in data_exp_design_dict['colnames']]
    median_list_isfinite = [ median for median in data_exp_design_dict['medians']  if np.isfinite(median) ]

    temp = pd.DataFrame(median_list_isfinite).describe()
    iqr = temp.loc['75%'] - temp.loc['25%']
    outlier_up = temp.loc['75%'] + 1.5*iqr
    outlier_down = temp.loc['25%'] - 1.5*iqr
    data_exp_design_dict['outliers'] = [ median for median in data_exp_design_dict['medians'] if (median > outlier_up[0] or median < outlier_down[0]) ]

    return data_exp_design_dict, data

# 讀取 select_condition_median_outlier(data, option) 產生的data_exp_design_dict，並畫出分布圖
# return selected_col_id_list 選擇的欄位編號
def draw_distogram_after_select_condition(rowdata_log2, data_exp_design_dict):
    data_colname = data_exp_design_dict['colnames']
    median_list  = data_exp_design_dict['medians']
    outlier_list = data_exp_design_dict['outliers']
    with st.expander("Select the set of columns you want to analyze"):
        # 把圖分成有選擇與沒選擇分別顯示
        cols = st.columns(3)
        # 有選擇的欄位編號
        selected_col_id_list = []
        for i in range(0, len(data_colname)):
            check_value = True if median_list[i] not in outlier_list  else False

            colname_checkobx = cols[ i % 3 ].checkbox(data_colname[i], value=check_value)
            if colname_checkobx :
                selected_col_id_list.append(i)

            data_isfinite = rowdata_log2[np.isfinite(rowdata_log2[data_colname[i]])]

            cols[i % 3].text(f"median: {round(median_list[i], 3)}")

            data_distri_fig, ax = plt.subplots()
            ax.hist(data_isfinite[data_colname[i]], bins=10)
            cols[i % 3].pyplot(data_distri_fig)

        plt.close('all')

        return selected_col_id_list

# 選定要留下的欄位後，先自動幫欄位分組，分組方式: 從A開始編號，從2開始看是否能整除所選的欄位數量
# return condition_replicate_list (type: list of list)
def generate_default_group(num_selected_col):
    condition_replicate_list = []
    num_category = 1

    for i in range(2, num_selected_col):
        if num_selected_col % i == 0:
            num_category = i # 找到要分幾組 (condition)
            break
    a_category_colnum = int(num_selected_col / num_category) # 一組有幾個欄位 (replicate)
    for i in range(0, num_category):
        for j in range(1, a_category_colnum+1):
            # 從A開始編號
            condition_replicate_list.append([chr(65+i), j])
    return condition_replicate_list

# 確保字串符合規則 (只由字母、數字、dot、底線組成，且開頭為字母或dot後面不加數字)
# 參考: R的function => make.names: Make Syntactically Valid Names
def make_name(string):
    # 只允許有 letters, numbers and the dot or underline characters
    string = re.sub(r'[^a-zA-Z0-9._]', '.', string)
    # starts with a letter or the dot not followed by a number
    allowed_pattern = r'([a-zA-Z][a-zA-Z0-9_]*|\.[^0-9].*)'
    match = re.match(allowed_pattern, string)
    if not match:
        string = "X" + string
    return string

# 根據 selected_col_id_list，呼叫 generate_default_group，預設分組與輸入值，生成 exp_design.csv並顯示在網頁
# return exp_design (選擇欄位的 condition, replicate)
def generate_exp_design_inputCondition(selected_col_id_list, data_colname):
    exp_design = pd.DataFrame(columns=["state","label","condition","replicate"], index=range(1, 11)).fillna("N") #創一個dataframe, 預設填N
    # 根據欄位數，預設分組與輸入值
    condition_replicate_list = generate_default_group(len(selected_col_id_list))

    with st.sidebar.expander("Select property columns"):
        for index, col_id in enumerate(selected_col_id_list):
            exp_design.at[col_id + 1, 'state'] = "T"
            exp_design.at[col_id + 1, 'label'] = data_colname[col_id]
            st.write(data_colname[col_id])
            condition_textInput = st.text_input('condition', key= f"condition{col_id}", value=condition_replicate_list[index][0])
            replicate_textInput = st.number_input('replicate', key= f"replicate{col_id}", value=condition_replicate_list[index][1], step=1, format="%d")

            if condition_textInput:
                condition_textInput = make_name(condition_textInput)
                exp_design.at[col_id + 1, 'condition'] = condition_textInput
            if replicate_textInput:
                exp_design.at[col_id + 1, 'replicate'] = replicate_textInput

    exp_design = exp_design[(exp_design.state != "N")]
    exp_design = exp_design.drop("state", axis=1)
    exp_design.to_csv(f"./file/{session_id}/exp_design.csv", index=False)
    st.write("exp_design: ",exp_design)

    return exp_design

# 執行 dep1_to_1_5.r 畫DEP部分的第一張圖至第五張圖，並顯示
@st.cache_data
def cache_DEP_data1(filename_py, session_id, colname_geneNames_py, colname_proteinIDs_py, filer_missing_values, normalizeOption_py):
    os.system(f"Rscript dep1_to_1_5.r {filename_py} {session_id} {colname_geneNames_py} {colname_proteinIDs_py} {normalizeOption_py} {filer_missing_values['option']} {filer_missing_values['nThr']} {filer_missing_values['min']}")

    # plot1_1
    try:
        st.write(" *Plot a barplot of the protein identification overlap between samples")
        st.image(Image.open(f'./file/{session_id}/image/plot1_1.png'))
    except Exception as e:
        st.error(f"error! {e}")
    # plot1_2
    try:
        st.write(" *Plot a barplot of the number of identified proteins per samples")
        st.image(Image.open(f'./file/{session_id}/image/plot1_2.png'))
    except Exception as e:
        st.error(f"error! {e}")
    # plot1_3
    try:
        st.write(" *Plot a barplot of the protein identification overlap between samples")
        st.image(Image.open(f'./file/{session_id}/image/plot1_3.png'))
    except Exception as e:
        st.error(f"error! {e}")
    # plot1_4
    try:
        st.header("3. Normalization")
        st.write(" Visualize normalization by boxplots for all samples before and after normalization")
        st.image(Image.open(f'./file/{session_id}/image/plot1_4.png'))
    except Exception as e:
        st.error(f"error! {e}")
    # plot1_5
    image_label_list = ["*Plot a heatmap of proteins with missing values",
                        "*Plot intensity distributions and cumulative fraction of proteins with and without missing values",
                        "*Plot intensity distributions before and after imputation"]
    st.header("4. Impute data for missing values: ")
    for i in range(3):
        try:
            st.write(image_label_list[i])
            st.image(Image.open(f'./file/{session_id}/image/plot1_5_{i+1}.png'))
        except Exception as e:
            st.error(f"error! {e}")

# 執行  dep2_1_6_to_1_8_1_10.r 畫DEP部分的第六、八、十張圖並顯示 (1_9因為可能會一直改變所以另外畫)
@st.cache_data
def cache_DEP_data2(filename_py, session_id, colname_geneNames_py, colname_proteinIDs_py, filer_missing_values, normalizeOption_py, control_py, alpha_py, lfc_py, contrast_py):
    os.system(f"Rscript dep2_1_6_to_1_8_1_10.r {filename_py} {session_id} {colname_geneNames_py} {colname_proteinIDs_py} {normalizeOption_py} {filer_missing_values['option']} {filer_missing_values['nThr']} {filer_missing_values['min']} {control_py} {contrast_py} {alpha_py} {lfc_py} ")
    # plot1_6
    image_label_list = ["*Plot the first and second principal components", "*Plot the Pearson correlation matrix"]
    for i in range(1, 3):
        try:
            st.write(image_label_list[i-1])
            st.image(Image.open(f'./file/{session_id}/image/plot1_6_{i}.png'))
        except Exception as e:
            st.error(f"error! {e}")
    # plot1_7
    image_label_list = ["*Plot a heatmap of all significant proteins with the data centered per protein",
                        "*Plot a heatmap of all significant proteins (rows) and the tested contrasts (columns)"]
    for i in range(1, 3):
        try:
            st.write(image_label_list[i-1])
            st.image(Image.open(f'./file/{session_id}/image/plot1_7_{i}.png'))
        except Exception as e:
            st.error(f"error! {e}")

    # plot1_8
    st.header("6. Volcano plots of specific contrasts:")
    try:
        st.image(Image.open(f"./file/{session_id}/image/plot1_8.png"))
    except Exception as e:
        st.error(f"error! {e}")

# 執行  dep3_1_9.r 畫DEP部分的第九張圖並顯示
@st.cache_data
def cache_DEP_data3(filename_py, session_id, colname_geneNames_py, colname_proteinIDs_py, filer_missing_values, normalizeOption_py, control_py, alpha_py, lfc_py, num_protein, proteinData_py):
    protein_type = "contrast" if num_protein > 1 else "centered"

    if len(proteinData_py) == len(set(proteinData_py)):
        try:
            with open(f"./file/{session_id}/proteinData.txt", "w") as f:
                for protein in proteinData_py:
                    f.write(protein)
                    f.write("\n")
            os.system(f"Rscript dep3_1_9.r {filename_py} {session_id} {colname_geneNames_py} {colname_proteinIDs_py} {normalizeOption_py} {filer_missing_values['option']} {filer_missing_values['nThr']} {filer_missing_values['min']} {control_py} {alpha_py} {lfc_py} {protein_type}")
            st.image(Image.open(f"./file/{session_id}/image/plot1_9.png"))
        except Exception as e:
            st.error(str(e), icon="🚨")
    else:
        st.error("Error!! The proteins selected cannot be duplicated", icon="🚨")

# 執行 dose1_config.r，做物種、名稱的轉換，將結果存為 dep_output_result.csv
@st.cache_data
def DOSE_data_config(exp_design, filer_missing_values, normalizeOption_py, control_py, alpha_py, lfc_py, session_id, species, ratioName_py):
    os.system(f"Rscript dose1_config.r {session_id} {species} {ratioName_py}")
    try:
        with open(f"./file/{session_id}/dep_output_result.csv", "r") as file:
            print("dep_output_result.csv created!")
    except Exception as e:
        st.error(f'與uriprot.org網站連接失敗(Entrez -> GeneID), 請重新載入網頁, {str(e)}', icon="🚨")
        sys.exit()

# 執行 dose2_plot.r， 畫DOSE部分的第1~7、9-10張圖 (1_8因為可能會有pvalue_2_8_py參數，另外畫)
@st.cache_data
def draw_DOSE_pic(exp_design, filer_missing_values, normalizeOption_py, control_py, alpha_py, lfc_py, ratioName_py, de_up_down_py, range_py, enrichment_analysis_methods_py, universal_enrichment_category_py, universal_enrichment_subcategory_py):
    os.system(f"Rscript dose2_plot.r {session_id} {ratioName_py} {de_up_down_py} {range_py} {enrichment_analysis_methods_py} {universal_enrichment_category_py} {universal_enrichment_subcategory_py}")

    st.header("10. Bar Plot (DOSE)")
    try:
        st.image(Image.open(f"./file/{session_id}/image/plot2_1.png"))
    except:
        st.error("No significant terms were enriched")
    st.header("11. Dot plot")
    try:
        st.image(Image.open(f"./file/{session_id}/image/plot2_2_1.png"))
        st.image(Image.open(f"./file/{session_id}/image/plot2_2_2.png"))
    except:
        st.error("No significant terms were enriched")
    st.header("12. Gene-Concept Network")
    try:
        for i in range(1, 6):
            st.image(Image.open(f"./file/{session_id}/image/plot2_3_{i}.png"))
    except:
        st.error("No significant terms were enriched")
    st.header("13. Heatmap-like functional classification")
    try:
        for i in range(1, 3):
            st.image(Image.open(f'./file/{session_id}/image/plot2_4_{i}.png'))
            download_button(f"./file/{session_id}/image/heatplot_{i}.png", f"Download heatplot_{i}.png")
    except:
        st.error("No significant terms were enriched")
    st.header("14. Enrichment Map")
    try:
        for i in range(1,5):
            st.image(Image.open(f"./file/{session_id}/image/plot2_5_{i}.png"))
    except:
        st.error("No significant terms were enriched")
    st.header("15. Biological theme comparison")
    try:
        for i in range(1, 5):
            st.image(Image.open(f"./file/{session_id}/image/plot2_6_{i}.png"))
    except Exception as e:
        st.error(e)
    st.header("16. UpSet Plot")
    try:
        st.image(Image.open(f"./file/{session_id}/image/plot2_7.png"))
    except:
        st.error("No significant terms were enriched")

# 執行 dose3_plot2_8.r， 畫DOSE部分的第8張圖
@st.cache_data
def draw_DOSE_plot2_8(exp_design, filer_missing_values, normalizeOption_py, control_py, alpha_py, lfc_py, ratioName_py, pvalue_2_8_py):
    os.system(f"Rscript dose3_plot2_8.r {session_id} {ratioName_py} {pvalue_2_8_py}")

# ---------------------------------------------------------設定網頁--------------------------------------------------------
# 設定網頁標題
st.title('Protein Analysis')
st.sidebar.header("0. Upload File")
st.header("1. Configure data for analysis (DEP)")
session_id = get_script_run_ctx().session_id
filename_py, data = upload_file()

# 刪除之前分析的資料與結果
delete_file(data)


# 分析資料要抓檔案中的 "Gene names"與 "Protein IDs"欄位，但不是每個檔都是叫一樣的名稱，
# 可能因為沒有空格、大小寫就抓不到，所以如果找不到欄位，就讓使用者自己選。
colname_geneNames_index = 0  if not "Gene.names"  in data.columns else data.columns.get_loc("Gene.names")
colname_proteinIDs_index = 0 if not "Protein.IDs" in data.columns else data.columns.get_loc("Protein.IDs")

colname_proteinIDs_py = st.sidebar.selectbox(label = "select 'Protein IDs'", options = data.columns, index = colname_proteinIDs_index)
colname_geneNames_py  = st.sidebar.selectbox(label = "select 'Gene names'", options = data.columns, index = colname_geneNames_index)

if colname_proteinIDs_py != "Protein.IDs" or colname_geneNames_py != "Gene.names":
    st.sidebar.warning("請確認'Protein IDs'與'Gene names'的欄位名稱選擇是否正確")

with st.expander("see data"):
    st.write(data)

st.sidebar.header("1. Configure data for analysis")
species_py, condition_py, numCondition_py = initialize_DEP_data(data)

select = False if numCondition_py > 1 else True
option = st.sidebar.selectbox(options=condition_py, label="Select condition", disabled=select)

data_exp_design_dict, data_log2 = select_condition_median_outlier(data, option) # colnames, medians, outliers
selected_col_id_list = draw_distogram_after_select_condition(data_log2, data_exp_design_dict)


# 控制 RUN、RERUN button
if 'CONFIG' not in st.session_state:
    st.session_state.CONFIG =  False

def change_configure_state(status):
    st.session_state.CONFIG  = status

# 確保當參數改變時，真的有重新畫圖 (因為可能會有 一開始參數為A，改成 B後重新畫B的圖，但之後又設回A，卻沒有重新畫A的圖)
def clear_cache():
    try:
        st.cache_data.clear()
    except:
        pass
def clear_cache_DEP_data2():
    try:
        cache_DEP_data2.clear()
    except:
        pass
def clear_cache_draw_dose_pic():
    try:
        draw_DOSE_pic.clear()
    except:
        pass

def config_data():
    # 製作 exp_design.csv
    exp_design = generate_exp_design_inputCondition(selected_col_id_list, data_exp_design_dict['colnames'])
    exp_design_condition_py = sorted(set(exp_design['condition']))
    # 選擇control、contrast
    control_py =  st.sidebar.selectbox(options=exp_design_condition_py, on_change=clear_cache_DEP_data2, label="Control:")

    contrastOption_py = list(exp_design_condition_py) #Str Object 轉為list
    contrastOption_py.remove(control_py) #移除control
    contrast_py =  st.sidebar.selectbox(options=contrastOption_py, on_change=clear_cache_DEP_data2, label="Contrast:")

    control_py = str(control_py)
    contrast_py = str(contrast_py)
    # 選擇 ratioColname
    ratioColname_py = [ f"{contrast}_vs_{control_py}_ratio" for contrast in contrastOption_py ]
    ratioName_py = st.sidebar.selectbox(label= "Ratio: ", on_change=clear_cache_draw_dose_pic, options= ratioColname_py)
    # 選擇物種
    species_list = ["mouse", "rat", "human"]
    species_index = [index for index, species in enumerate(species_list) if species == species_py]
    species_py_new = st.sidebar.selectbox(options=species_list, on_change=clear_cache_draw_dose_pic, index=species_index[0], label="Species:")
    # 選擇正規化方式
    normalizeOption_py = st.sidebar.selectbox(options=["Log2", "Median", "Mean", "VSN", "Quantile", "Cyclic Loess", "RLR", "Global Intensity"], on_change=clear_cache, label="Normalize")

    st.sidebar.subheader("1-2 Filter on missing values:")
    # Filter for proteins that are identified in all replicates of at least one condition
    # 選擇 filer的方式、參數
    maxReplicate_py = max(exp_design['replicate'])
    filer_missing_values = {'option': "condition", "nThr": 0, "min": 0.0}

    filer_missing_values['option'] = st.sidebar.selectbox(options=["condition", "complete", "fraction"], on_change=clear_cache, label="Filer options")
    if (filer_missing_values['option']  == "condition"):
        filer_missing_values['nThr']  = st.sidebar.slider(' the threshold for the allowed number of missing values in at least one condition:',
                                          min_value = 0,max_value = maxReplicate_py, on_change=clear_cache ,value = 1, step=1, format="%d")
    elif (filer_missing_values['option']  == "fraction"):
        filer_missing_values['min']  = st.sidebar.number_input("the threshold for the minimum fraction of valid values allowed for any protein:",
                                               min_value=0.0, max_value=1.0, on_change=clear_cache, value=0.66, step=0.1)

    st.sidebar.subheader("1-5. Differential enrichment analysis")
    alpha_py = st.sidebar.number_input("alpha:", on_change=clear_cache_DEP_data2, min_value = 0.0,max_value = 1.0 ,value = 1.0, step=0.01)
    lfc_py = st.sidebar.number_input("lfc = log2(value):", on_change=clear_cache_DEP_data2, min_value=0.0, max_value=2.0, value=1.5, step=0.1)
    # 待使用者選定參數並按下RUN按鈕後再開始畫圖
    with st.sidebar:
        with st.form(key="config_data_form"):
            c1, c2 = st.columns(2)
            with c1:
                if st.form_submit_button("🟢 DONE", on_click=change_configure_state, args=(True,), disabled=st.session_state.CONFIG):
                    pass
            with c2:
                if st.form_submit_button("🔴 RERUN", on_click=change_configure_state, args=(False,), disabled=not st.session_state.CONFIG):
                    pass

    if st.session_state.CONFIG:

        st.header('2. Filter on missing values:')
        # plot 1_1~1_5
        cache_DEP_data1(filename_py, session_id, colname_geneNames_py, colname_proteinIDs_py, filer_missing_values, normalizeOption_py)

        st.header("5. Differential enrichment analysis")
        # 1_6~1_8, 1_10
        cache_DEP_data2(filename_py, session_id, colname_geneNames_py, colname_proteinIDs_py, filer_missing_values, normalizeOption_py, control_py, alpha_py, lfc_py, contrast_py)

        st.header("7. Barplots of a protein of interest: ")
        st.sidebar.subheader("1-7. Barplots of a protein of interest")
        # 讀data_unique.csv(由dep1_to_1_5.r產生)，讓使用者選要看的基因名稱
        data_unique = pd.read_csv(f"./file/{session_id}/data_unique.csv", sep = ",", usecols=[colname_geneNames_py], error_bad_lines=False, dtype=object)

        dataGeneName = data_unique[colname_geneNames_py].dropna() #robjects.r("data_unique[, colname_geneNames]")
        dataGeneName = list(dataGeneName)
        dataGeneName = [x for x in dataGeneName if x != '']

        num_protein = st.sidebar.slider('Barplots of a protein of interest: ', min_value = 1, max_value = 20 ,value = 1, step=1, format="%d")
        proteinData_py = []

        for i in range(0, num_protein):
            proteinData_py.append( st.sidebar.selectbox(options=dataGeneName, index=i, label="protein", key=f"protein{i}") )

        cache_DEP_data3(filename_py, session_id, colname_geneNames_py, colname_proteinIDs_py, filer_missing_values, normalizeOption_py, control_py, alpha_py, lfc_py, num_protein, proteinData_py)


        st.header("8. Frequency plot of significant proteins and overlap of conditions")
        try :
            st.image(Image.open(f'./file/{session_id}/image/plot1_10.png'))
            st.image(Image.open(f"./file/{session_id}/image/plot1_10_2.png"))
        except:
            st.error('condition要有三組以上', icon="🚨")

        # DOSE
        DOSE_data_config(exp_design, filer_missing_values, normalizeOption_py, control_py, alpha_py, lfc_py, session_id, species_py_new, ratioName_py)

        #------------------------------------- 2. DOSE -------------------------------------
        st.sidebar.subheader("2. DOSE")
        #選擇 enrichment analysis方式
        enrichment_analysis_methods_py = st.sidebar.selectbox(label = "select enrichment_analysis: ",on_change= clear_cache_draw_dose_pic, options = ["DGN", "KEGG", "WikiPathways", "Reactome","Disease", "Universal"])
        universal_enrichment_category_py = 0
        universal_enrichment_subcategory_py = 0

        if (enrichment_analysis_methods_py == "Universal"):
            msigdbr_collections = pd.read_csv(f"./file/{session_id}/msigdbr_collections.csv", sep=",", dtype=object).fillna("")

            universal_enrichment_category_py = st.sidebar.selectbox(label = "Category: ", on_change= clear_cache_draw_dose_pic, options = list(dict.fromkeys(msigdbr_collections["gs_cat"])))
            subcategory_options = msigdbr_collections[msigdbr_collections["gs_cat"] == universal_enrichment_category_py]["gs_subcat"]
            subcategory_select_disable = False
            if len(subcategory_options) == 1:
                subcategory_select_disable = True
            universal_enrichment_subcategory_py = st.sidebar.selectbox(label = "Subcategory: ", disabled = subcategory_select_disable,on_change= clear_cache_draw_dose_pic, options = subcategory_options)

        de_up_down_py = st.sidebar.selectbox(options=['de','up', 'down'], on_change= clear_cache_draw_dose_pic, index=0, label="geneList dataset:")
        range_py = st.sidebar.number_input('range: ', on_change= clear_cache_draw_dose_pic, min_value = -20.0, max_value = 20.0, value = 1.5, step=0.01)
        if de_up_down_py == "down" and range_py > 0:
            st.sidebar.warning("the value of fold change should be less than zero.")

        st.sidebar.subheader("16. UpSet Plot pvalue")
        pvalue_2_8_py = st.sidebar.number_input("pvalue: ", min_value=0.001, max_value=1.0, step=0.001, value=0.05)
        draw_DOSE_pic(exp_design, filer_missing_values, normalizeOption_py, control_py, alpha_py, lfc_py, ratioName_py, de_up_down_py, range_py, enrichment_analysis_methods_py, universal_enrichment_category_py, universal_enrichment_subcategory_py)
        draw_DOSE_plot2_8(exp_design, filer_missing_values, normalizeOption_py, control_py, alpha_py, lfc_py, ratioName_py, pvalue_2_8_py)
        try:
            st.image(Image.open(f"./file/{session_id}/image/plot2_8.png"))
        except:
            st.warning("no term enriched under specific pvalueCutoff")
        st.header("17. ridgeline plot for expression distribution of GSEA result")
        try:
            st.image(Image.open(f"./file/{session_id}/image/plot2_9.png"))
        except Exception as e:
            st.error(e)

        st.header("18. running score and preranked list of GSEA result")
        for i in range(1, 11):
            try:
                st.image(Image.open(f"./file/{session_id}/image/plot2_10_{i}.png"))
            except Exception as e:
                st.error(e)

        with st.sidebar:
            st.subheader("19. Download result file")
            download_button(f"./file/{session_id}/dep_output_result.csv", "Download dep_output_result.csv")

        st.success('DONE!', icon="✅")
    # RERUN把cache清掉
    else:
        print("CLEAR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        st.cache_data.clear()
        # st.write(st.session_state)

config_data()