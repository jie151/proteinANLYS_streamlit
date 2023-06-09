import os

os.environ['R_HOME']= 'C:\\Program Files\\R\\R-4.2.3'
os.environ['PATH'] += os.pathsep + 'C:\\Program Files\\R\\R-4.2.3\\bin\\X64\\'
os.environ['PATH'] += os.pathsep + 'C:\\Program Files\\R\\R-4.2.3\\'

from rpy2 import robjects
import streamlit as st
import pandas as pd
import os
from PIL import Image
from contextlib import contextmanager, redirect_stdout
from io import StringIO
import streamlit.components.v1 as components
import base64
import uuid
import re

# Load package in r
robjects.r('''
    library(cowplot) # save_plot
    library(dplyr)
    library(ggplot2)
    library(httr)
    library(clusterProfiler)
    library(DEP)
    library(DOSE)
    library(enrichplot)
    library(NormalyzerDE)
    library(SummarizedExperiment)
    library(biomaRt)
''')

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


# 儲存檔案
def save_uploadedfile(uploadedfile):
    with open(os.path.join("uploadFile",uploadedfile.name),"wb") as f:
        f.write(uploadedfile.getbuffer())

# 上傳檔案
def upload_file():
    # 上傳檔案
    uploaded_file = st.sidebar.file_uploader('Upload a text file(.csv or .txt)', type=['csv', 'txt'], accept_multiple_files=False)
    if uploaded_file is not None:
        # 將上傳的檔案儲存下來 (由R開啟)
        filename = "./file/uploadFile/uploadFile.txt"
        with open(filename,"wb") as f:
            f.write( uploaded_file.getbuffer())

        sep_word = "," if uploaded_file.type == "text/csv" else " "
        data = pd.read_csv(uploaded_file, sep = sep_word, error_bad_lines=False)
    else:
        filename = "./file/uploadFile/proteinGroups_HsinYuan_Rat.txt"
        data = pd.read_csv(filename, sep=" ")
    return filename, data

def r_initialize_data(filename_py):
    robjects.r.assign('filename', filename_py)
    robjects.r('''
        #判斷物種,在DOSE時會使用
        MusMusculus_OX <- "OX=10090" #mouse
        RattusNorvegicus_OX <- "OX=10116" #rat

        data <- read.csv(file=filename, header=TRUE, fileEncoding ="UTF-8", sep = '\t')

        cat("file's row * column =", dim(data), "\n")

        if ( "Reverse" %in% names(data) ){ #如果有的話，要做篩選
            data <- filter(data, Reverse != "+")
        }
        if ( "Potential.contaminant" %in% names(data) ){
            data <- filter(data, Potential.contaminant != "+")
        }
        # We filter for contaminant proteins and decoy database hits, which are indicated by "+" in the columns "Potential.contaminants" and "Reverse", respectively.

        if ( grepl(MusMusculus_OX, data[1,"Fasta.headers"], fixed = TRUE) ) {
            species <- "mouse"
        }else if ( grepl(RattusNorvegicus_OX, data[1,"Fasta.headers"], fixed = TRUE) ) {
            species <- "rat"
        }else {
            species <- "human"
        }

        df <- data[ , grepl( "Reporter.intensity.corrected" , names( data) )]
        conditionColumn <- colnames(df)
        condition <- gsub("Reporter.intensity.corrected.[0-9]+", "", conditionColumn)
        condition <- unique(condition) #去掉重複的condition
        numCondition <- length(condition)
    ''')
    numCondition_py = robjects.r("numCondition")
    condition_py = robjects.r("condition")
    return numCondition_py, condition_py

def r_select_condition(numCondition_py, condition_py):
    select = False if numCondition_py[0] > 1 else True
    option = st.sidebar.selectbox(options=condition_py, label="Select condition", disabled=select)

    robjects.r.assign('conditionName', option)
    robjects.r('''
        df <- data[ , grepl( "Reporter.intensity.corrected" , names( data) )]
        df <- df[, grepl(conditionName, names(df))]
        df = log2(df)
        colname <- colnames(df)
        numColname <- length(colname)
        # 找離群值
        median_list <- apply( df[, 1: numColname ], 2, median)
        outliers_val <-  boxplot(median_list)$out

        draw_hist <- function(id) { #1~10, 例如Reporter.intensity.corrected.1-10.LA
            return (hist(x=df[,id], breaks=25,
                    xlim=c(0,max(df)), main=colname[id], # 圖片的名稱
                    xlab="", ylab="" ))
        }
    ''')
    df_colname_py = robjects.r("colname") #python的變數 df_colname_py = colname
    df_ncol_py = robjects.r("numColname")
    median_list_py = robjects.r("median_list")

    return df_colname_py, df_ncol_py, median_list_py

def generate_default_group(num_selected_col):
    condition_replicate_list = []
    num_category = 1

    for i in range(2, num_selected_col):
        if num_selected_col % i == 0:
            num_category = i
            break
    a_category_colnum = int(num_selected_col / num_category)
    for i in range(0, num_category):
        for j in range(1, a_category_colnum+1):
            # 從A開始編號
            condition_replicate_list.append([chr(65+i), j])
    return condition_replicate_list

def r_experimental_design_file(experimental_design):
    experimental_design = experimental_design[(experimental_design.state != "N")]
    experimental_design = experimental_design.drop("state", axis=1)
    experimental_design.to_csv("./file/experimental_design.csv", index=False)
    st.write("experimental_design: ",experimental_design)

    with st.expander("data Info. "):
        output = st.empty()

        with st_capture(output.code):
            robjects.r('''
                # Are there any duplicated gene names?
                cat('* Are there any duplicated gene names? ', data$Gene.names %>% duplicated() %>% any(), "\n")

                if ( data$Gene.names %>% duplicated() %>% any() ){
                    # Make a table of duplicated gene names
                    write.table(data %>% group_by(Gene.names) %>% summarize(frequency = n()) %>%
                        arrange(desc(frequency)) %>% filter(frequency > 1), file = "./file/my_data1.txt", row.names =FALSE)
                }
                cat("a table of duplicated gene names: (table ", dim(table), "\n")
                table <- read.csv(file="./file/my_data1.txt", header=TRUE, fileEncoding ="UTF-8", sep = ' ')
                print(head(table, 7))

                # Make unique names using the annotation in the "Gene.names" column as primary names and the annotation in "Protein.IDs" as name for those that do not have an gene name.
                data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")

                # Generate a SummarizedExperiment object using an experimental design
                LFQ_columns <- grep("Reporter.intensity.corrected.", colnames(data_unique)) # get LFQ column numbers

                experimental_design <- read.csv('./file/experimental_design.csv',header=TRUE ,fileEncoding ="UTF-8")
                experimental_design$label = gsub(" ", ".", experimental_design$label)

                maxReplicate <- max(experimental_design$replicate) # max replicate
                experimental_design_condition <- unique(experimental_design$condition)

                data_se <- make_se(data_unique, LFQ_columns, experimental_design)

                # Let's have a look at the SummarizedExperiment object
                cat("* the SummarizedExperiment object:  \n")
                print(data_se)
            ''')

def draw_distogram_after_select_condition(df_colname_py, df_ncol_py, median_list_py):
    with st.expander("Select the set of columns you want to analyze"):
        # 把圖分成有選擇與沒選擇分別顯示
        cols = st.columns(3)
        outliers_val_py = robjects.r("outliers_val")
        # 有選擇的欄位編號
        selected_col_id_list = []
        for i in range(0, df_ncol_py[0]):
            check_value = True if median_list_py[i] not in outliers_val_py else False

            colname_checkobx = cols[ i % 3 ].checkbox(df_colname_py[i], value=check_value)
            if colname_checkobx :
                selected_col_id_list.append(i)

            robjects.r.assign("id", df_colname_py[i]) #將檔案名稱python -> R
            robjects.r('''
                    png(file="./file/image/histogram.png", width=250, height=250)
                    draw_hist(id)
                    dev.off()
                ''')
            cols[i % 3].text(f"median: {round(median_list_py[i], 3)}")
            cols[ i % 3 ].image(Image.open('./file/image/histogram.png'))
        return selected_col_id_list

def generate_experimental_design_inputCondition(selected_col_id_list):
    experimental_design = pd.DataFrame(columns=["state","label","condition","replicate"], index=range(1, 11)).fillna("N") #創一個dataframe, 預設填N
    # 根據欄位數，預設分組與輸入值
    condition_replicate_list = generate_default_group(len(selected_col_id_list))

    with st.sidebar.expander("Select property columns"):
        for index, col_id in enumerate(selected_col_id_list):
            experimental_design.at[col_id + 1, 'state'] = "T"
            experimental_design.at[col_id + 1, 'label'] = df_colname_py[col_id]
            st.write(df_colname_py[col_id])
            condition_textInput = st.text_input('condition', key= f"condition{col_id}", value=condition_replicate_list[index][0])
            replicate_textInput = st.text_input('replicate', key= f"replicate{col_id}", value=condition_replicate_list[index][1])
            if condition_textInput:
                experimental_design.at[col_id + 1, 'condition'] = condition_textInput
            if replicate_textInput:
                experimental_design.at[col_id + 1, 'replicate'] = replicate_textInput

    return experimental_design

def r_normalize_function():
    robjects.r('''
        logNorm <- function(dat) {
            logInt <- log2(dat)
            #logInt <- replace(is.infinite(logInt), NA)
            logInt[is.infinite(as.matrix(logInt))] <- NA
            return(as.matrix(logInt))
        }

        medianNorm <- function(logDat) {
            # Find medians of each sample
            # Divide by median
            # Multiply by mean of medians
            sampleMed <- apply(logDat, 2, median, na.rm=TRUE)
            meanMed <- mean(sampleMed, na.rm=TRUE)
            out <- t(t(logDat) / sampleMed)
            out <- out * meanMed
            return(as.matrix(out))
        }

        meanNorm <- function(logDat) {
            # Find means of each sample
            # Divide by mean
            # Multiply by mean of means
            sampleMean <- apply(logDat, 2, mean, na.rm=TRUE)
            meanMean <- mean(sampleMean, na.rm=TRUE)
            out <- t(t(logDat) / sampleMean)
            out <- out * meanMean
            return(as.matrix(out))
        }

        vsnNorm <- function(dat) {
            vsnNormed <- suppressMessages(vsn::justvsn(as.matrix(dat)))
            colnames(vsnNormed) <- colnames(dat)
            row.names(vsnNormed) <- rownames(dat)
            return(as.matrix(vsnNormed))
        }

        quantNorm <- function(logDat) {
            quantNormed <- preprocessCore::normalize.quantiles(as.matrix(logDat), copy=FALSE)
            colnames(quantNormed) <- colnames(logDat)
            row.names(quantNormed) <- rownames(logDat)
            return(as.matrix(quantNormed))
        }

        cycLoessNorm <- function(logDat) {
            cycLoessNormed <- limma::normalizeCyclicLoess(as.matrix(logDat), method="fast")
            colnames(cycLoessNormed) <- colnames(logDat)
            row.names(cycLoessNormed) <- rownames(logDat)
            return(as.matrix(cycLoessNormed))
        }

        rlrNorm <- function(logDat) {
            rlrNormed <- NormalyzerDE::performGlobalRLRNormalization(as.matrix(logDat), noLogTransform=TRUE)
            colnames(rlrNormed) <- colnames(logDat)
            row.names(rlrNormed) <- rownames(logDat)
            return(as.matrix(rlrNormed))
        }

        giNorm <- function(logDat) {
            giNormed <- NormalyzerDE::globalIntensityNormalization(as.matrix(logDat), noLogTransform=TRUE)
            colnames(giNormed) <- colnames(logDat)
            row.names(giNormed) <- rownames(logDat)
            return(as.matrix(giNormed))
        }

        normalize_proteiNorm <- function(se, normalizeOption) {
            # Show error if inputs are not the required classes
            assertthat::assert_that(inherits(se, "SummarizedExperiment"))

            # Variance stabilization transformation on assay data
            se_vsn <- se

            if(normalizeOption == "VSN")

                assay(se_vsn)  <- vsnNorm(2 ^ assay(se_vsn))

            else if(normalizeOption == "Log2")

                assay(se_vsn)  <- logNorm(2 ^ assay(se_vsn))

            else if(normalizeOption == "Median")

                assay(se_vsn)  <- medianNorm(logNorm(2 ^ assay(se_vsn)))

            else if(normalizeOption == "Mean")

                assay(se_vsn)  <- meanNorm(logNorm(2 ^ assay(se_vsn)))

            else if(normalizeOption == "Quantile")

                assay(se_vsn)  <- quantNorm(logNorm(2 ^ assay(se_vsn)))

            else if(normalizeOption == "Cyclic Loess")

                assay(se_vsn)  <- cycLoessNorm(logNorm(2 ^ assay(se_vsn)))

            else if(normalizeOption == "RLR")

                assay(se_vsn) <- rlrNorm(logNorm(2 ^ assay(se_vsn)))

            else
                assay(se_vsn)  <- giNorm(logNorm(2 ^ assay(se_vsn)))

            return(se_vsn)
            }
    ''')


# Plot a barplot of the protein identification overlap between samples
@st.cache_data
def r_plot_frequency_1_1():
    robjects.r(''' save_plot("./file/image/plot1_1.png", plot =  plot_frequency(data_se), , base_height = 4, base_width = 4.5) ''')
    st.write(" *Plot a barplot of the protein identification overlap between samples")

    st.image(Image.open('./file/image/plot1_1.png'))
    with st.expander("data"):
        output = st.empty()
        with st_capture(output.code):
            robjects.r(''' print(plot_frequency(data_se, plot= FALSE)) ''')

# Plot a barplot of the number of identified proteins per samples
@st.cache_data
def r_plot_numbers_filter_missval_1_2(nThr_py):
    robjects.r.assign("nThr", nThr_py)
    robjects.r('''
        data_filt <- filter_missval(data_se, thr = nThr) #讓使用者選0~4(重複)
        save_plot("./file/image/plot1_2.png", plot =  plot_numbers(data_filt), base_height = 4, base_width = 4.5)
    ''')
    st.write(" *Plot a barplot of the number of identified proteins per samples")
    st.image(Image.open('./file/image/plot1_2.png'))

# Plot a barplot of the protein identification overlap between samples
@st.cache_data
def r_plot_coverage_1_3():
    robjects.r(''' save_plot("./file/image/plot1_3.png", plot =  plot_coverage(data_filt), base_height = 4, base_width = 4) ''')
    st.write(" *Plot a barplot of the protein identification overlap between samples")
    st.image(Image.open('./file/image/plot1_3.png'))

@st.cache_data
def r_plot_normalization_1_4():

    robjects.r('''
        data_norm <- normalize_proteiNorm(data_filt, normalizeOption)
        pic1 <- plot_normalization(data_filt, data_norm)
        save_plot("./file/image/plot1_4.png", plot =  pic1, base_height = 4, base_width = 4)
    ''')

    st.write(" Visualize normalization by boxplots for all samples before and after normalization")
    st.image(Image.open('./file/image/plot1_4.png'))

# Plot a heatmap of proteins with missing values
@st.cache_data
def r_plot_heatmap_1_5():
    robjects.r('''
        png(file = "./file/image/plot1_5_1.png")
        pic1 <- plot_missval(data_filt)
        dev.off()

        # Plot intensity distributions and cumulative fraction of proteins with and without missing values
        png(file = "./file/image/plot1_5_2.png")
        pic1 <- plot_detect(data_filt)
        dev.off()

        # All possible imputation methods are printed in an error, if an invalid function name is given.
        #impute(data_norm, fun = "")

        # Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
        data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)

        # Impute missing data using random draws from a manually defined left-shifted Gaussian distribution (for MNAR)
        #data_imp_man <- impute(data_norm, fun = "man", shift = 1.8, scale = 0.3)

        # Impute missing data using the k-nearest neighbour approach (for MAR)
        #data_imp_knn <- impute(data_norm, fun = "knn", rowmax = 0.9)

        # Plot intensity distributions before and after imputation
        pic1 <- plot_imputation(data_norm, data_imp)
        save_plot("./file/image/plot1_5_3.png", plot =  pic1, base_width = 4, base_height = 4)
    ''')
    image_label_list = ["*Plot a heatmap of proteins with missing values",
                        "*Plot intensity distributions and cumulative fraction of proteins with and without missing values",
                        "*Plot intensity distributions before and after imputation"]
    for i in range(1, 4):
        st.write(f"{image_label_list[i-1]}")
        st.image(Image.open(f'./file/image/plot1_5_{i}.png'))

@st.cache_data
def r_plot_pca_1_6(control_py, alpha_py, lfc_py):
    robjects.r.assign("control_r", control_py)
    robjects.r.assign("alpha_r", alpha_py)
    robjects.r.assign("lfc_r", lfc_py)
    robjects.r('''
        # Test every sample versus control
        data_diff <- test_diff(data_imp, type = "control", control = control_r)
        # Denote significant proteins based on user defined cutoffs
        dep <- add_rejections(data_diff, alpha = alpha_r, lfc = log2(lfc_r)) #alpha、lfc調整
        pic1 <- plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4)
        save_plot("./file/image/plot1_6_1.png", pic1)

        # Plot the Pearson correlation matrix
        png(file="./file/image/plot1_6_2.png")
        pic1 <- plot_cor(dep, significant = FALSE, lower = 0.9, upper = 1, pal = "Reds")
        dev.off()
    ''')
    image_label_list = ["*Plot the first and second principal components",
                        "*Plot the Pearson correlation matrix"]
    for i in range(1, 3):
        st.write(image_label_list[i-1])
        st.image(Image.open(f'./file/image/plot1_6_{i}.png'))

@st.cache_data
def r_plot_heatmap_dep_1_7():
    robjects.r('''
        # Plot a heatmap of all significant proteins with the data centered per protein
        png(file="./file/image/plot1_7_1.png")
        pic1 <- plot_heatmap(dep, type = "centered", kmeans = TRUE,
                    k = 6, col_limit = 4, show_row_names = FALSE,
                    indicate = c("condition", "replicate"))
        dev.off()

        # Plot a heatmap of all significant proteins (rows) and the tested contrasts (columns)
        png(file="./file/image/plot1_7_2.png")
        pic1 <- plot_heatmap(dep, type = "contrast", kmeans = TRUE,
                    k = 6, col_limit = 10, show_row_names = FALSE)
        dev.off()
    ''')
    image_label_list = ["*Plot a heatmap of all significant proteins with the data centered per protein",
                        "*Plot a heatmap of all significant proteins (rows) and the tested contrasts (columns)"]
    for i in range(1, 3):
        st.write(image_label_list[i-1])
        st.image(Image.open(f'./file/image/plot1_7_{i}.png'))

@st.cache_data
def r_plot_volcano_1_8(contrast_py):
    robjects.r.assign("contrast_r", contrast_py)
    robjects.r('''
        contrastSample <- paste(contrast_r, "_vs_", control_r)
        contrastSample <- gsub(" ", "", contrastSample, fixed = TRUE)
        print(contrastSample)

        # Plot a volcano plot for the contrast "Ubi6 vs Ctrl""
        pic <- plot_volcano(dep, contrast = contrastSample, label_size = 2, add_names = TRUE)
        save_plot("./file/image/plot1_8.png", pic)
    ''')
    st.image(Image.open('./file/image/plot1_8.png'))

def r_plot_single_1_9():
    dataGeneName = robjects.r("data_unique$Gene.names")
    dataGeneName = list(dataGeneName)
    dataGeneName = [x for x in dataGeneName if x != '']

    num_protein = st.sidebar.slider('Barplots of a protein of interest: ', min_value = 1, max_value = 20 ,value = 1, step=1, format="%d")
    protein_type = "contrast" if num_protein > 1 else "centered"

    proteinData_py = []
    for i in range(0, num_protein):
        proteinData_py.append( st.sidebar.selectbox(options=dataGeneName, index=i, label="protein", key=f"protein{i}") )

    # 確認選擇的蛋白質有無重複
    if len(proteinData_py) == len(set(proteinData_py)):
        robjects.r.assign("protein_type_r", protein_type)
        robjects.r.assign("proteinData", proteinData_py)
        robjects.r('''
            pic1 <- plot_single(dep, proteins = unlist( proteinData) , type = protein_type_r)
            save_plot("./file/image/plot1_9.png", pic1)
        ''')
        st.image(Image.open('./file/image/plot1_9.png'))
    else:
        st.error("Error!! The proteins selected cannot be duplicated", icon="🚨")

@st.cache_data
def r_plot_cond_1_10():
    robjects.r('''
        # Generate a results table
        data_results <- get_results(dep)
        # Number of significant proteins
        data_results %>% filter(significant) %>% nrow()
        write.csv(data_results,"./file/dep_output.csv", row.names = FALSE, quote=F)

        # Plot a frequency plot of significant proteins for the different conditions
        png(file="./file/image/plot1_10.png")
        pic1 <- plot_cond(dep)
        dev.off()

    ''')
    st.image(Image.open('./file/image/plot1_10.png'))

@st.cache_data
def r_uniprotAPI():
    robjects.r('''
        # uriprot.org API
        isJobReady <- function(jobId) {
            pollingInterval = 15
            nTries = 60
            for (i in 1:nTries) {
                url <- paste("https://rest.uniprot.org/idmapping/status/", jobId, sep = "")
                r <- GET(url = url, accept_json())
                status <- httr::content(r, as = "parsed")
                if (!is.null(status[["results"]]) || !is.null(status[["failedIds"]])) {
                    return(TRUE)
                }
                if (!is.null(status[["messages"]])) {
                    print(status[["messages"]])
                    return (FALSE)
                }
                    Sys.sleep(pollingInterval)
            }
                return(FALSE)
        }
        getResultsURL <- function(redirectURL) {
            if (grepl("/idmapping/results/", redirectURL, fixed = TRUE)) {
                url <- gsub("/idmapping/results/", "/idmapping/stream/", redirectURL)
            } else {
                url <- gsub("/results/", "/results/stream/", redirectURL)
            }
        }

        generate_uniprot_resultTable <- function(ID_string){
            result_1 <- POST("https://rest.uniprot.org/idmapping/run",
                        add_headers("user-agent" = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/110.0.0.0 Safari/537.36"),
                        body = list(ids =ID_string , 'from'="UniProtKB_AC-ID",'to'="GeneID")) #c("C0HK80,C0HKD9")
            submission <- httr::content(result_1, as = "parsed")
            url <- paste("https://rest.uniprot.org/idmapping/details/", submission[["jobId"]], sep = "")
            r <- GET(url = url, accept_json())
            details <- httr::content(r, as = "parsed")
            url <- getResultsURL(details[["redirectURL"]])

            # Using TSV format see: https://www.uniprot.org/help/api_queries#what-formats-are-available
            url <- paste(url, "?format=tsv", sep = "")
            r <- GET(url = url, accept_json())
            resultsTable = read.table(text = httr::content(r), sep = "\t", header=TRUE, col.names=c('ID_uniprot', 'GeneID_entrez'), encoding = "UTF-8")
            return (resultsTable)
        }
    ''')

@st.cache_data
def r_init_DOSE_data():

    robjects.r('''
        dep_output_data = read.csv(file = "./file/dep_output.csv", header=TRUE, fileEncoding ="UTF-8")
        dep_output_data[,'ID'] <- sub("-.*", "", dep_output_data$ID)

        ratioCol <- grep( "ratio$" , colnames(dep_output_data) )
        ratioColname <- colnames(dep_output_data[ratioCol])
        numRatio <- length(ratioCol)
    ''')

    try:
        robjects.r('''
            # convert uniprot -> entrez
            ID_string <- paste(dep_output_data$ID, collapse=",")
            resultsTable <- generate_uniprot_resultTable(ID_string)
        ''')
    except:
        st.error('與uriprot.org網站連接失敗(Entrez -> GeneID), 請重新載入網頁', icon="🚨")

    robjects.r('''
        uniprot_entrez <- resultsTable[-c(1),]
        #去除重複 (因為uniprot可能會有多個entrez => 排序後選最小的)
        uniprot_entrez[, 2] <- as.numeric(uniprot_entrez[, 2])
        uniprot_entrez <- uniprot_entrez [ order(uniprot_entrez$ID_uniprot, uniprot_entrez$GeneID_entrez),]
        uniprot_entrez <- uniprot_entrez[!duplicated(uniprot_entrez$ID_uniprot),]
        write.csv(uniprot_entrez, file="./file/uniprot_entrez.csv")
    ''')

@st.cache_data
def r_convert_species_gene(rationName_py):
    robjects.r.assign("ratioName", rationName_py)
    robjects.r('''
        if(species != "human"){
            human <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
            print("******not human******")
            if(species == "mouse"){
                mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
                human_mouse_rat <- getLDS(mart=mouse, attributes=c("entrezgene_id"), filters="entrezgene_id" , values=uniprot_entrez[,"GeneID_entrez"], attributesL= c("entrezgene_id"), martL = human, uniqueRows=T ) #values=uniprots_ID
                print("******species => mouse******")
            }else{
                rat <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
                human_mouse_rat <- getLDS(mart=rat, attributes=c("entrezgene_id"), filters="entrezgene_id" , values=uniprot_entrez[,"GeneID_entrez"], attributesL= c("entrezgene_id"), martL = human, uniqueRows=T ) #values=uniprots_ID
                print("******species => rat******")
            }
            colnames(human_mouse_rat) <- c('GeneID_entrez','human')
            human_mouse_rat <-human_mouse_rat[ order(human_mouse_rat$GeneID_entrez, human_mouse_rat$human),]
            human_mouse_rat <- human_mouse_rat[!duplicated(human_mouse_rat$GeneID_entrez),]

            temp <- inner_join( uniprot_entrez, human_mouse_rat, by="GeneID_entrez")
            colnames(temp) <- c('ID', 'GeneID', 'human_entrez')
            temp[, "ID"] <- as.character(temp[, "ID"])

            data <- inner_join(temp, dep_output_data, by="ID")
            print(dim(data))
            print(head(data))
        }else{
            print("species => human******")
            colnames(uniprot_entrez) <- c("ID", "human_entrez")
            data <- inner_join(uniprot_entrez, dep_output_data, by="ID")
        }
        write.csv(data, file="./file/dep_output_result.csv")

        ## feature 1: numeric vector
        geneList = data[,ratioName] #-log(pValue)

        ## feature 2: named vector
        names(geneList) = as.character(data[, 'human_entrez'])

        ## feature 3: decreasing orde
        geneList = sort(geneList, decreasing = TRUE)
        print(head(geneList))

        de <- names(geneList)[abs(log(geneList)) > 1]
        edo <- enrichDGN(de)
        edo2 <- gseDO(geneList, pvalueCutoff=1)
        edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')

    ''')

@st.cache_data
def r_test_file():
    robjects.r('''

        data(geneList)
        de <- names(geneList)[abs(geneList) > 2]
        edo <- enrichDGN(de)
        edo2 <- gseDO(geneList, pvalueCutoff=1)
        edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')

    ''')

@st.cache_data
def r_plot_barplot_2_1():
    robjects.r('''
        pic1 <- barplot(edo,x = "Count", color="p.adjust", showCategory=20)+ xlab("Count")
        save_plot("./file/image/plot2_1.png", pic1, base_height = 10, base_aspect_ratio = 1)
    ''')
    st.image(Image.open('./file/image/plot2_1.png'))

@st.cache_data
def r_plot_dotplot_2_2():
    robjects.r('''
        pic11_1 <- dotplot(edo, showCategory=30) + ggtitle("dotplot for ORA")
        pic11_2 <- dotplot(edo2, showCategory=30) + ggtitle("dotplot for GSEA")
        save_plot("./file/image/plot2_2_1.png", pic11_1, base_height = 10, base_aspect_ratio = 1)
        save_plot("./file/image/plot2_2_2.png", pic11_2, base_height = 10, base_aspect_ratio = 1)
    ''')
    st.image(Image.open('./file/image/plot2_2_1.png'))
    st.image(Image.open('./file/image/plot2_2_2.png'))

@st.cache_data
def r_plot_cnetplot_2_3():
    robjects.r('''
        rescale.AsIs <- function(x, ...){
        dropAsis <- function(x){
            cls <- class(x)
            structure(x, class = setdiff(cls, "AsIs"))
        }
        scales:::rescale(dropAsis(x), ...)
        }

        pic12_1 <- cnetplot(edox, foldChange=geneList)
        ## categorySize can be scaled by 'pvalue' or 'geneNum'
        pic12_2 <- cnetplot(edox, categorySize="pvalue", foldChange=geneList)
        pic12_3 <- cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
        pic12_4 <- cnetplot(edox, node_label="category")
        pic12_5 <- cnetplot(edox, node_label="gene")
        pic12_6 <- cnetplot(edox, node_label="all")
        pic12_7 <- cnetplot(edox, node_label="none")

        save_plot("./file/image/plot2_3_1.png", pic12_1, base_height = 10, base_aspect_ratio = 1)
        save_plot("./file/image/plot2_3_2.png", pic12_2, base_height = 10, base_aspect_ratio = 1)
        save_plot("./file/image/plot2_3_3.png", pic12_3, base_height = 10, base_aspect_ratio = 1)
        save_plot("./file/image/plot2_3_4.png", pic12_4, base_height = 10, base_aspect_ratio = 1)
        save_plot("./file/image/plot2_3_5.png", pic12_5, base_height = 10, base_aspect_ratio = 1)
        save_plot("./file/image/plot2_3_6.png", pic12_6, base_height = 10, base_aspect_ratio = 1)
        save_plot("./file/image/plot2_3_7.png", pic12_7, base_height = 10, base_aspect_ratio = 1)

    ''')
    for i in range(1, 8):
        st.image(Image.open(f"./file/image/plot2_3_{i}.png"))

@st.cache_data
def r_plot_heatplot_2_4():
    robjects.r('''
        pic13_1 <- heatplot(edox)
        pic13_2 <- heatplot(edox, foldChange=geneList)
        save_plot("./file/image/plot2_4_1.png", pic13_1)
        save_plot("./file/image/plot2_4_2.png", pic13_2)

        save_plot("./file/image/heatplot_1.png", pic13_1, base_height = 10, base_aspect_ratio = 3,limitsize = FALSE)
        save_plot("./file/image/heatplot_2.png", pic13_2, base_height = 10, base_aspect_ratio = 3,limitsize = FALSE)
    ''')
    for i in range(1, 3):
        st.image(Image.open(f'./file/image/plot2_4_{i}.png'))
        download_button(f"./file/image/heatplot_{i}.png", f"Download heatplot_{i}.png")

@st.cache_data
def r_plotenrichment_map_2_5():
    robjects.r('''
        edo <- pairwise_termsim(edo)
        pic14_1 <- emapplot(edo)
        pic14_2 <- emapplot(edo, cex_category=1.5)
        pic14_3 <- emapplot(edo, layout="kk")
        pic14_4 <- emapplot(edo, cex_category=1.5,layout="kk")
        save_plot("./file/image/plot2_5_1.png", pic14_1, base_height = 10, base_aspect_ratio = 1)
        save_plot("./file/image/plot2_5_2.png", pic14_2, base_height = 10, base_aspect_ratio = 1)
        save_plot("./file/image/plot2_5_3.png", pic14_3, base_height = 10, base_aspect_ratio = 1)
        save_plot("./file/image/plot2_5_4.png", pic14_4, base_height = 10, base_aspect_ratio = 1)
    ''')
    for i in range(1, 5):
        st.image(Image.open(f'./file/image/plot2_5_{i}.png'))

@st.cache_data
def r_plot_emapplot_2_6():
    output = st.empty()
    with st_capture(output.code):
        robjects.r('''
            xx <- compareCluster(data, fun="enrichKEGG",
                                organism="hsa", pvalueCutoff=0.05) #pvalueCutoff正常為0.05
            xx <- pairwise_termsim(xx)
            pic15_1 <- emapplot(xx)
            pic15_2 <- emapplot(xx, legend_n=2)
            pic15_3 <- emapplot(xx, pie="count")
            pic15_4 <- emapplot(xx, pie="count", cex_category=1.5, layout="kk")

            save_plot("./file/image/plot2_6_1.png", pic15_1, base_height = 10, base_aspect_ratio = 1)
            save_plot("./file/image/plot2_6_2.png", pic15_2, base_height = 10, base_aspect_ratio = 1)
            save_plot("./file/image/plot2_6_3.png", pic15_3, base_height = 10, base_aspect_ratio = 1)
            save_plot("./file/image/plot2_6_4.png", pic15_4, base_height = 10, base_aspect_ratio = 1)
    ''')
    for i in range(1, 5):
        st.image(Image.open(f'./file/image/plot2_6_{i}.png'))

@st.cache_data
def r_plot_upseplot_2_7():
    robjects.r('''
        save_plot("./file/image/plot2_7.png", upsetplot(edo), base_height = 10, base_aspect_ratio = 1.5)
    ''')
    st.image(Image.open("./file/image/plot2_7.png"))

@st.cache_data
def r_plot_upsetplot_with_splider_2_8():
    robjects.r('''
        kk2 <- gseKEGG(geneList = geneList,
                       organism = 'hsa',
                        nPerm = 1000,
                       minGSSize = 120,
                       pvalueCutoff = pvalue_2_8,
                       verbose = FALSE)
        if (length(kk2[,2]) < 1) {
            print("no term enriched under specific pvalueCutoff")
        }else {
            pic16_2 <- upsetplot(kk2)
            save_plot("./file/image/plot2_8.png", pic16_2, base_height = 10, base_aspect_ratio = 1.5)
        }
    ''')
    st.image(Image.open("./file/image/plot2_8.png"))

@st.cache_data
def r_plot_ridgeplot_2_9():
    robjects.r('''
        pic17 <- ridgeplot(edo2) + xlab("expression distributions of enriched genes (log2FC)")
        save_plot("./file/image/plot2_9.png", pic17, base_height = 12, base_aspect_ratio = 0.65)
    ''')
    st.image(Image.open("./file/image/plot2_9.png"))

@st.cache_data
def r_plot_gseaplot_2_10():
    robjects.r('''

        pic18_1 <- gseaplot(edo2, geneSetID = 1, by = "runningScore", title = edo2$Description[1])
        pic18_2 <- gseaplot(edo2, geneSetID = 1, by = "preranked", title = edo2$Description[1])
        pic18_3 <- gseaplot(edo2, geneSetID = 1, title = edo2$Description[1])

        #cowplot::plot_grid(pic9_1, pic9_2, pic9_3, ncol=1, labels=LETTERS[1:3])
        pic18_4 <- gseaplot2(edo2, geneSetID = 1, title = edo2$Description[1])

        pic18_5 <- gseaplot2(edo2, geneSetID = 1:3)
        pic18_6 <- gseaplot2(edo2, geneSetID = 1:3, pvalue_table = TRUE,
                color = c("#E495A5", "#86B875", "#7DB0DD"), ES_geom = "dot")

        pic18_7 <- gseaplot2(edo2, geneSetID = 1:3, subplots = 1) + xlab("Rank")
        pic18_8 <- gseaplot2(edo2, geneSetID = 1:3, subplots = 1:2) + xlab("Rank")


        pic18_9 <- gsearank(edo2, 1, title = edo2[1, "Description"])

        png("./test.png")
        pic18_10 <- lapply(1:3, function(i){
            anno <- edo2[i, c("NES", "pvalue", "p.adjust")]
            lab <- paste0(names(anno), "=",  round(anno, 3), collapse="\n")
            gsearank(edo2, i, edo2[i, 2]) + xlab("Rank")+
            annotate("text", 5000, edo2[i, "enrichmentScore"]* .75, label =lab, hjust=0, vjust=0)})
        dev.off()
        save_plot("./file/image/plot2_10_1.png", pic18_1, base_height = 10, base_aspect_ratio = 1)
        save_plot("./file/image/plot2_10_2.png", pic18_2, base_height = 10, base_aspect_ratio = 1)
        save_plot("./file/image/plot2_10_3.png", pic18_3, base_height = 10, base_aspect_ratio = 1)
        save_plot("./file/image/plot2_10_4.png", pic18_4, base_height = 10, base_aspect_ratio = 1)
        save_plot("./file/image/plot2_10_5.png", pic18_5, base_height = 10, base_aspect_ratio = 1)
        save_plot("./file/image/plot2_10_6.png", pic18_6, base_height = 10, base_aspect_ratio = 1)
        save_plot("./file/image/plot2_10_7.png", pic18_7, base_height = 10, base_aspect_ratio = 1)
        save_plot("./file/image/plot2_10_8.png", pic18_8, base_height = 10, base_aspect_ratio = 1)
        save_plot("./file/image/plot2_10_9.png", pic18_9, base_height = 10, base_aspect_ratio = 1)
        #save_plot("./file/image/plot2_10_10.png", pic18_10, base_height = 10, base_aspect_ratio = 1)
    ''')
    for i in range(1, 11):
        st.image(Image.open(f"./file/image/plot2_10_{i}.png"))


# 設定網頁標題
st.title('Protein Analysis')
st.header("1. Configure data for analysis (DEP)")
st.sidebar.header("0. Upload File")
filename_py, data = upload_file()
with st.expander("see data"):
    st.write(data)


css = r'''
   <style>
        [data-testid="stForm"] {border: 0px}
    </style>
'''
st.markdown(css, unsafe_allow_html=True)


st.sidebar.header("1. Configure data for analysis")
numCondition_py, condition_py = r_initialize_data(filename_py)
df_colname_py, df_ncol_py, median_list_py = r_select_condition(numCondition_py, condition_py)
selected_col_id_list = draw_distogram_after_select_condition(df_colname_py, df_ncol_py, median_list_py)

def all():
#with st.form("Configure data for analysis"):

    experimental_design = generate_experimental_design_inputCondition(selected_col_id_list)
    r_experimental_design_file(experimental_design)
    experimental_design_condition_py = robjects.r("experimental_design_condition")
    control_py =  st.sidebar.selectbox(options=experimental_design_condition_py, label="Control:")
    contrastOption_py = list(experimental_design_condition_py) #Str Object 轉為list
    contrastOption_py.remove(control_py[0]) #移除control
    contrast_py =  st.sidebar.selectbox(options=contrastOption_py, label="Contrast:")

    ratioColname_py = [ f"{contrast}_vs_{control_py}_ratio" for contrast in contrastOption_py ]
    rationName_py = st.sidebar.selectbox(label= "Ratio: ", options= ratioColname_py)

    species_py = robjects.r("species")
    species_list = ["mouse", "rat", "human"]
    species_index = [index for index, species in enumerate(species_list) if species == species_py[0]]
    species_py_new = st.sidebar.selectbox(options=species_list, index=species_index[0], label="Species:")
    robjects.r.assign("species", species_py_new)

    normalizeOption_py = st.sidebar.selectbox(options=["Log2", "Median", "Mean", "VSN", "Quantile", "Cyclic Loess", "RLR", "Global Intensity"], label="Normalize")
    robjects.r.assign("normalizeOption", normalizeOption_py)

    st.sidebar.subheader("1-2 Filter on missing values:")
    # Filter for proteins that are identified in all replicates of at least one condition
    maxReplicate_py = robjects.r("maxReplicate")
    nThr_py = st.sidebar.slider('Filter for proteins that are identified in all replicates of at least one condition: ', min_value = 0,max_value = maxReplicate_py[0] ,value = 0, step=1, format="%d")

    st.sidebar.subheader("1-5. Differential enrichment analysis")
    alpha_py = st.sidebar.slider("alpha: ",min_value = 0.0,max_value = 1.0 ,value = 1.0, step=0.01, format="%f")
    lfc_py = st.sidebar.slider("lfc = log2(value): ", min_value=0.0, max_value=2.0, value=1.5, step=0.1, format="%f")

    with st.sidebar:
        #submitted = st.form_submit_button("Done")
        submitted = st.button("Done")

    if submitted:

        r_normalize_function()
        st.header('2. Filter on missing values:')
        r_plot_frequency_1_1()
        r_plot_numbers_filter_missval_1_2(nThr_py)
        r_plot_coverage_1_3()
        st.header("3. Normalization")

        r_plot_normalization_1_4()

        st.header("4. Impute data for missing values: ")
        r_plot_heatmap_1_5()

        st.header("5. Differential enrichment analysis")

        r_plot_pca_1_6(control_py, alpha_py, lfc_py)

        try :
            r_plot_heatmap_dep_1_7()
        except:
            st.error("Error! 需要調整alpha, lfc的值 (如: alpha = 1, lfc = 1)",icon="🚨")

        st.header("6. Volcano plots of specific contrasts:")
        r_plot_volcano_1_8(contrast_py)

        st.header("7. Barplots of a protein of interest: ")
        st.sidebar.subheader("7. Barplots of a protein of interest")
        r_plot_single_1_9()

        try :
            st.header("8. Frequency plot of significant proteins and overlap of conditions")
            r_plot_cond_1_10()
        except:
            st.error('condition要有三組以上', icon="🚨")
        # DOSE
        #r_uniprotAPI()
        #st.sidebar.subheader("9. Configure data for analysis (DOSE)")

        #r_init_DOSE_data()
        #r_convert_species_gene(rationName_py)
        r_test_file() # 測試用
        st.header("10. Bar Plot (DOSE)")
        r_plot_barplot_2_1()
        st.header("11. Dot plot")
        r_plot_dotplot_2_2()
        st.header("12. Gene-Concept Network")
        r_plot_cnetplot_2_3()
        st.header("13. Heatmap-like functional classification")
        r_plot_heatplot_2_4()

        st.header("14. Enrichment Map")
        r_plotenrichment_map_2_5()

        st.header("15. Biological theme comparison")
        #st.warning('有問題，待處理 error compareCluster(data, fun = "enrichKEGG", organism = "hsa", pvalueCutoff = 0.05) : No enrichment found in any of gene cluster, please check your input.')
        r_plot_emapplot_2_6() #error compareCluster(data, fun = "enrichKEGG", organism = "hsa", pvalueCutoff = 0.05) : No enrichment found in any of gene cluster, please check your input.

        st.header("16. UpSet Plot")
        r_plot_upseplot_2_7()
        #st.warning("有問題，待處理 (gseKEGG: Error in check_gene_id(geneList, geneSets) : --> No gene can be mapped....)")

        pvalue_2_8_py = st.sidebar.slider("pvalue: ", min_value=0.001, max_value=1.0, step=0.001, value=0.05)
        robjects.r.assign("pvalue_2_8", pvalue_2_8_py)
        r_plot_upsetplot_with_splider_2_8() #error

        st.header("17. ridgeline plot for expression distribution of GSEA result")
        r_plot_ridgeplot_2_9()

        st.header("18. running score and preranked list of GSEA result")
        st.warning("有問題, 待處理")
        #r_plot_gseaplot_2_10()
        with st.sidebar:
            st.subheader("19. Download result file")
            download_button("./file/dep_output.csv", "Download dep_output.csv")
            download_button("./file/uniprot_entrez.csv", "Download uniprot_entrez.csv")
            download_button("./file/dep_output_result.csv", "Download dep_output_result.csv")

all()