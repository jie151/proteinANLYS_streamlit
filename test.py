from streamlit.runtime.scriptrunner.script_run_context import get_script_run_ctx
import streamlit as st
import os
os.environ['R_HOME']= 'C:\\Program Files\\R\\R-4.2.3'
os.environ['PATH'] += os.pathsep + 'C:\\Program Files\\R\\R-4.2.3\\bin\\X64\\'
os.environ['PATH'] += os.pathsep + 'C:\\Program Files\\R\\R-4.2.3\\'

import os
import rpy2.robjects as robjects
import streamlit as st
import pandas as pd
from PIL import Image
from contextlib import contextmanager, redirect_stdout
from io import StringIO
import base64
import uuid
import re
import sys
import time

# def upload_file():
#     # 上傳檔案
#     uploaded_file = st.sidebar.file_uploader('Upload a TXT/CSV file', type=['csv', 'txt'], accept_multiple_files=False)
#     if uploaded_file is not None:
#         # 將上傳的檔案儲存下來 (由R開啟)
#         filename = "./uploadFile.txt"
#         with open(filename,"wb") as f:
#             f.write( uploaded_file.getbuffer())

#         sep_word = "," if uploaded_file.type == "text/csv" else "\t"
#         data = pd.read_csv(uploaded_file, sep = sep_word, error_bad_lines=False, dtype=object)
#     else:
#         filename = "./proteinGroups_HsinYuan_Rat.txt"
#         data = pd.read_csv(filename, sep="\t", dtype=object)
#     return filename, data



def get_session_id():
    session_id = get_script_run_ctx().session_id
    return session_id

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

# filename, data = upload_file()


# id = get_session_id()
# robjects.r.assign("id_r", id)
# robjects.r('''
#     print(id_r)
# ''')

# for i in range(10):

#     st.write(f"{i} python: {id}")
#     output = st.empty()
#     with st_capture(output.code):
#         robjects.r('''
#             print( id_r)
#         ''')

#     time.sleep(1)


import concurrent.futures
from concurrent.futures import ProcessPoolExecutor
import time

import streamlit as st


if 'save' not in st.session_state:
    st.session_state.save = []


def task(v):
    """session state does not work here"""
    robjects.r.assign("id_r", v)
    print("python ", v)
    for i  in range(10):
        robjects.r('''print(id_r)''')
    time.sleep(1)
    return v


if __name__ == '__main__':
    num_workers = 2
    id = get_session_id()
    jobs = [id, 2, 3, 4, 5, 6, 7, 8, 9]
    processed_jobs = []

    start = st.button('start work')

    if start:
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            for j in jobs:
                pj = executor.submit(task, j)
                processed_jobs.append(pj)

            for future in concurrent.futures.as_completed(processed_jobs):
                try:
                    res = future.result()
                    st.write(f'res: {res}')

                    # Incrementally save the completed task so far.
                    st.session_state.save.append(res)

                except concurrent.futures.process.BrokenProcessPool as ex:
                    raise Exception(ex)

    if len(st.session_state.save):
        st.write('#### Completed Jobs')
        st.write(f'{st.session_state.save}')