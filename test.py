from streamlit.runtime.scriptrunner.script_run_context import get_script_run_ctx
import streamlit as st
import os
os.environ['R_HOME']= 'C:\\Program Files\\R\\R-4.2.3'
os.environ['PATH'] += os.pathsep + 'C:\\Program Files\\R\\R-4.2.3\\bin\\X64\\'
os.environ['PATH'] += os.pathsep + 'C:\\Program Files\\R\\R-4.2.3\\'

import os

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

import threading


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



def long_r_function() :
    session_id = get_script_run_ctx().session_id
    from rpy2 import robjects
    robjects.r.assign("id", session_id)

    for i in range(1, 10):
        time.sleep(1)

        st.write(f"{i}. {session_id}")

        output = st.empty()

        with st_capture(output.code):
            robjects.r("print(id)")

long_r_function()