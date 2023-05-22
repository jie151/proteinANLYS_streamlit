import streamlit as st
import rpy2.robjects as robjects
from streamlit.server.server import Server

# 创建一个全局变量存储每个用户的R实例
global_r_instances = {}

# 获取当前用户的R实例
def get_r_instance():
    session_id = st.report_thread.get_report_ctx().session_id
    if session_id not in global_r_instances:
        # 初始化R实例
        r = robjects.r
        global_r_instances[session_id] = r
    return global_r_instances[session_id]

# Streamlit应用程序
def streamlit_app():
    st.write("# Rpy2 with Streamlit Example")

    # 获取当前用户的R实例
    r = get_r_instance()

    # 在R实例上执行操作
    result = r('1 + 2')
    python_result = int(result[0])

    # 在Streamlit中展示结果
    st.write("R 结果：", result)
    st.write("Python 结果：", python_result)

# 主函数
def main():
    # 使用SessionState模块创建独立会话状态
    session_state = SessionState.get(report_context=None)

    # 检查是否在Streamlit服务器上运行
    if isinstance(session_state.report_context.request, Server):
        streamlit_app()

if __name__ == "__main__":
    main()
