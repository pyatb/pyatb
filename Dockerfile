# Dockerfile
FROM quay.io/pypa/manylinux2014_x86_64

# 安装必要的库
RUN yum install -y openblas-devel lapack lapack-devel

# 配置环境变量以确保链接器可以找到库
ENV LD_LIBRARY_PATH=/usr/lib64:/usr/local/lib:$LD_LIBRARY_PATH

# 将容器默认命令设置为进入bash shell，以便后续调试
CMD ["/bin/bash"]
