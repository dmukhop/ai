FROM ubuntu

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    build-essential \
    curl \
    git \
    #python3-dev \
    python-pip \
    python-setuptools \
    python-dev \
  && rm -rf /var/lib/apt/lists/*
#RUN apt-get install python-pip
RUN pip install --upgrade pip
RUN pip install h5py pandas
#RUN pip install theano
RUN pip install xlrd
#RUN pip install --upgrade -I setuptools \
#  && pip install --upgrade \
#    keras

RUN pip install  \
    matplotlib \
    seaborn

RUN pip install scikit-learn tables
#RUN pip install --upgrade pip
RUN pip install 'ipython<6'
#RUN pip3 install --user --upgrade tensorflow  # install in $HOME

RUN pip install jupyter

VOLUME /w210
WORKDIR /w210
EXPOSE 8888

#ENV KERAS_BACKEND=TensorFlow

#  CMD jupyter notebook --no-browser --ip=0.0.0.0 --NotebookApp.token= --allow-root
CMD jupyter notebook --no-browser --ip=0.0.0.0 --allow-root
