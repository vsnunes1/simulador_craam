FROM ubuntu

WORKDIR /code

ENV PORT 80

COPY . /code

RUN apt-get update
RUN apt install python2 -y
RUN apt install python-tk -y
RUN apt install ffmpeg -y
RUN ln -s /usr/bin/ffmpeg /usr/bin/avconv
RUN python2 /code/get-pip.py
RUN pip install django
RUN pip install scipy
RUN pip install numpy
RUN pip install matplotlib 

CMD ["python2", "/code/manage.py", "runserver", "0.0.0.0:80"]