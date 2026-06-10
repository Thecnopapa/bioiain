FROM debian:trixie-slim AS build

RUN uname -m > /arch

WORKDIR /app


RUN apt-get update && apt-get install -y python3 python3-pip python3.13-venv

RUN apt-get clean

RUN python3 -m venv venv

ENV VIRTUAL_ENV="./venv"
ENV PATH="./venv/bin:$PATH"
RUN which python
RUN pwd
RUN ls

RUN pip3 install -r ./requirements.txt

COPY ./app .


RUN chmod +x ./start-server.sh

RUN echo "executing: start-server.sh ..."
CMD ["./start-server.sh"]