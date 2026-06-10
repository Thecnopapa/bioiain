FROM debian:trixie-slim AS build

RUN uname -m > /arch

WORKDIR /app


RUN apt-get update && apt-get install -y python3 python3-pip python3.13-venv

RUN apt-get clean

RUN python3 -m venv venv

ENV VIRTUAL_ENV="./venv"
ENV PATH="./venv/bin:$PATH"
RUN which python

RUN pip3 install bioiain
RUN pip3 install pdoc3


COPY ./app .

RUN chmod +x ./build-wiki.sh
RUN ls

RUN echo "executing: build-wiki.sh ..."
CMD ["./build-wiki.sh"]