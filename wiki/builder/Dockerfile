FROM debian:trixie-slim AS build

RUN uname -m > /arch

WORKDIR /app


RUN apt-get update && apt-get install -y python3 python3-pip python3.13-venv

RUN apt-get clean

RUN python3 -m venv venv

ENV VIRTUAL_ENV="./venv"
ENV PATH="./venv/bin:$PATH"
RUN which python

RUN pip3 install torch torchvision --index-url https://download.pytorch.org/whl/cpu
RUN pip3 install bioiain[ml,aleph]
RUN pip3 install pdoc==16.0.0
RUN pip3 cache purge
RUN which pdoc

RUN echo "$(find "./venv"|grep 'pdoc/__init__.py')"

COPY ./app .
#RUN echo "./pdoc-init-fixed.py --> $PDOC_INIT"
COPY ./pdoc-init-fixed.py .
RUN mv ./pdoc-init-fixed.py "./venv/lib/python3.13/site-packages/pdoc/__init__.py"
#RUN mv ./pdoc-init-fixed.py "$(find "./venv"|grep 'pdoc/__init__.py')"

RUN chmod +x ./build-wiki.sh

RUN echo "executing: build-wiki.sh ..."
CMD ["./build-wiki.sh"]