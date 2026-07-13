#!/bin/bash


# PROJECT SPECIFIC AFTER THIS

testing(){
  cd $FOLDER_PATH/test && python test.py "$@"
}

upload(){
  cd $FOLDER_PATH | exit 0
  python3 -m build
  python3 -m twine upload --repository pypi dist/* --verbose --skip-existing
  rm dist/*
}
upload-test(){
  cd $FOLDER_PATH | exit 0
  python3 -m build
  python3 -m twine upload --repository testpypi dist/* --verbose --skip-existing
  rm dist/*
}

upload-gcloud(){
  cd $FOLDER_PATH | exit 0
  python3 -m build
  python3 -m twine upload --repository-url  https://europe-west1-python.pkg.dev/iainvisa/python/ dist/* --verbose
  rm dist/*
}

upload-all(){
  upload
  upload-gcloud
}


update(){
    echo '$ pip install bioiain -U'
    pip install bioiain -U
}
update-test(){
    echo '$ pip install -i https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ bioiain -U'
    pip install -i https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ bioiain -U
}

update-gcloud(){
    gcloud config set project iainvisa
    echo '$ pip install -i https://europe-west1-python.pkg.dev/iainvisa/python/simple/ --extra-index-url https://pypi.org/simple/ bioiain -U'
    pip install -i  https://europe-west1-python.pkg.dev/iainvisa/python/simple/ --extra-index-url https://pypi.org/simple/ bioiain -U
}

update-all(){
    pip bioiain -U
    pip bioiain[ml] -U
    pip bioiain[conda] -U
}
board(){
  tensorboard --logdir=. &
}

profile(){
  python -m cProfile -s time test.py -t --epochs 1 > performance.txt
}
