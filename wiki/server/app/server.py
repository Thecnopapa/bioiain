import flask, requests, os, sys, json
from werkzeug.utils import secure_filename

app = flask.Flask(__name__)
app.secret_key = os.environ['FLASK_KEY']

version_list=sorted(os.listdir("/docs"), reverse=True)
latest_version=version_list[0]

@app.route('/<version>/<path:path>')
def page(version=None, path=None):
    path = secure_filename(path)
    version = version or latest_version
    print(version, path)
    if version not in version_list:
        return flask.abort(404)

    return flask.send_from_directory('/docs/{version}', path)

if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=int(os.environ.get("PORT", 8080)))