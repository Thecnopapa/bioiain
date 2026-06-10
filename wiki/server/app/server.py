import flask, requests, os, sys, json
from werkzeug.utils import secure_filename
from werkzeug.exceptions import NotFound

app = flask.Flask(__name__)
app.secret_key = os.environ['FLASK_KEY']

version_list=sorted(os.listdir("/docs"), reverse=True)

latest_version=version_list[0]
print("Version list: {}".format(version_list))
print("Latest version: {}".format(latest_version))

@app.route('/')
def index():
    return flask.redirect(f"/{latest_version}/index.html")

@app.route('/<version>/')
def redirect_version_to_index(version):
    return flask.redirect(f"/{version}/index.html")

@app.route('/<version>/<path:path>/')
def redirect_to_index(version, path):
    return flask.redirect(f"/{version}/{path}/index.html")

@app.route('/<version>/<filename>')
@app.route('/<version>/<path:path>/<filename>')
def page(version=None, path=None, filename=None):
    if version is None:
        version = latest_version
    else:
        version = secure_filename(version)
    print(version)
    if version not in version_list:
        print(f"Version ({version})not found")
        return flask.abort(404)
    filepath = f"/docs/{version}"

    print(path)
    if path is None or path == '':
        pass
    else:
        path = secure_filename(path)
        filepath = os.path.join(filepath, path)

    if filename is None or filename == '':
        filename = "index.html"

    print(version, filepath, filename)

    try:
        resp = flask.send_from_directory(filepath, filename)
    except NotFound:
        print(f"File ({filepath}/{filename}) not found")
        return flask.abort(404)
    return resp

if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=int(os.environ.get("PORT", 8080)))