import flask, requests, os, sys, json
from werkzeug.utils import secure_filename

app = flask.Flask(__name__)
app.secret_key = os.environ['FLASK_KEY']

version_list=sorted(os.listdir("/docs"), reverse=True)
latest_version=version_list[0]
print("Version list: {}".format(version_list))
print("Latest version: {}".format(latest_version))

@app.route('/')
def index():
    return flask.redirect(f"/{latest_version}/index.html")

@app.route('/<version>/<path:path>/<filename>')
def page(version=None, path=None, filename=None):
    print(version)
    if version not in version_list:
        return flask.abort(404)
    filepath = f"/docs/{version}"
    path = secure_filename(path)
    print(path)
    if path is None or path == '':
        pass
    else:
        filepath = os.path.join(filepath, path)

    if filename is None or filename == '':
        filename = "index.html"
    version = version or latest_version
    print(version, filepath, filename)

    return flask.send_from_directory(filepath, filename)

if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=int(os.environ.get("PORT", 8080)))