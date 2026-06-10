import flask, requests, os, sys, json
from werkzeug.utils import secure_filename
from werkzeug.exceptions import NotFound

app = flask.Flask(__name__)
app.secret_key = os.environ['FLASK_KEY']

version_list=sorted(os.listdir("/docs"), reverse=True)

latest_version=version_list[0]
print("Version list: {}".format(version_list))
print("Latest version: {}".format(latest_version))


@app.route("/test")
def test():
    print("Test")
    return "Test"
@app.route('/')
def index():
    return flask.redirect(f"/{latest_version}/index.html")

@app.route('/<version>/')
@app.route('/<version>/<path:path>/')
@app.route('/<version>/<filename>')
@app.route('/<version>/<path:path>/<filename>')
def debugpage(version=None, path=None, filename=None):
    return f"{version}/{path}/{filename}"

def page(version=None, path=None, filename=None):
    print("Fetching page...")
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

    filepath = os.path.join(filepath, filename)

    print(filepath)
    if not os.path.exists(filepath):
        print(f"File ({filepath}/{filename}) not found")
        return flask.abort(404)

    with open(filepath) as f:
        resp = flask.make_response(f.read())
    return resp


if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=int(os.environ.get("PORT", 8080)))