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
def page(version=None, path=None, filename=None):
    global version_list
    print("Fetching page...")
    if version is None:
        version_list=sorted(os.listdir("/docs"), reverse=True)
        print(version, version in version_list)
        latest_version=version_list[0]
        version = latest_version
    else:
        version = secure_filename(version)

    if version not in version_list:
        return f"Version ({version}) not found"

    docspath = os.path.join("/docs", version)
    #print(filepath)
    if path is None or path == "":
        print("No path provided")
        filepath = docspath
    else:
        #print(path)
        path = "/".join([secure_filename(pp) for pp in path.split("/")])
        #print(path)
        print("path:"path)
        filepath = os.path.join(docspath, path)
    


    if filename is None or filename == '':
        filename = "index.html"
    else:
        filename = secure_filename(filename)

    print("filename:", filename)
    filepath = os.path.join(filepath, filename)
    print("final_path:", filepath)
    if not os.path.exists(filepath):
        return f"File not found: {filepath}"

    with open(filepath) as f:
        resp = flask.make_response(f.read())

    return resp



if __name__ == "__main__":
    app.run(debug=False, host="0.0.0.0", port=int(os.environ.get("PORT", 8080)))
