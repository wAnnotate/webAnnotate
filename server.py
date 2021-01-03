from flask import Flask,flash, request, jsonify, render_template,Response,redirect,send_from_directory



app = Flask(__name__)

@app.route("/")
def index():
    return Response(response="Hello! This is Annotator")


if __name__ == "__main__":
    app.run()



