from flask import Flask, flash, request, jsonify, render_template, Response, redirect, send_from_directory
from perlfunc import perl5lib, perlfunc, perlreq

app = Flask(__name__)


@perlfunc
@perlreq('modules/annotate_variation.pl')
def annotate_variation():
    pass


@perlfunc
@perlreq('modules/coding_change.pl')
def coding_change():
    pass


@perlfunc
@perlreq('modules/convert2annovar.pl')
def convert2annovar():
    pass


@perlfunc
@perlreq('modules/retrieve_seq_from_fasta.pl')
def retrieve_seq_from_fasta():
    pass


@perlfunc
@perlreq('modules/table_annovar.pl')
def table_annovar():
    pass


@perlfunc
@perlreq('modules/variants_reduction.pl')
def variants_reduction():
    pass


@app.route("/")
def index():
    return Response(response="Hello! This is Annotator")


if __name__ == "__main__":
    app.run()
