import os
from flask import Flask
from flask_cors import CORS
from api import patients, drugs, predict, query

def create_app(test_config=None):
    """Create and configure an instance of the Flask application."""
    app = Flask(__name__, instance_relative_config=True)
    CORS(app)

    app.config.from_mapping(
        # a default secret that should be overridden by instance config
        SECRET_KEY="dev",
        # store the database in the instance folder
        DATABASE=os.path.join(app.instance_path, "flaskr.sqlite"),
    )

    if test_config is None:
        # load the instance config, if it exists, when not testing
        app.config.from_pyfile("config.py", silent=True)
    else:
        # load the test config if passed in
        app.config.update(test_config)

    # ensure the instance folder exists
    try:
        os.makedirs(app.instance_path)
    except OSError:
        pass

    #keep for testing
    @app.route("/hello")
    def hello():
        return "Hello, World!"

    # register the database commands
    #from flaskr import db
    #db.init_app(app)

    # apply the blueprints to the app
    app.register_blueprint(patients.bp)
    app.register_blueprint(drugs.bp)

    app.register_blueprint(predict.bp)
    app.register_blueprint(query.bp)
    return app
