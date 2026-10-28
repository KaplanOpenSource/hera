"""The app is assembled from the routers, with the catch-all last."""
from args_settings import parse_args_settings
from server import create_app


def _paths():
    # The OpenAPI schema is the documented view of what the app serves, in the
    # order the routes were added.
    return list(create_app(parse_args_settings([])).openapi()["paths"].keys())


def test_every_router_is_included():
    paths = _paths()

    for path in [
        "/healthz",
        "/cors",
        "/init",
        "/exec",
        "/start_workflow",
        "/workflow_status/{token}",
        "/jupyter/ensure",
        "/node-catalog",
        "/node-params/validate",
        "/file/{file_path}",
        "/",
    ]:
        assert path in paths, f"{path} is not served"


def test_the_spa_catch_all_comes_last():
    # Its /{full_path} route matches anything, so a route added after it would
    # never be reached.
    paths = _paths()

    assert paths[-1] == "/{full_path}"


def test_cors_is_the_outer_middleware():
    # The error reporter must sit inside CORS, so its 500 answers carry the
    # CORS header. Last added is outermost, so CORS is added last.
    app = create_app(parse_args_settings([]))

    assert "CORSMiddleware" in str(app.user_middleware[0])
