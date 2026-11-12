"""What the repository layer refuses, and the router turns into 404 and 400."""


class NoSuchRepository(Exception):
    """Nothing readable at this path, or what is there is not a repository."""


class Refused(Exception):
    """The request asks for something the repository layer will not do."""
