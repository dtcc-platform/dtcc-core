import os
import ssl

import requests


def download_ssl_context():
    """Return the verifying SSL context for data downloads.

    It trusts Python's default CA locations. Unless ``SSL_CERT_FILE`` or
    ``SSL_CERT_DIR`` selects the trusted CAs, it also trusts the bundle that
    ``requests`` uses for lookups, so a Python with no system CA store can
    download whatever it can look up.
    """
    context = ssl.create_default_context()
    if not (os.environ.get("SSL_CERT_FILE") or os.environ.get("SSL_CERT_DIR")):
        context.load_verify_locations(cafile=requests.certs.where())
    return context
