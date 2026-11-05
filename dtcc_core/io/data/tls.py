import os
import ssl

import requests


def download_ssl_context():
    """Return the verifying SSL context for data downloads.

    It trusts Python's default CA locations, which ``SSL_CERT_FILE`` and
    ``SSL_CERT_DIR`` select, and the bundle ``requests`` uses for lookups:
    ``REQUESTS_CA_BUNDLE`` or ``CURL_CA_BUNDLE`` when one is set, otherwise
    the ``certifi`` bundle unless ``SSL_CERT_FILE`` or ``SSL_CERT_DIR`` is set.
    """
    context = ssl.create_default_context()
    bundle = os.environ.get("REQUESTS_CA_BUNDLE") or os.environ.get("CURL_CA_BUNDLE")
    if bundle:
        if os.path.isdir(bundle):
            context.load_verify_locations(capath=bundle)
        else:
            context.load_verify_locations(cafile=bundle)
    elif not (os.environ.get("SSL_CERT_FILE") or os.environ.get("SSL_CERT_DIR")):
        context.load_verify_locations(cafile=requests.certs.where())
    return context
