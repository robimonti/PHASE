"""Updater HTTPS must work with missing default roots without weakening TLS."""
import hashlib
import io
import json
import ssl
import sys
from types import SimpleNamespace
import pytest
import phase_update


def test_certifi_supplies_roots_when_python_default_store_is_empty(monkeypatch):
    pytest.importorskip('certifi')
    empty = ssl.SSLContext(ssl.PROTOCOL_TLS_CLIENT)
    assert empty.cert_store_stats()['x509_ca'] == 0
    monkeypatch.setattr(phase_update.ssl, 'create_default_context', lambda: empty)
    context = phase_update.tls_context()
    assert context is empty
    assert context.cert_store_stats()['x509_ca'] > 0
    assert context.check_hostname
    assert context.verify_mode == ssl.CERT_REQUIRED


def test_system_trust_remains_verified_without_certifi(monkeypatch):
    default = ssl.create_default_context()
    monkeypatch.setattr(phase_update.ssl, 'create_default_context', lambda: default)
    monkeypatch.setitem(sys.modules, 'certifi', None)
    assert phase_update.tls_context() is default
    assert default.check_hostname and default.verify_mode == ssl.CERT_REQUIRED


def test_invalid_ca_bundle_is_not_silently_ignored(monkeypatch, tmp_path):
    monkeypatch.setitem(sys.modules, 'certifi', SimpleNamespace(where=lambda: str(tmp_path / 'missing.pem')))
    with pytest.raises(OSError):
        phase_update.tls_context()


def test_metadata_and_download_use_verified_context_and_keep_digest_check(monkeypatch, tmp_path):
    context = ssl.create_default_context()
    monkeypatch.setattr(phase_update, 'tls_context', lambda: context)
    payload = b'PHASE archive fixture'
    digest = 'sha256:' + hashlib.sha256(payload).hexdigest()
    url = 'https://github.com/robimonti/PHASE/releases/download/v7.1.12/phase7-engine.zip'
    metadata = {'tag_name': 'v7.1.12', 'assets': [{
        'name': 'phase7-engine.zip', 'browser_download_url': url, 'digest': digest,
    }]}
    calls = []

    def urlopen(request, *, timeout, context):
        calls.append((request.full_url, timeout, context))
        return io.BytesIO(json.dumps(metadata).encode() if request.full_url == phase_update.API_URL else payload)

    monkeypatch.setattr(phase_update.urllib.request, 'urlopen', urlopen)
    release = phase_update.release_info()
    output = tmp_path / 'archive.zip'
    phase_update.download_verified(release['url'], release['digest'], output)
    assert output.read_bytes() == payload
    assert calls == [(phase_update.API_URL, 15, context), (url, 60, context)]
    with pytest.raises(RuntimeError, match='checksum'):
        phase_update.download_verified(url, 'sha256:' + '0' * 64, output)
