#!/usr/bin/env python3
"""
results-receiver.py -- stores catalogue.html's DCP results on people, next to the dataset.

Workers POST one application/x-www-form-urlencoded form per slice (DCP remote result storage:
the job's params plus contentType, element and content). The job params say what it is:

  kind=fit                 content is a runInitialFit report carrying galaxyId
                           -> galaxies/<id>/initial-fit.json; a new fit deletes that galaxy's
                              old bootstraps, since they were measured against the old fit
  kind=bootstrap&galaxyId  content is a runBootstrapRealization report
                           -> galaxies/<id>/bootstrap/<realizationIndex>.json

The reply, which DCP hands back to the job as the slice result, is a short acknowledgement
naming the stored file, so the page fetches the fit from people rather than through DCP.
index.json (rewritten at most every 2 s) summarises every galaxy for the page's status column
and resume.

Only POST <prefix>/<token> is accepted, where <token> is the content of TOKEN_FILE (reread
when it changes). CORS is answered here, so the nginx location only needs to proxy.
Listens on 127.0.0.1 only. Python 3 standard library only.
"""
import hmac
import json
import os
import re
import shutil
import threading
import time
from datetime import datetime, timezone
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs

HOME = os.path.expanduser('~')
STORE = os.environ.get('KIDNAS_RESULTS_DIR', os.path.join(HOME, 'public_html', '3kidnas-results'))
TOKEN_FILE = os.environ.get('KIDNAS_RESULTS_TOKEN_FILE', os.path.join(HOME, '3kidnas-results', 'token'))
PORT = int(os.environ.get('KIDNAS_RESULTS_PORT', '8661'))
MAX_BODY = 64 * 1024 * 1024

GALAXY_RE = re.compile(r'^WALLABY_J[0-9]{6}[+-][0-9]{6}$')
FIT_SUMMARY = ('FITAchieved', 'INCLINATION', 'POSITIONANGLE', 'VSYS', 'error')

lock = threading.Lock()
index = {'galaxies': {}}
dirty = threading.Event()
token_cache = {'mtime': None, 'token': None}


def now():
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


def current_token():
    mtime = os.stat(TOKEN_FILE).st_mtime_ns
    if mtime != token_cache['mtime']:
        with open(TOKEN_FILE) as f:
            token_cache['token'] = f.read().strip()
        token_cache['mtime'] = mtime
    token = token_cache['token']
    return token if token and len(token) >= 32 else None


def write_json(path, value):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(value, f, separators=(',', ':'))
    os.chmod(tmp, 0o644)
    os.replace(tmp, path)


def first(v):
    return v[0] if isinstance(v, list) and v else v


def load_index():
    """Rebuild the index from the files on disk, so a restart loses nothing."""
    galaxies = {}
    root = os.path.join(STORE, 'galaxies')
    for gid in sorted(os.listdir(root)) if os.path.isdir(root) else []:
        entry = galaxies[gid] = {'bootstraps': []}
        fit_path = os.path.join(root, gid, 'initial-fit.json')
        if os.path.isfile(fit_path):
            with open(fit_path) as f:
                fit = json.load(f)
            entry['fit'] = {k: first(fit.get(k)) for k in FIT_SUMMARY}
            entry['fit']['at'] = datetime.fromtimestamp(os.path.getmtime(fit_path), timezone.utc).isoformat(timespec='seconds')
        bs_dir = os.path.join(root, gid, 'bootstrap')
        if os.path.isdir(bs_dir):
            entry['bootstraps'] = sorted(int(n[:-5]) for n in os.listdir(bs_dir) if re.match(r'^\d+\.json$', n))
    index['galaxies'] = galaxies


def index_writer():
    while True:
        dirty.wait()
        time.sleep(2)
        dirty.clear()
        with lock:
            snapshot = json.loads(json.dumps(index))
        snapshot['updatedAt'] = now()
        write_json(os.path.join(STORE, 'index.json'), snapshot)


def store_fit(report):
    gid = report.get('galaxyId')
    if not isinstance(gid, str) or not GALAXY_RE.match(gid):
        raise ValueError('bad galaxyId')
    folder = os.path.join(STORE, 'galaxies', gid)
    with lock:
        shutil.rmtree(os.path.join(folder, 'bootstrap'), ignore_errors=True)
        write_json(os.path.join(folder, 'initial-fit.json'), report)
        summary = {k: first(report.get(k)) for k in FIT_SUMMARY}
        summary['at'] = now()
        index['galaxies'][gid] = {'fit': summary, 'bootstraps': []}
    dirty.set()
    return {'stored': True, 'kind': 'fit', 'galaxyId': gid, 'FITAchieved': bool(report.get('FITAchieved')),
            'error': str(report.get('error') or '')[:300] or None,
            'url': f'galaxies/{gid}/initial-fit.json'}


def store_bootstrap(gid, report):
    if not isinstance(gid, str) or not GALAXY_RE.match(gid):
        raise ValueError('bad galaxyId')
    n = report.get('realizationIndex')
    if not isinstance(n, int) or not 0 <= n < 100000:
        raise ValueError('bad realizationIndex')
    with lock:
        write_json(os.path.join(STORE, 'galaxies', gid, 'bootstrap', f'{n}.json'), report)
        entry = index['galaxies'].setdefault(gid, {'bootstraps': []})
        if n not in entry['bootstraps']:
            entry['bootstraps'].append(n)
    dirty.set()
    return {'stored': True, 'kind': 'bootstrap', 'galaxyId': gid, 'realizationIndex': n,
            'sofiaFailed': bool(report.get('sofiaFailed')), 'url': f'galaxies/{gid}/bootstrap/{n}.json'}


class Handler(BaseHTTPRequestHandler):
    server_version = '3kidnas-results'

    def log_request(self, code='-', size='-'):
        # The default line includes the path, and the path holds the upload token.
        self.log_message('%s %s', self.command, code)

    def cors(self):
        self.send_header('Access-Control-Allow-Origin', '*')
        self.send_header('Access-Control-Allow-Methods', 'POST')
        self.send_header('Access-Control-Allow-Headers', 'Content-Type')
        self.send_header('Access-Control-Max-Age', '86400')

    def reply(self, status, body):
        data = json.dumps(body).encode()
        self.send_response(status)
        self.cors()
        self.send_header('Content-Type', 'application/json')
        self.send_header('Content-Length', str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def do_OPTIONS(self):
        self.send_response(204)
        self.cors()
        self.end_headers()

    def do_GET(self):
        # Lets the page check that the upload endpoint is reachable, without the token.
        if self.path.rstrip('/').endswith('/ping'):
            return self.reply(200, {'ok': True})
        self.reply(405, {'error': 'POST only'})

    def do_POST(self):
        try:
            token = current_token()
            supplied = self.path.split('?')[0].rstrip('/').rsplit('/', 1)[-1]
            if not token or not hmac.compare_digest(supplied, token):
                return self.reply(404, {'error': 'not found'})
            length = int(self.headers.get('Content-Length') or 0)
            if length <= 0 or length > MAX_BODY:
                return self.reply(413, {'error': 'missing or oversized body'})
            form = parse_qs(self.rfile.read(length).decode('ascii'), max_num_fields=32)
            kind = form.get('kind', [''])[0]
            if 'content' not in form:
                return self.reply(400, {'error': 'slice failed', 'detail': form.get('error', [''])[0][:300]})
            report = json.loads(form['content'][0])
            if kind == 'fit':
                ack = store_fit(report)
            elif kind == 'bootstrap':
                ack = store_bootstrap(form.get('galaxyId', [''])[0], report)
            else:
                return self.reply(400, {'error': 'unknown kind'})
            self.log_message('stored %s %s', kind, ack['galaxyId'])
            self.reply(200, ack)
        except Exception as e:  # never echo request content back
            self.log_message('rejected: %s', type(e).__name__)
            self.reply(400, {'error': type(e).__name__})


if __name__ == '__main__':
    os.makedirs(os.path.join(STORE, 'galaxies'), exist_ok=True)
    os.chmod(STORE, 0o755)
    load_index()
    write_json(os.path.join(STORE, 'index.json'), dict(index, updatedAt=now()))
    threading.Thread(target=index_writer, daemon=True).start()
    ThreadingHTTPServer(('127.0.0.1', PORT), Handler).serve_forever()
