#!/usr/bin/env python3
"""
results-receiver.py -- stores catalogue.html's DCP results on people, next to the dataset.

Workers POST one application/x-www-form-urlencoded form per slice (DCP remote result storage:
the job's params plus contentType, element and content). The job params say what it is:

  kind=fit                 content is a runInitialFit report carrying galaxyId
                           -> galaxies/<id>/initial-fit.json, and its bootstrapPayload (built by
                              the worker) -> galaxies/<id>/bootstrap-payload.json, which is all
                              the page fetches to launch the bootstrap job; a new fit deletes
                              that galaxy's old bootstraps, so a rerun overwrites
  kind=bootstrap&galaxyId  content is a runBootstrapRealization report
                           -> one line appended to galaxies/<id>/bootstrap.ndjson (a retried
                              slice can add a second line for the same realizationIndex;
                              readers keep the last)

Every stored result also appends one line to events.ndjson (arrival time, kind, galaxy,
realization, the worker's own timings, size), which is what throughput analysis reads.

The reply, which DCP hands back to the job as the slice result, is a short acknowledgement
naming the stored file, so the page fetches the fit from here rather than through DCP.
index.json (rewritten at most every INDEX_EVERY s) has each galaxy's fit summary and
realization count, for the page's status column and resume.

Only POST <prefix>/<token> is accepted, where <token> is any line of TOKEN_FILE (reread
when it changes), so a new key can be added without breaking jobs still using an older one.
At most BUSY_MAX posts are parsed at once, which bounds memory (a big fit is ~125 MB); nginx
buffers the rest. Results are refused once the disk has less than MIN_FREE_GB free.
CORS is answered here, so the nginx location only needs to proxy.
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
MAX_BODY = 256 * 1024 * 1024
INDEX_EVERY = 30
BUSY_MAX = 2
MIN_FREE_GB = float(os.environ.get('KIDNAS_RESULTS_MIN_FREE_GB', '5'))

GALAXY_RE = re.compile(r'^WALLABY_J[0-9]{6}[+-][0-9]{6}$')
FIT_SUMMARY = ('FITAchieved', 'INCLINATION', 'POSITIONANGLE', 'VSYS', 'error')

lock = threading.Lock()
index = {'galaxies': {}}
dirty = threading.Event()
token_cache = {'mtime': None, 'tokens': []}
busy = threading.BoundedSemaphore(BUSY_MAX)
log_lock = threading.Lock()


def now():
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


def current_tokens():
    mtime = os.stat(TOKEN_FILE).st_mtime_ns
    if mtime != token_cache['mtime']:
        with open(TOKEN_FILE) as f:
            token_cache['tokens'] = [t.strip() for t in f if len(t.strip()) >= 32]
        token_cache['mtime'] = mtime
    return token_cache['tokens']


def token_ok(supplied):
    # Compare against every key so the time taken does not depend on which one matched.
    return sum(hmac.compare_digest(supplied, t) for t in current_tokens()) > 0


def disk_full():
    st = os.statvfs(STORE)
    return st.f_bavail * st.f_frsize < MIN_FREE_GB * 1e9


def append_line(path, value):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'a') as f:
        f.write(json.dumps(value, separators=(',', ':')) + '\n')


def log_event(kind, gid, k, report, nbytes, run=None):
    t = report.get('timings') or report.get('Timing') or {}
    event = {'t': now_ms(), 'run': run, 'kind': kind, 'galaxyId': gid, 'k': k, 'bytes': nbytes,
             'totalMs': t.get('totalMs') if isinstance(t, dict) else None,
             'ok': not report.get('sofiaFailed') if kind == 'bootstrap' else bool(report.get('FITAchieved'))}
    with log_lock:
        append_line(os.path.join(STORE, 'events.ndjson'), event)


def now_ms():
    return int(time.time() * 1000)


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
        entry = galaxies[gid] = {'bootstraps': 0}
        fit_path = os.path.join(root, gid, 'initial-fit.json')
        if os.path.isfile(fit_path):
            with open(fit_path) as f:
                fit = json.load(f)
            entry['fit'] = {k: first(fit.get(k)) for k in FIT_SUMMARY}
            entry['fit']['at'] = datetime.fromtimestamp(os.path.getmtime(fit_path), timezone.utc).isoformat(timespec='seconds')
        entry['_seen'] = set()
        bs_path = os.path.join(root, gid, 'bootstrap.ndjson')
        if os.path.isfile(bs_path):
            with open(bs_path) as f:
                for line in f:
                    try:
                        entry['_seen'].add(int(json.loads(line)['realizationIndex']))
                    except (ValueError, KeyError, TypeError):
                        pass  # a line cut short by a crash
        entry['bootstraps'] = len(entry['_seen'])
    index['galaxies'] = galaxies


def index_writer():
    while True:
        dirty.wait()
        time.sleep(INDEX_EVERY)
        dirty.clear()
        with lock:
            snapshot = {'galaxies': {gid: {k: v for k, v in e.items() if k != '_seen'} for gid, e in index['galaxies'].items()}}
        snapshot['updatedAt'] = now()
        write_json(os.path.join(STORE, 'index.json'), snapshot)


def store_fit(report, nbytes, run=None):
    gid = report.get('galaxyId')
    if not isinstance(gid, str) or not GALAXY_RE.match(gid):
        raise ValueError('bad galaxyId')
    folder = os.path.join(STORE, 'galaxies', gid)
    bs_payload = report.pop('bootstrapPayload', None)
    if bs_payload is not None:
        write_json(os.path.join(folder, 'bootstrap-payload.json.incoming'), bs_payload)
        del bs_payload
    write_json(os.path.join(folder, 'initial-fit.json.incoming'), report)
    with lock:
        shutil.rmtree(os.path.join(folder, 'bootstrap'), ignore_errors=True)
        try:
            os.remove(os.path.join(folder, 'bootstrap.ndjson'))
        except FileNotFoundError:
            pass
        os.replace(os.path.join(folder, 'initial-fit.json.incoming'), os.path.join(folder, 'initial-fit.json'))
        has_payload = os.path.exists(os.path.join(folder, 'bootstrap-payload.json.incoming'))
        if has_payload:
            os.replace(os.path.join(folder, 'bootstrap-payload.json.incoming'), os.path.join(folder, 'bootstrap-payload.json'))
        else:
            try:
                os.remove(os.path.join(folder, 'bootstrap-payload.json'))
            except FileNotFoundError:
                pass
        summary = {k: first(report.get(k)) for k in FIT_SUMMARY}
        summary['at'] = now()
        index['galaxies'][gid] = {'fit': summary, 'bootstraps': 0, '_seen': set()}
    dirty.set()
    log_event('fit', gid, None, report, nbytes, run)
    ack = {'stored': True, 'kind': 'fit', 'galaxyId': gid, 'FITAchieved': bool(report.get('FITAchieved')),
           'error': str(report.get('error') or report.get('bootstrapPayloadError') or '')[:300] or None,
           'url': f'galaxies/{gid}/initial-fit.json',
           'payloadUrl': f'galaxies/{gid}/bootstrap-payload.json' if has_payload else None}
    ack.update({k: first(report.get(k)) for k in ('INCLINATION', 'POSITIONANGLE', 'VSYS')})
    return ack


def store_bootstrap(gid, report, nbytes, run=None):
    if not isinstance(gid, str) or not GALAXY_RE.match(gid):
        raise ValueError('bad galaxyId')
    n = report.get('realizationIndex')
    if not isinstance(n, int) or not 0 <= n < 100000:
        raise ValueError('bad realizationIndex')
    with lock:
        append_line(os.path.join(STORE, 'galaxies', gid, 'bootstrap.ndjson'), report)
        entry = index['galaxies'].setdefault(gid, {'bootstraps': 0, '_seen': set()})
        entry.setdefault('_seen', set()).add(n)
        entry['bootstraps'] = len(entry['_seen'])
    dirty.set()
    log_event('bootstrap', gid, n, report, nbytes, run)
    # The reply is the slice result, so it carries what the page needs to draw the dot.
    return {'stored': True, 'kind': 'bootstrap', 'galaxyId': gid, 'realizationIndex': n,
            'sofiaFailed': bool(report.get('sofiaFailed')), 'converged': report.get('converged'),
            'error': str(report.get('error') or '')[:300] or None,
            'url': f'galaxies/{gid}/bootstrap.ndjson'}


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
            supplied = self.path.split('?')[0].rstrip('/').rsplit('/', 1)[-1]
            if not token_ok(supplied):
                return self.reply(404, {'error': 'not found'})
            length = int(self.headers.get('Content-Length') or 0)
            if length <= 0 or length > MAX_BODY:
                self.log_message('rejected: body of %d bytes', length)
                return self.reply(413, {'error': 'missing or oversized body'})
            if disk_full():
                self.log_message('rejected: under %d GB free', MIN_FREE_GB)
                return self.reply(507, {'error': 'results disk full'})
            with busy:
                self.process(length)
        except Exception as e:  # never echo request content back
            self.log_message('rejected: %s', type(e).__name__)
            self.reply(400, {'error': type(e).__name__})

    def process(self, length):
        try:
            form = parse_qs(self.rfile.read(length).decode('ascii'), max_num_fields=32)
            kind = form.get('kind', [''])[0]
            if 'content' not in form:
                return self.reply(400, {'error': 'slice failed', 'detail': form.get('error', [''])[0][:300]})
            report = json.loads(form['content'][0])
            form_gid = form.get('galaxyId', [''])[0]
            form_run = form.get('run', [''])[0][:40] or None
            del form
            if kind == 'fit':
                ack = store_fit(report, length, form_run)
            elif kind == 'bootstrap':
                ack = store_bootstrap(form_gid, report, length, form_run)
            else:
                return self.reply(400, {'error': 'unknown kind'})
            self.reply(200, ack)
        except Exception as e:  # never echo request content back
            self.log_message('rejected: %s', type(e).__name__)
            self.reply(400, {'error': type(e).__name__})


if __name__ == '__main__':
    os.makedirs(os.path.join(STORE, 'galaxies'), exist_ok=True)
    os.chmod(STORE, 0o755)
    load_index()
    write_json(os.path.join(STORE, 'index.json'), {'galaxies': {g: {k: v for k, v in e.items() if k != '_seen'} for g, e in index['galaxies'].items()}, 'updatedAt': now()})
    threading.Thread(target=index_writer, daemon=True).start()
    ThreadingHTTPServer(('127.0.0.1', PORT), Handler).serve_forever()
