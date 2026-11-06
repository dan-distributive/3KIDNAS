#!/usr/bin/env python3
"""Add ?v=<content hash> to catalogue.html's local script tags, so browsers fetch changed files
instead of reusing cached copies (a normal reload does not revalidate scripts). Run before
copying catalogue.html to the server."""
import hashlib, os, re, sys
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
page = os.path.join(root, 'catalogue.html')
html = open(page).read()
def stamp(m):
    path = m.group(1)
    digest = hashlib.sha256(open(os.path.join(root, path), 'rb').read()).hexdigest()[:10]
    return f'<script src="{path}?v={digest}"></script>'
html = re.sub(r'<script src="((?!https?:)[^"?]+)(?:\?v=[0-9a-f]+)?"></script>', stamp, html)
open(page, 'w').write(html)
print('\n'.join(re.findall(r'<script src="((?!https?:)[^"]+)"', html)))
