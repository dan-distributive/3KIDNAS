#!/usr/bin/env python3
"""
build-wallaby-dataset.py -- turn WALLABY source releases into a web dataset for catalogue.html.

Reads every <release>.fits catalogue and its <release>_sources.tar in SRC_DIR and writes OUT_DIR:

  catalogue.json        one row per galaxy: catalogue columns plus beam, cube shape, systemic
                        velocity and the PA / inclination estimates 3KIDNAS starts the fit from
  galaxies/<name>.json  { galaxyId, observedCubeRawFitsB64, observedMaskRawFitsB64 }, the slice
                        a worker fetches by URL for runInitialFit

The cubes are released on a frequency axis; runInitialFit only takes velocity axes. Each cube's
spectral header cards are rewritten to optical velocity (VOPT, m/s), linearised at channel
NAXIS3 // 2 with c = 2.9979245e8, which reproduces the _VelCube.fits files in 3KIDNASTests to
1e-12. The data bytes are not touched. Masks are passed through as they are.

PA and inclination follow GeometryEstimates.py with IncMethod 'WALLABY_Like' and the catalogue's
kin_pa, as SoFiA_Driver.py does. A galaxy in several releases is taken from the latest one.
Cubes without BMAJ are skipped unless --assume-beam-arcsec is given.

Python 3 standard library only.
  python3 build-wallaby-dataset.py SRC_DIR OUT_DIR [--assume-beam-arcsec 30] [--include REGEX]
"""
import argparse
import base64
import glob
import hashlib
import json
import math
import os
import re
import struct
import sys
import tarfile
from datetime import datetime, timezone

C_MS = 2.9979245e8
HI_REST_HZ = 1.42040575179e9


def read_header(data, start=0):
    """FITS header cards from `start` up to END, and the offset of the data that follows."""
    cards, pos = [], start
    while True:
        block = data[pos:pos + 2880].decode('ascii', 'replace')
        pos += 2880
        for i in range(0, 2880, 80):
            cards.append(block[i:i + 80])
            if block[i:i + 8] == 'END     ':
                return cards, pos


def card_value(card):
    raw = card[10:]
    if raw.lstrip().startswith("'"):
        return raw.lstrip()[1:].split("'")[0].strip()
    raw = raw.split('/')[0].strip()
    try:
        return int(raw)
    except ValueError:
        try:
            return float(raw)
        except ValueError:
            return raw


def header_dict(cards):
    return {c[:8].strip(): card_value(c) for c in cards if len(c) > 9 and c[8] == '='}


def format_card(key, value):
    if isinstance(value, str):
        quoted = "'" + value.ljust(8) + "'"
        text = f'{key:<8}= {quoted:<20}'
    else:
        text = f'{key:<8}= {repr(float(value)) if isinstance(value, float) else str(value):>20}'
    return text.ljust(80)[:80]


def set_card(cards, key, value):
    for i, c in enumerate(cards):
        if c[:8].strip() == key:
            cards[i] = format_card(key, value)
            return
    end = next(i for i, c in enumerate(cards) if c[:8] == 'END     ')
    cards.insert(end, format_card(key, value))


def rebuild(cards, data, data_start):
    text = ''.join(cards)
    text += ' ' * (-len(text) % 2880)
    return text.encode('ascii') + data[data_start:]


def to_velocity_cube(data, assume_beam_deg):
    """Rewrite a FREQ cube's spectral axis to VOPT; returns (new bytes, header dict)."""
    cards, data_start = read_header(data)
    h = header_dict(cards)
    if str(h.get('CTYPE3', '')).strip() != 'FREQ':
        raise ValueError(f"CTYPE3 is {h.get('CTYPE3')!r}, expected FREQ")
    f0 = float(h.get('RESTFRQ') or h.get('RESTFREQ') or HI_REST_HZ)
    ref = h['NAXIS3'] // 2
    freq = lambda p: h['CRVAL3'] + (p - h['CRPIX3']) * h['CDELT3']
    vel = lambda p: C_MS * (f0 / freq(p) - 1)
    set_card(cards, 'CTYPE3', 'VOPT')
    set_card(cards, 'CRPIX3', ref)
    set_card(cards, 'CRVAL3', vel(ref))
    set_card(cards, 'CDELT3', vel(ref + 1) - vel(ref))
    set_card(cards, 'CUNIT3', 'm/s')
    if 'BMAJ' not in h:
        if not assume_beam_deg:
            raise ValueError('no BMAJ in header')
        set_card(cards, 'BMAJ', assume_beam_deg)
        set_card(cards, 'BMIN', assume_beam_deg)
        set_card(cards, 'BPA', 0.0)
    out = rebuild(cards, data, data_start)
    return out, header_dict(read_header(out)[0])


def read_table(path):
    """Rows of the first binary-table extension of a FITS catalogue."""
    data = open(path, 'rb').read()
    _, pos = read_header(data)
    cards, pos = read_header(data, pos)
    h = header_dict(cards)
    width, n, nf = h['NAXIS1'], h['NAXIS2'], h['TFIELDS']
    sizes = {'A': 1, 'B': 1, 'L': 1, 'I': 2, 'J': 4, 'E': 4, 'K': 8, 'D': 8}
    codes = {'B': 'B', 'L': '?', 'I': 'h', 'J': 'i', 'E': 'f', 'K': 'q', 'D': 'd'}
    fields, off = [], 0
    for j in range(1, nf + 1):
        form = h[f'TFORM{j}']
        count, kind = int(form[:-1] or 1), form[-1]
        fields.append((h[f'TTYPE{j}'], off, count, kind))
        off += count * sizes[kind]
    rows = []
    for r in range(n):
        rec = data[pos + r * width: pos + (r + 1) * width]
        row = {}
        for name, o, count, kind in fields:
            raw = rec[o:o + count * sizes[kind]]
            if kind == 'A':
                row[name] = raw.decode('ascii', 'replace').strip('\x00 ')
            else:
                row[name] = struct.unpack('>' + codes[kind], raw[:sizes[kind]])[0]
        rows.append(row)
    return rows


def geometry_estimate(row, beam_pix):
    """GeometryEstimates.GetGeometryEstimates with IncMethod 'WALLABY_Like'."""
    min2 = (2 * row['ell_min']) ** 2 - beam_pix ** 2
    maj2 = (2 * row['ell_maj']) ** 2 - beam_pix ** 2
    if maj2 < 0:
        maj2 = beam_pix ** 2
    inc = 89.0 if min2 <= 0 else math.degrees(math.acos(math.sqrt(min2) / math.sqrt(maj2)))
    pa = row['ell_pa'] if math.isnan(row['kin_pa']) else row['kin_pa']
    return pa, inc


def finite(x):
    return None if isinstance(x, float) and not math.isfinite(x) else x


def release_date(name):
    m = re.search(r'(\d{4}-\d{2}-\d{2})', name)
    return m.group(1) if m else ''


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('src')
    ap.add_argument('out')
    ap.add_argument('--assume-beam-arcsec', type=float, default=None)
    ap.add_argument('--include', default='.', help='regex on release name')
    args = ap.parse_args()
    assume_beam_deg = args.assume_beam_arcsec / 3600 if args.assume_beam_arcsec else None
    os.makedirs(os.path.join(args.out, 'galaxies'), exist_ok=True)

    releases = []
    for cat in sorted(glob.glob(os.path.join(args.src, '*.fits'))):
        name = os.path.basename(cat)[:-5]
        tar = os.path.join(args.src, name + '_sources.tar')
        if ' (' in name or not os.path.exists(tar) or not re.search(args.include, name):
            continue
        releases.append((release_date(name), name, cat, tar))
    releases.sort(reverse=True)  # latest first, so it wins on duplicates

    catalogue, seen, skipped = [], set(), {}
    for _, release, cat, tar in releases:
        rows = {r['name'].replace(' ', '_'): r for r in read_table(cat)}
        field = re.sub(r'_(DR\d+_)?\d{4}-\d{2}-\d{2}$', '', release).replace('_', ' ')
        pending = {}
        with tarfile.open(tar) as tf:
            for member in tf:
                m = re.search(r'/(WALLABY_J[0-9+-]+)_(cube|mask)\.fits$', member.name)
                if not m or not member.isfile():
                    continue
                gid, kind = m.groups()
                if gid in seen or gid not in rows:
                    continue
                pending.setdefault(gid, {})[kind] = tf.extractfile(member).read()
                if len(pending[gid]) < 2:
                    continue
                files = pending.pop(gid)
                try:
                    cube, h = to_velocity_cube(files['cube'], assume_beam_deg)
                except ValueError as e:
                    skipped[str(e)] = skipped.get(str(e), 0) + 1
                    continue
                row = rows[gid]
                beam_pix = abs(h['BMAJ'] / h['CDELT1'])
                pa, inc = geometry_estimate(row, beam_pix)
                body = json.dumps({
                    'galaxyId': gid,
                    'observedCubeRawFitsB64': base64.b64encode(cube).decode(),
                    'observedMaskRawFitsB64': base64.b64encode(files['mask']).decode(),
                }, separators=(',', ':')).encode()
                with open(os.path.join(args.out, 'galaxies', gid + '.json'), 'wb') as f:
                    f.write(body)
                seen.add(gid)
                catalogue.append({
                    'id': gid, 'name': row['name'], 'field': field, 'release': release,
                    'ra': finite(row['ra']), 'dec': finite(row['dec']),
                    'vsys': round(C_MS * (HI_REST_HZ / row['freq'] - 1) / 1000, 1),
                    'fSum': finite(row['f_sum']), 'rms': finite(row['rms']),
                    'ellMaj': finite(row['ell_maj']), 'ellMin': finite(row['ell_min']),
                    'ellPa': finite(row['ell_pa']), 'kinPa': finite(row['kin_pa']),
                    'flag': int(row['flag']),
                    'beamPix': round(beam_pix, 3),
                    'sizeBeams': round(2 * row['ell_maj'] / beam_pix, 2),
                    'beamAssumed': bool(assume_beam_deg) and 'BMAJ' not in header_dict(read_header(files['cube'])[0]),
                    'shape': [h['NAXIS1'], h['NAXIS2'], h['NAXIS3']],
                    'paEst': round(pa % 360, 2), 'incEst': round(inc, 2),
                    'file': 'galaxies/' + gid + '.json', 'bytes': len(body),
                    'sha256': hashlib.sha256(body).hexdigest(),
                })
        print(f'{release}: {sum(1 for c in catalogue if c["release"] == release)} galaxies', flush=True)

    catalogue.sort(key=lambda c: (c['field'], c['ra'] or 0))
    meta = {
        'created': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'releases': [r[1] for r in releases],
        'skipped': skipped,
        'spectralAxis': 'VOPT m/s, linearised at NAXIS3//2, c=2.9979245e8 (matches 3KIDNASTests _VelCube.fits)',
        'geometry': "GeometryEstimates WALLABY_Like, PA = kin_pa",
        'galaxies': catalogue,
    }
    tmp = os.path.join(args.out, 'catalogue.json.tmp')
    with open(tmp, 'w') as f:
        json.dump(meta, f, separators=(',', ':'))
    os.replace(tmp, os.path.join(args.out, 'catalogue.json'))
    print(f'{len(catalogue)} galaxies, skipped {skipped}')


if __name__ == '__main__':
    sys.exit(main())
