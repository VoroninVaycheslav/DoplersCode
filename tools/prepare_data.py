#!/usr/bin/env python3
"""Download a pinned ENDFtk evaluation and extract unchanged ENDF sections."""
import argparse
import hashlib
import json
import pathlib
import urllib.request

COMMIT = '6596073f17b21952407f6f16b5b812ca3b4a8d67'
URL = f'https://raw.githubusercontent.com/njoy/ENDFtk/{COMMIT}/python/examples/resources/n-U238.endf'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=pathlib.Path, help='Use an already downloaded evaluation')
    parser.add_argument('--output', type=pathlib.Path, default=pathlib.Path('data/u238'))
    args = parser.parse_args()
    manifest = pathlib.Path(__file__).resolve().parents[1] / 'data/u238/provenance.json'
    expected = json.loads(manifest.read_text())['source_sha256']
    raw = args.source.read_bytes() if args.source else urllib.request.urlopen(URL, timeout=120).read()
    actual = hashlib.sha256(raw).hexdigest()
    if actual != expected:
        raise SystemExit(f'Source checksum mismatch: {actual}')
    sections = {(1, 451): [], (3, 1): [], (3, 2): [], (4, 2): []}
    for line in raw.splitlines(keepends=True):
        if len(line) < 75:
            continue
        key = (int(line[70:72]), int(line[72:75]))
        if key in sections:
            sections[key].append(line)
    names = {(1, 451): 'evaluation_header.endf', (3, 1): 'total.endf',
             (3, 2): 'elastic.endf', (4, 2): 'angular.endf'}
    args.output.mkdir(parents=True, exist_ok=True)
    for key, lines in sections.items():
        if not lines:
            raise SystemExit(f'Missing section {key}')
        content = b''.join(lines)
        path = args.output / names[key]
        path.write_bytes(content)
        print(path, len(content), hashlib.sha256(content).hexdigest())


if __name__ == '__main__':
    main()
