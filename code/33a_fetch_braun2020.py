#!/usr/bin/env python3
"""Download Braun et al. Nat Med 2020 Supplementary Tables (MOESM2 xlsx) with parallel byte ranges.

Keeps an existing contiguous prefix and fills the rest. Verifies the final size against Content-Range.
"""
import os
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

URL = ("https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-020-0839-y/"
       "MediaObjects/41591_2020_839_MOESM2_ESM.xlsx")
P2 = "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
DEST = os.path.join(P2, "data/external_validation/braun2020/41591_2020_839_MOESM2_ESM.xlsx")
CHUNK = 2 * 1024 * 1024
WORKERS = 8
RETRIES = 5


def fetch(start, end):
    req = urllib.request.Request(URL, headers={"Range": f"bytes={start}-{end}"})
    with urllib.request.urlopen(req, timeout=120) as r:
        data = r.read()
    if len(data) != end - start + 1:
        raise IOError(f"short read {start}-{end}: {len(data)}")
    return start, data


def total_size():
    req = urllib.request.Request(URL, headers={"Range": "bytes=0-0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return int(r.headers["Content-Range"].split("/")[-1])


total = total_size()
have = os.path.getsize(DEST) if os.path.exists(DEST) else 0
have = min(have, total)
have -= have % CHUNK
print("total", total, "keep prefix", have)

with open(DEST, "r+b" if os.path.exists(DEST) else "wb") as f:
    f.truncate(total)

ranges = [(s, min(s + CHUNK, total) - 1) for s in range(have, total, CHUNK)]
pending = ranges
for attempt in range(RETRIES):
    failed = []
    with ThreadPoolExecutor(WORKERS) as ex, open(DEST, "r+b") as f:
        futs = {ex.submit(fetch, s, e): (s, e) for s, e in pending}
        for i, fu in enumerate(as_completed(futs), 1):
            try:
                s, data = fu.result()
                f.seek(s)
                f.write(data)
            except Exception as exc:  # noqa: BLE001
                failed.append(futs[fu])
                print("fail", futs[fu], exc, file=sys.stderr)
            if i % 10 == 0:
                print(f"attempt {attempt + 1}: {i}/{len(pending)}", flush=True)
    if not failed:
        break
    pending = failed

if pending and failed:
    sys.exit(f"incomplete: {len(failed)} ranges failed")
assert os.path.getsize(DEST) == total
print("ALL DONE", total)
