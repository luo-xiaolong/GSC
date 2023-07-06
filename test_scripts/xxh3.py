import sys
# sys.path.append('/home/pdx/.local/lib/python3.8/site-packages/')

import xxhash

import glob
import os


def cal_hash_from_file(fname):
    if not os.path.exists(fname):
        return False
    x = xxhash.xxh3_64()
    with open(fname, 'rb') as F:
        for chunk in iter(lambda: F.read(8192), b''):
            x.update(chunk)
    return x.hexdigest()


def cal_hash_from_stdin():
    x = xxhash.xxh3_64()
    while True:
        chunk = sys.stdin.buffer.read(8192)
        if not chunk:
            break
        x.update(chunk)
    return x.hexdigest()


if len(sys.argv) > 1:
    for i in sys.argv[1:]:
        for f in glob.iglob(i):
            print('%s\t%s' % (i, cal_hash_from_file(i)))
else:
    print(cal_hash_from_stdin())
