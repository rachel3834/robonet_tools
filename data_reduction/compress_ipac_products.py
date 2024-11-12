from os import path, walk
import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('data_dir', help='Path to IPAC lightcurves top-level directory')
args = parser.parse_args()

nlc = 0
for (root,dirs,files) in walk(args.data_dir):
    for f in files:
        if '.gz' not in f:
            file_path = path.join(root, f)
            result = subprocess.run(["gzip", file_path], capture_output=True, text=True)
            if '.fits' in f: nlc += 1

print('Found ' + str(nlc) + ' lightcurves')
