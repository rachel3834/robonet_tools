from pyDANDIA import crossmatch
import argparse

def rm_star_norm_tables(args, xmatch):
    if hasattr(xmatch, 'normalizations'):
        delattr(xmatch, 'normalizations')
        xmatch.save(args.xmatch_file)

def load_xmatch_table(args):
    xmatch = crossmatch.CrossMatchTable()
    xmatch.load(args.xmatch_file, log=None)
    return xmatch

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('xmatch_file', help='Path to crossmatch table file')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = get_args()
    xmatch = load_xmatch_table(args)
    rm_star_norm_tables(args, xmatch)
