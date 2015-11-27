#!/usr/local/bin/python
"""
plot_results.py

"""
import matplotlib.pyplot as plt

def parse_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--pkl-group', action="append", nargs='+')
    return parser.parse_args()

def main(args):
    return

if __name__ == "__main__":
    main(parse_args())
