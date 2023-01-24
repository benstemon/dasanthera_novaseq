import os
import argparse

#BE CAREFUL with this script. It deletes files permanently. You may not want to use it.

parser = argparse.ArgumentParser(description='Delete files with less than a specified number of ">" symbols')
parser.add_argument('directory', type=str, help='The directory to search')
parser.add_argument('min_gt', type=int, help='Minimum number of ">" symbols')
args = parser.parse_args()

path = args.directory
min_gt = args.min_gt

for root, dirs, files in os.walk(path):
    for file in files:
        file_path = os.path.join(root, file)
        with open(file_path, 'r') as f:
            contents = f.read()
        if contents.count(">") < min_gt:
            os.remove(file_path)
