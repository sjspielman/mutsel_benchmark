## Script to replace placeholder with correct working directory.
## Using python rather than sed because I want to replace placeholder w/ a variable like $SCRATCH or $HOME, and those contain slashes :(


import sys

assert(len(sys.argv) == 3), "Usage: python replace_placeholder.py <batchfile> <replacement_string>."

batchfile = sys.argv[1]
rep_string = sys.argv[2]

with open(batchfile, "r") as f:
    file = f.read()  

file = file.replace("placeholder", rep_string)

with open(batchfile, "w") as f:
    f.write(file)