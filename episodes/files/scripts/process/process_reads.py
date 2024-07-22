#!/usr/bin/env python
import gzip
import sys
reads = 0
bases = 0

# Read gzipped fastq file
with gzip.open(sys.argv[1], 'rb') as read:
    for id in read:
      seq = next(read)
      reads += 1
      bases += len(seq.strip())
      next(read)
      next(read)

print("reads", reads)
print("bases", bases)