import sys

# Parameters
SPACER = 10000
INPUT_FILE = sys.argv[1]
OUTPUT_FILE = INPUT_FILE.replace(".txt", "") + "_spaced.txt"

# Dictionary to track last written position per chromosome
last_positions = {}

with open(INPUT_FILE, 'r') as infile, open(OUTPUT_FILE, 'w') as outfile:
    header = infile.readline()
    outfile.write(header)  # Write header as-is

    for line in infile:
        fields = line.strip().split('\t')
        if len(fields) < 2:
            continue  # Skip malformed lines

        # Extract chromosome and position from second-to-last field
        try:
            chrom_pos = fields[-2]
            chrom, pos = chrom_pos.rsplit("_", 1)
            pos = int(pos)
        except ValueError:
            continue  # Skip lines that don't have proper format

        # Check if SNP should be written based on spacing
        if chrom not in last_positions or pos - last_positions[chrom] >= SPACER:
            outfile.write(line)
            last_positions[chrom] = pos
        elif pos < last_positions[chrom]:  # Optional: update if out of order
            last_positions[chrom] = pos
            outfile.write(line)
