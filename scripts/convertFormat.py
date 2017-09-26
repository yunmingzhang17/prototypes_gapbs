import sys

# converts float weights to integer weights by flooring (GAPBS only supports integer weights for binary format for now)
def convert_float_to_int_weight(filename):
    out_file = open(filename + "_int_wt", 'w')
    with open(filename) as f:
        for line in f:
            line_items = line.split()
            if len(line_items) > 1:
                #print line_items
                src = int(line_items[0])
                dst = int(line_items[1])
                wt = int(float(line_items[2]))
                out_line = str(src) + " " + str(dst) + " " + str(wt) + "\n"
                out_file.write(out_line)


def find_largest_node(filename):
    max_src = 0
    max_dst = 0
    with open(filename) as f:
        for line in f:
            line_items = line.split()

            if len(line_items) > 1:
                src = int(line_items[0])
                dst = int(line_items[1])
                if src > max_src:
                    max_src = src
                if dst > max_dst:
                    max_dst = dst
    print (" max src: " + str(max_src) + "\n")
    print (" max dst: " + str(max_dst) + "\n")

filename = sys.argv[1]
#convert_float_to_int_weight(filename)
find_largest_node(filename)
