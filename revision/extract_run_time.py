
import os


seqs = ['hifi', 'ont']
coverages = ['5x', '10x','20x', '40x']

input_path = "/mnt/d/workspace/svision-pro/paper_revision/run_time/pc_8_s5"
# input_path = "/mnt/d/workspace/svision-pro/paper_revision/run_time/cluster_24_s5"


for seq in seqs:
    for coverage in coverages:

        file_path = os.path.join(input_path, "{}.{}.svision-pro.cmd.log".format(seq, coverage))

        time = ""
        memory = ""

        with open(file_path) as fin:
            for line in fin:

                if "Elapsed (wall clock) time" in line:

                    time = line.strip().split(" ")[-1]

                    time = time.split(":")

                    if len(time) == 3:

                        time = int(time[0]) * 60 + int(time[1])

                    elif len(time) == 2:

                        time = int(time[0])

                if "Maximum resident set size" in line:

                    memory = round(int(line.strip().split(" ")[-1]) / 1024 / 1024, 2)


        print("{}\t{}\t{}\t{}".format(seq, coverage, time, memory))



