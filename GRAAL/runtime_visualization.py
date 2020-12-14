import matplotlib.pyplot as plt
import numpy as np
import argparse

def main( filename):

    n_procs = [];
    runtimes = [];
    
    with open(filename, "r+") as file:

        count_name = 0;
        for line in file:

            if (count_name != 0):
                count = 0;
                for word in line.split():
                    print(word);
                    if (count_name == 1 and count == 2):
                        graph_name = word;
                    if (count == 1):
                        n_procs.append(int(word));
                    if (count == 3 and count_name == 1):
                        graph_size = int(word);
                    if (count == 4):
                        runtimes.append(float(word));
                    count += 1;
            count_name += 1

    # Perform analysis to average runs
    # Start printing array forms of stored lists above with mpl
    print(n_procs);
    print(runtimes);
    plt.plot(np.array(n_procs), np.array(runtimes), 'bo');
    plt.xlabel("Number of Processes");
    plt.ylabel("Runtime");
    plt.xlim(0, 4);
    plt.ylim(0, 1);
    plt.savefig( "fido_runtime_per_nprocs.png",
                 bbox_inches="tight",
                 pad_inches=0.25
               );

if __name__ == "__main__":
    desc = "Generates figure showing runtime data"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("filename", 
                        help="Path to file storing runtime data")
    args = parser.parse_args()
    main( args.filename )
