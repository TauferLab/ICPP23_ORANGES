import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression

import argparse

def main( filename):

    n_procs = [];
    runtimes = [];
    
    with open(filename, "r+") as file:

        count_name = 0;
        for line in file:

            if (count_name != 0):
                count = 0;
                use_line = False;
                for word in line.split():
                    print(word);
                    if ( (count_name == 1) and (count == 2) ):
                        graph_name1 = word;
                    if ( (count_name == 1) and (count == 3) ):
                        graph_name2 = word;
                    if (count == 1):
                        if (int(word) < 17):
                            n_procs.append(int(word));
                            use_line = True;
                    if ( (count == 4) and (count_name == 1) ):
                        graph_size = int(word);
                    if ( (count == 5) and (use_line == True) ):
                        runtimes.append(float(word));
                    count += 1;
            count_name += 1

    # Perform analysis to average runs
    # Start printing array forms of stored lists above with mpl
    print(n_procs);
    print(runtimes);

    # Prep variables for plotting
    x = np.array(n_procs);
    y = np.array(runtimes);
    p = np.poly1d(np.polyfit(x, y, 1))

    print(p);
    # Create figure and plot data and linear regression to it
    plt.figure(figsize=(4, 3))
    ax = plt.axes()
    ax.scatter(x, y, c='b', marker='o')
    ax.plot(x, p(x))
    plt.xlabel("Number of Processes");
    plt.ylabel("Runtime");
    ax.axis('tight');
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
