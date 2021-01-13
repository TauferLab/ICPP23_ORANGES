import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression

import argparse

def main( filename):

    gdvs = []

    with open(filename, "r+") as file:

        count_line = 0;
        for line in file:
            gdv = [];
            count_word = 0;
            for word in line.split():
                gdv.append(int(word));
                count_word += 1;
            gdvs.append(gdv)

    print(gdvs);

    avg_gdvs = np.mean(gdvs, axis = 0);
    print(avg_gdvs)
    print(np.shape(avg_gdvs))

    x_vals = np.array(range(1,len(avg_gdvs)+1));
    print(np.shape(x_vals))
    plt.bar(x_vals, avg_gdvs)

    plt.xlabel('Orbit Number')
    plt.ylabel('Average Graphlet Degree')
    plt.axis([0, 23, 0, 20]);
    plt.title('Average Graphlet Degree for each Orbit');
    plt.savefig( "avg_gdv_bar.png", bbox_inches="tight", pad_inches=0.25);

if __name__ == "__main__":
    desc = "Generates figure showing runtime data"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("filename",
                        help="Path to file storing runtime data")
    args = parser.parse_args()
    main( args.filename )






