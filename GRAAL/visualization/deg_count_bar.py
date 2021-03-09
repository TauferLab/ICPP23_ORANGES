import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression

import argparse

def main( filename , procs ):

    deg_count_file_data = []*int(procs);

    #filename = filepath + "/run_00" + str(run_num) + "/runtimes_rec.txt"
    with open(filename, "r") as file:
        
        count_line = 0;
        for line in file:
            file_data_line = [];
            count_word = 0;
            for word in line.split():
                if (count_word == 0):
                    file_data_line.append(int(word));
                if (count_word == 1):
                    file_data_line.append(float(word));
                if (count_word == 2):
                    file_data_line.append(int(word));
                count_word += 1;
            count_line += 1;
            deg_count_file_data.append(file_data_line);

    #print(runtimes_total);
    #print(runtimes_prior_gather);
    #print(runtimes_post_gather);

    avg_deg = [];
    num_degs = [];
    print(len(deg_count_file_data));
    for i in range(len(deg_count_file_data)):
        for j in range(len(deg_count_file_data)):
          if (i == deg_count_file_data[j][0]):
              avg_deg.append(deg_count_file_data[j][1]);
              num_degs.append(deg_count_file_data[j][2]);

    print(avg_deg)
    print(num_degs)
    x_vals = np.array(range(0, len(deg_count_file_data)));
    degs = np.array(avg_deg);
    num_per_proc = np.array(num_degs);


    fig2 = plt.figure();
    bar1 = plt.bar(x_vals, degs);
    plt.xticks(x_vals);

    proc_count = 0;
    for box in bar1:
        height = box.get_height();
        num_deg = num_degs[proc_count];
        #plt.text(box.get_x() + box.get_width()/2.0, height, '%d' % int(num_deg), ha='center', va='bottom');
        proc_count += 1;

    plt.xlabel('MPI Rank x')
    plt.ylabel('Average Node Degree Assigned to x')
    #plt.legend(['Prior Gather', 'Post Gather']);
    plt.axis('tight')
    plt.title('Amount of work assigned to each MPI process based on average node degree');
    plt.savefig( "png_files/avg_deg_bar.png", bbox_inches="tight", pad_inches=0.25);


if __name__ == "__main__":
    desc = "Generates figure showing average node degree assigned to each MPI process"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("filename",
                        help="Path to run directories storing average degree count data")
    parser.add_argument("-p", "--procs", required=True,
                        help="Number of processes used in the run of fido")
    args = parser.parse_args()
    main( args.filename, args.procs )






