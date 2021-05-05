import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression

import argparse

def main( filepath , run_max ):

    runtimes_prior_gather_2d = []*int(run_max);
    runtimes_post_gather_2d = []*int(run_max);
    runtimes_total_2d = []*int(run_max);

    for run_num in range(1, int(run_max)+1):
        if (run_num < 10):
            filename = filepath + "/run_00" + str(run_num) + "/runtimes_rec.txt";
        if (run_num >= 10 and run_num < 100):
            filename = filepath + "/run_0" + str(run_num) + "/runtimes_rec.txt";
        with open(filename, "r") as file:

            runtimes_prior_gather = [];                                                                                                                                                                                              
            runtimes_post_gather = [];
            runtimes_total = [];

            count_line = 0;
            for line in file:
                if (count_line != 0):
                    count_word = 0;
                    #runtimes_prior_gather = [];
                    #runtimes_post_gather = [];
                    #runtimes_total = [];
                    for word in line.split():
                        if (count_word == 1):
                            runtimes_prior_gather.append(float(word));
                        if (count_word == 2):
                            runtimes_post_gather.append(float(word));
                        if (count_word == 3):
                            runtimes_total.append(float(word));
                        count_word += 1;
                count_line += 1;
            runtimes_prior_gather_2d.append(runtimes_prior_gather);
            runtimes_post_gather_2d.append(runtimes_post_gather);
            runtimes_total_2d.append(runtimes_total);

    #print(runtimes_total);
    #print(runtimes_prior_gather);
    #print(runtimes_post_gather);

    #x_vals = np.array(range(0, len(runtimes_total)));
    #print(np.shape(x_vals))
    prior_gather = np.array(runtimes_prior_gather_2d);
    post_gather = np.array(runtimes_post_gather_2d);
    total_time = np.array(runtimes_total_2d);

    print(np.shape(total_time));
    ax = 0;
    prior_gather_avg = np.mean(prior_gather, axis=ax);
    post_gather_avg = np.mean(post_gather, axis=ax);
    total_time_avg = np.mean(total_time, axis=ax);
    x_vals = np.array(range(0, len(runtimes_total)));
    print(total_time_avg);
    gather_addition = np.subtract(post_gather_avg, prior_gather_avg);
    print(gather_addition)
    print(x_vals)

    fig1 = plt.figure();
    plt.bar(x_vals, prior_gather_avg, color='b');
    plt.bar(x_vals, gather_addition, bottom=prior_gather_avg, color='r');
    plt.xticks(x_vals);

    plt.xlabel('MPI Rank x')
    plt.ylabel('Runtime (s.) on MPI Rank x')
    plt.legend(['Prior Gather', 'Post Gather']);
    plt.axis('tight')
    plt.title('Time Taken on Each Rank of GDV Vector Calculation');
    plt.savefig( "png_files/load_imb_bar_vec.png", bbox_inches="tight", pad_inches=0.25);


    fig2 = plt.figure();
    plt.bar(x_vals, prior_gather_avg, color='b');
    plt.xticks(x_vals);

    plt.xlabel('MPI Rank x')
    plt.ylabel('Runtime (s.) on MPI Rank x')
    #plt.legend(['Prior Gather', 'Post Gather']);
    plt.axis('tight')
    plt.title('Time Taken on Each Rank of GDV Vector Calculation Prior to MPI_Gather');
    plt.savefig( "png_files/load_imb_bar_prior_gather.png", bbox_inches="tight", pad_inches=0.25);


    fig3 = plt.figure();
    plt.bar(x_vals, post_gather_avg, color='b');
    plt.xticks(x_vals);

    plt.xlabel('MPI Rank x')
    plt.ylabel('Runtime (s.) on MPI Rank x')
    #plt.legend(['Prior Gather', 'Post Gather']);
    plt.axis('tight')
    plt.title('Time Taken on Each Rank of GDV Vector Calculation');
    plt.savefig( "png_files/load_imb_bar_post_gather.png", bbox_inches="tight", pad_inches=0.25);


    fig4 = plt.figure();
    plt.bar(x_vals, total_time_avg, color='b');
    #plt.bar(x_vals, post_gather, bottom=prior_gather);
    plt.xticks(x_vals)

    plt.xlabel('MPI Rank x')
    plt.ylabel('Runtime (s.) on MPI Rank x')
    plt.axis('tight')
    plt.title('Time Taken on Each Rank of Full Similarity Calculation');
    plt.savefig( "png_files/load_imb_bar_full.png", bbox_inches="tight", pad_inches=0.25);


if __name__ == "__main__":
    desc = "Generates figure showing runtime load imbalance data"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("filepath",
                        help="Path to run directories storing runtime data")
    parser.add_argument("-r", "--run_max", required=False, default=1,
                        help="Number of run directories to use")
    args = parser.parse_args()
    main( args.filepath, args.run_max )






