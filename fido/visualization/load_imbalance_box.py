import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression

import argparse

def main( filepath , run_max, n_procs ):

    runtimes_graph_1_2d = []*int(run_max);
    runtimes_graph_2_2d = []*int(run_max);
    runtimes_total_2d = []*int(run_max);
    runtimes_x = []*3*int(run_max);
    runtimes_y = []*3*int(run_max);

    # Intake Overhead Runtimes
    for run_num in range(0, int(run_max):
        if (run_num < 10):
            filename = filepath + "/run_00" + str(run_num) + "/runtime_data/runtimes_rec_over.txt";
        if (run_num >= 10 and run_num < 100):
            filename = filepath + "/run_0" + str(run_num) + "/runtime_data/runtimes_rec_over.txt";
        with open(filename, "r") as file:

            runtimes_graph_1 = [];                                                                                                                                                                                              
            runtimes_graph_2 = [];
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
                            runtimes_graph_1.append(float(word));
                        if (count_word == 2):
                            runtimes_graph_2.append(float(word));
                        if (count_word == 3):
                            runtimes_total.append(float(word));
                        count_word += 1;
                count_line += 1;
            runtimes_graph_1_2d.append(runtimes_graph_1);
            runtimes_graph_2_2d.append(runtimes_graph_2);
            runtimes_total_2d.append(runtimes_total);

    proc_list = [1, 2, 4, 8, 16, 32];

    # Intake graph 1 computation data
    for run_num in range(0, int(run_max)):

        runtimes_x_run = []*3

        for rank_num in n_procs:
            if (run_num < 10):
                filename = filepath + "/run_00" + str(run_num) + "/runtime_data/runtimes_rec_x_" + rank_num + ".txt";
            if (run_num >= 10 and run_num < 100):
                filename = filepath + "/run_0" + str(run_num) + "/runtime_data/runtimes_rec_x_" + rank_num + ".txt";
            with open(filename, "r") as file:

                runtimes_x_rank = [];
                count_line = 0;

                for line in file:
                    count_word = 0;
                    for word in line.split():
                        if (count_word == 0):
                            runtimes_x_rank.append(int(word));
                        if (count_word == 1):
                            runtimes_x_rank.append(int(word));
                        if (count_word == 2):
                            runtimes_x_rank.append(float(word));
                        count_word += 1;
                    count_line += 1;
                runtimes_x_run.append(runtimes_x_rank);
        runtimes_x_run = sorted(runtimes_x_run, key=lambda x: x[0])
        runtimes_x.append(runtimes_x_run);

    # Intake graph 2 computation data
    for run_num in range(0, int(run_max)):

        runtimes_y_run = []*3

        for rank_num in n_procs:
            if (run_num < 10):
                filename = filepath + "/run_00" + str(run_num) + "/runtime_data/runtimes_rec_y_" + rank_num + ".txt";
            if (run_num >= 10 and run_num < 100):
                filename = filepath + "/run_0" + str(run_num) + "/runtime_data/runtimes_rec_y_" + rank_num + ".txt";
            with open(filename, "r") as file:

                runtimes_y_rank = [];
                count_line = 0;

                for line in file:
                    count_word = 0;
                    for word in line.split():
                        if (count_word == 0):
                            runtimes_y_rank.append(int(word));
                        if (count_word == 1):
                            runtimes_y_rank.append(int(word));
                        if (count_word == 2):
                            runtimes_y_rank.append(float(word));
                        count_word += 1;
                    count_line += 1;
                runtimes_y_run.append(runtimes_y_rank);
        runtimes_y_run = sorted(runtimes_y_run, key=lambda x: x[0])
        runtimes_y.append(runtimes_y_run);


    #print(runtimes_total);
    #print(runtimes_prior_gather);
    #print(runtimes_post_gather);

    #x_vals = np.array(range(0, len(runtimes_total)));
    #print(np.shape(x_vals))
    graph_1_comm = np.array(runtimes_graph_1_2d);
    graph_2_comm = np.array(runtimes_graph_2_2d);
    total_time = np.array(runtimes_total_2d);
    graph_1_comp = np.array(runtimes_x);
    graph_2_comp = np.array(runtimes_y);

    print(np.shape(total_time));
    print(np.shape(graph_1_comm));
    ax = 0;
    graph_1_comm_avg = np.mean(graph_1_comm, axis=ax);
    graph_2_comm_avg = np.mean(graph_2_comm, axis=ax);
    total_time_avg = np.mean(total_time, axis=ax);
    x_vals = np.array(range(0, len(runtimes_total)));
    print(total_time_avg);
    #gather_addition = np.subtract(post_gather_avg, prior_gather_avg);
    #print(gather_addition)
    print(x_vals)

    #fig1 = plt.figure();
    #plt.bar(x_vals, prior_gather_avg, color='b');
    #plt.bar(x_vals, gather_addition, bottom=prior_gather_avg, color='r');
    #plt.xticks(x_vals);

    #plt.xlabel('MPI Rank x')
    #plt.ylabel('Runtime (s.) on MPI Rank x')
    #plt.legend(['Prior Gather', 'Post Gather']);
    #plt.axis('tight')
    #plt.title('Time Taken on Each Rank of GDV Vector Calculation');
    #plt.savefig( "png_files/load_imb_bar_vec.png", bbox_inches="tight", pad_inches=0.25);


    c1 = "blue"
    c2 = "red"
    c3 = "purple"
    fig2 = plt.figure();
    #plt.bar(x_vals, prior_gather_avg, color='b');
    plt.boxplot(total_time[0], positions=np.array([0]), notch=True, patch_artist=True,
            boxprops=dict(facecolor=c1, color=c1),
            capprops=dict(color=c1),
            whiskerprops=dict(color=c1),
            flierprops=dict(color=c1, markeredgecolor=c1),
            medianprops=dict(color=c1));
    plt.boxplot(graph_1_comm + graph_2_comm, positions=x_vals, notch=True, patch_artist=True,
            boxprops=dict(facecolor=c2, color=c2),
            capprops=dict(color=c2),
            whiskerprops=dict(color=c2),
            flierprops=dict(color=c2, markeredgecolor=c2),
            medianprops=dict(color=c2));
    plt.boxplot(graph_1_comm, positions=x_vals, notch=True, patch_artist=True,
            boxprops=dict(facecolor=c3, color=c3),
            capprops=dict(color=c3),
            whiskerprops=dict(color=c3),
            flierprops=dict(color=c3, markeredgecolor=c3),
            medianprops=dict(color=c3));    
    plt.xticks(x_vals);

    plt.xlabel('MPI Rank x')
    plt.ylabel('Runtime (s.) on MPI Rank x')
    #plt.legend(['Prior Gather', 'Post Gather']);
    #plt.axis('tight')
    plt.ylim([0, 75]);
    plt.xticks(np.arange(0, len(runtimes_total), step=2), np.arange(0, len(runtimes_total), step=2));
    plt.title('Time taken in GDV Vector Calculation');
    plt.savefig( "png_files/load_imb_box_total_time_" + n_procs + ".png", bbox_inches="tight", pad_inches=0.25);



    #fig3 = plt.figure();
    #plt.bar(x_vals, post_gather_avg, color='b');
    #plt.xticks(x_vals);

    #plt.xlabel('MPI Rank x')
    #plt.ylabel('Runtime (s.) on MPI Rank x')
    ##plt.legend(['Prior Gather', 'Post Gather']);
    #plt.axis('tight')
    #plt.title('Time Taken on Each Rank of GDV Vector Calculation');
    #plt.savefig( "png_files/load_imb_bar_post_gather.png", bbox_inches="tight", pad_inches=0.25);


    #fig4 = plt.figure();
    #plt.bar(x_vals, total_time_avg, color='b');
    ##plt.bar(x_vals, post_gather, bottom=prior_gather);
    #plt.xticks(x_vals)

    #plt.xlabel('MPI Rank x')
    #plt.ylabel('Runtime (s.) on MPI Rank x')
    #plt.axis('tight')
    #plt.title('Time Taken on Each Rank of Full Similarity Calculation');
    #plt.savefig( "png_files/load_imb_bar_full.png", bbox_inches="tight", pad_inches=0.25);


if __name__ == "__main__":
    desc = "Generates figure showing runtime load imbalance data"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("filepath",
                        help="Path to run directories storing runtime data")
    parser.add_argument("-r", "--run_max", required=True,
                        help="Number of run directories to use")
    parser.add_argument("-p", "--n_procs", required=True,
                        help="Number of processes parallelized across")
    args = parser.parse_args()
    main( args.filepath, args.run_max )






