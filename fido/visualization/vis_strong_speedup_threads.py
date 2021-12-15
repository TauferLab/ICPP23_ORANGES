
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression

import argparse

def main( filepath , run_min , run_max , vis_type ):

    #runtimes_prior_gather_2d = []*int(run_max);
    #runtimes_post_gather_2d = []*int(run_max);
    #runtimes_total_2d = []*int(run_max);

    n_procs = 1
    thread_scales = [1, 2, 4, 8, 16, 32]
    n_runs = (int(run_max)+1)-int(run_min)
    runtime_data = np.zeros((n_runs , len(thread_scales)))

    run_index = 0
    for run_num in range(int(run_min), int(run_max)+1):
        proc_index = 0
        for num_threads in thread_scales:

            if (run_num < 10):
                filename = filepath + "/run_00" + str(run_num) + "/procs_" + str(n_procs) + "/threads_per_proc_" + str(num_threads) + "/runtime_data/runtime_over.txt"
            if (run_num >= 10 and run_num < 100):
                filename = filepath + "/run_0" + str(run_num) + "/procs_" + str(n_procs) + "/threads_per_proc_" + str(num_threads) + "/runtime_data/runtime_over.txt"
            with open(filename, "r") as file:

                #runtimes_prior_gather = [];                                                                                                                                                                                              
                #runtimes_post_gather = [];
                #runtimes_total = [];

                count_line = 0;
                for line in file:
                    if (count_line == 2):
                        count_word = 0;
                        #runtimes_prior_gather = [];
                        #runtimes_post_gather = [];
                        #runtimes_total = [];
                        for word in line.split():
                            #if (count_word == 1):
                                #runtimes_prior_gather.append(float(word));
                                #runtime_data[ run_index , proc_index ] = 
                            #if (count_word == 2):
                            #    runtimes_post_gather.append(float(word));
                            #if (count_word == 3):
                            #    runtimes_total.append(float(word));
                            if (count_word == 5):
                                runtime_data[ run_index , proc_index ] = float(word)
                            count_word += 1;
                    count_line += 1;
                #runtimes_prior_gather_2d.append(runtimes_prior_gather);
                #runtimes_post_gather_2d.append(runtimes_post_gather);
                #runtimes_total_2d.append(runtimes_total);

            proc_index += 1
        run_index += 1

    #print(runtimes_total);
    #print(runtimes_prior_gather);
    #print(runtimes_post_gather);

    #x_vals = np.array(range(0, len(runtimes_total)));
    #print(np.shape(x_vals))
    #prior_gather = np.array(runtimes_prior_gather_2d);
    #post_gather = np.array(runtimes_post_gather_2d);
    #total_time = np.array(runtimes_total_2d);

    print(np.shape(runtime_data))
    print(runtime_data)
    #print(np.shape(prior_gather));
    ax = 0;
    #prior_gather_avg = np.mean(prior_gather, axis=ax);
    #post_gather_avg = np.mean(post_gather, axis=ax);
    #total_time_avg = np.mean(total_time, axis=ax);
    #x_vals = np.array(range(0, len(runtimes_total)));
    x_vals = thread_scales
    #print(total_time_avg);
    #gather_addition = np.subtract(post_gather_avg, prior_gather_avg);
    #print(gather_addition)
    #print(x_vals)

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

    fig = plt.figure()
    ax = plt.axes()
    #plt.bar(x_vals, prior_gather_avg, color='b');
    #plt.boxplot(prior_gather, positions=x_vals);

    if ( vis_type == "runtime" ):
        parts = ax.violinplot(runtime_data, positions=x_vals, widths=2, showmedians=True, showextrema=True)
        for vp in parts['bodies']:
            vp.set_facecolor('olive')
            vp.set_edgecolor('black')
            vp.set_alpha(1)
        parts['cbars'].set_linewidths(1)
        parts['cbars'].set_edgecolors('black')
        parts['cmins'].set_linewidths(2)
        parts['cmins'].set_edgecolors('black')
        parts['cmaxes'].set_linewidths(1)
        parts['cmaxes'].set_edgecolors('black')
        parts['cmedians'].set_linewidths(1)
        parts['cmedians'].set_edgecolors('black')
        plt.xticks(x_vals)
        plt.xlabel('Number of MPI Processes Used')
        plt.ylabel('Runtime (s.)')
        plt.ylim(bottom=0)
        plt.savefig( "png_files/fido_strong_scaling_runtime.png", bbox_inches="tight", pad_inches=0.25)
    
    elif ( vis_type == "speedup" ):
        runtime_avg = np.mean(runtime_data , axis=0)
        plt.plot(run_scales, run_scales, c='k', marker='x')
        plt.plot(run_scales, runtime_avg[0]/runtime_avg, c='b', marker='o')
        plt.xticks(x_vals)
        plt.xlabel('Number of MPI Processes Used')
        plt.ylabel('Runtime Speedup over Lowest Process Count')
        plt.ylim(bottom=0)
        plt.legend(['Ideal Speedup', 'Actual Speedup']);
        plt.savefig( "png_files/fido_strong_scaling_speedup.png", bbox_inches="tight", pad_inches=0.25)
        
    elif ( vis_type == "loss" ):
        runtime_avg = np.mean(runtime_data , axis=0)
        loss_vals = 100*np.subtract(run_scales,runtime_avg[0]/runtime_avg)/run_scales
        bar1 = plt.bar(run_scales, loss_vals, color='b', width=np.array([0.5, 1, 2, 4, 8]));
        plt.xticks(x_vals)
        plt.xlabel('Number of MPI Processes Used')
        plt.ylabel('Speedup Loss as a Percentage of the Ideal')
        ax.set_xscale('log')
        ax.set_xticklabels(x_vals)
        plt.ylim(bottom=0)
        plt.savefig( "png_files/fido_strong_scaling_loss.png", bbox_inches="tight", pad_inches=0.25)

    #plt.legend(['Prior Gather', 'Post Gather']);
    #plt.axis('tight')
    #plt.ylim([0, 75]);
    #plt.xticks(np.arange(0, len(runtimes_total), step=2), np.arange(0, len(runtimes_total), step=2));
    #plt.title('Time Taken on Each Rank of GDV Vector Calc');
    #plt.savefig( "png_files/strong_scaling.png", bbox_inches="tight", pad_inches=0.25);


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
    parser.add_argument("--vis_type",
                        help="Type of visualization to generate (options: runtime, speedup, loss)")
    parser.add_argument("-rn", "--run_min", required=True,
                        help="Starting run directory index")
    parser.add_argument("-rx", "--run_max", required=True,
                        help="Final run directory index")
    args = parser.parse_args()
    main( args.filepath, args.run_min, args.run_max , args.vis_type )






