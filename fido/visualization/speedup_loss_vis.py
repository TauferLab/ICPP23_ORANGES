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
                        if ((int(word) & (int(word) - 1) == 0) and (int(word) != 0)):
                            n_procs.append(int(word));
                            use_line = True;
                    if ( (count == 4) and (count_name == 1) ):
                        graph_size = int(word);
                    if ( (count == 5) and (use_line == True) ):
                        runtimes.append(float(word));
                    count += 1;
            count_name += 1


    # Gather average sets of data
    rows, cols = (6, 10)
    times_per_proc = []*rows
    for power in range(rows):
        times_per_power = [];
        for i in range(len(n_procs)):
            #if (8*(2+power) == n_procs[i]):
            if ( n_procs[i] == 2**power ):
                times_per_power.append(runtimes[i]);
        times_per_proc.append(times_per_power);
    times = np.array(times_per_proc);
    avg_times = np.mean(times_per_proc, axis = 1);
    procs = np.array([1, 2, 4, 8, 16, 32]);
    #high_procs = np.array([16, 24, 32, 40, 48, 56, 64]);


    # Perform analysis to average runs
    # Start printing array forms of stored lists above with mpl
    print(n_procs);
    print(runtimes);
    print(times_per_proc);
    print(avg_times)


    # Prep variables for plotting
    x = np.array(n_procs);
    y = np.array(runtimes);


    # Create figure and plot data and linear regression to it
    plt.figure(figsize=(4, 3))
    ax = plt.axes()
    #print(avg_times[0]/avg_times)
    #print(avg_times[0])
    #plt.plot(procs, avg_times[0]/avg_times,  marker='o');
    #plt.plot(procs, procs, c='b'r);
    bar1 = plt.bar(procs, 100*np.subtract(procs,avg_times[0]/avg_times)/procs, color='b', width=np.array([0.25, 0.5, 1, 2, 4, 8]));


    # Add details to figure
    params = {'mathtext.default': 'regular' }          
    plt.rcParams.update(params)
    plt.xlabel("Number of Processes");
    plt.ylabel("Speedup Loss ($S_i$ - $S_m$)/$S_i$");
    #plt.axis([0, 64, 0, 100]);
    ax.set_xscale('log');
    ax.set_xticks(procs);
    ax.set_xticklabels(['1', '2', '4', '8', '16', '32']);
    plt.title('Fido Speedup Loss as a Percent of the Ideal');
    ax.axis('tight');
    ax.set_ylim(0, 100);
    ax.set_xlim(0, 64);
    #plt.legend(['Actual Speedup', 'Ideal Speedup']);

    #proc_count = 0;
    for box in bar1:
        height = box.get_height();
        #num_deg = num_degs[proc_count];
        plt.text(box.get_x() + box.get_width()/2.0, height, '%d' % int(height), ha='center', va='bottom');
        #proc_count += 1;

    plt.savefig( "png_files/fido_loss_per_nprocs.png",
                 bbox_inches="tight",
                 pad_inches=0.25
               );


if __name__ == "__main__":
    desc = "Generates figure showing speedup loss data"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("filename", 
                        help="Path to file storing runtime data")
    args = parser.parse_args()
    main( args.filename )
