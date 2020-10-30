

from xml.dom import minidom
import os

# Pass XML file path and filename. It creates graph file and saves it with the provided filename
def graph_creation_from_XML_to_txt(path,filename):
    f= open(filename,"w+")
    xmldoc = minidom.parse(path)
    itemlist = xmldoc.getElementsByTagName('edge') 
    for items in itemlist:
        src=items.attributes['source'].value
        trg=items.attributes['target'].value
        f.write(("{0} {1} {2}\r\n").format(int(src[1:]),int(trg[1:]),1))
    f.close()
    return

def graph_dictionary_generation(graph_parent_dictionary,run_list,slice_number):
    graph_dictionary={}
    for i in range(1,len(graph_parent_dictionary)+1):
        graph_dictionary[i-1]=run_list + "/"+graph_parent_dictionary[i] + "/slices/slice_"+str(slice_number) +".graphml.xml"
    return graph_dictionary

def initialize_dataset(run_list,slice_number):
    # run_list = 'C:/Users/lohit/Desktop/RA/Datasets/Naive-2/naive_reduce/system_quartz/n_procs_11/proc_placement_spread/msg_size_1'
    graph_parent_dictionary={}
    for filename in os.listdir(run_list):
        index = filename[4:].replace('_','')
        graph_parent_dictionary[int(index)]=filename
    # now for slice 4
    graph_dictionary = graph_dictionary_generation(graph_parent_dictionary,run_list,slice_number)

    return graph_dictionary

def main():
    # Converting a XML file to a text file
    xmlpath = r'C:\Users\lohit\Desktop\RA\Datasets\naive_reduce\naive_reduce\naive_reduce\msg_size_1\without_ninja\run_2\slices\slice_3.graphml'
    filename = xmlpath[-16:]
    name = filename.replace('\\','') + ".txt"
    graph_creation_from_XML_to_txt(xmlpath,filename)

    # Converting a group of files based on the folder structure. 
    # This helps us to get all graphs of a single slice as one output
    # For example, if we wanted to compare slice_2 in all runs (run1, 2,3, ...)
     
    run_list = 'C:/Users/lohit/Desktop/RA/Datasets/Naive-2/naive_reduce/system_quartz/n_procs_11/proc_placement_spread/msg_size_1'
    slice_number = 4
    # If we wanted a dictionary of all grpahs belonging to slice 4
    graph_dictionary = initialize_dataset(run_list,slice_number)
    print(graph_dictionary)


if __name__ == '__main__':
    main()

