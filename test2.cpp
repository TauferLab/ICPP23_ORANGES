#include<cstdlib>
#include<iostream>
#include <algorithm>
#include<vector>
using namespace std;

class graph_combination
{
   public:
   int parent;
   int child;
};

class subgraph
{
    public:
    int root;
    int no_of_nodes; 
    vector<graph_combination> list;
    vector<int> list_of_nodes;
    vector<int> end_nodes;
    bool is_subgraph;
};

int search_in_the_list(vector<graph_combination> &list,int j){
    for(int i = 0; i < list.size(); i++)
    {
        if(list[i].child == j)
        {
            return i;
        }
    }

}
int trace_path_till_root(int root, vector<graph_combination> &list, int value)
{
    int actual_value = list[value].child;
    int parent_value = list[value].parent;
    int current_node = parent_value;
    int minimum = list[value].child;
    int minimum_in_loop = list[value].parent;
    int iterator = value - 1;
    if(minimum < parent_value){
        while(current_node != root){
            while(list[iterator].child == current_node)
            {
            iterator--;
            }
            if(list[iterator].parent < minimum)
            {
                break;
            }
            if(minimum_in_loop > list[iterator].parent &&(list[iterator].parent < parent_value))
            {
                minimum_in_loop = list[iterator].parent;
                cout<< list[value].child << list[iterator].parent;
            }
            current_node = list[iterator].parent;
            iterator--;
        }
    }
            ADJ_Bundle node_temp = graph[list[value].child];
            for(int p = 0; p<node_temp.ListW.size(); p++)
            {
               if(node_temp.ListW[p].first == list[iterator].parent)
               {
                   cout<<"Node Exists"
               }
              
            }
    
}

int main()
{
    int no_of_nodes;
    vector<graph_combination> list;
    vector<int> list_of_nodes;
    cin >> no_of_nodes;
    int root;
    cout<<"Enter root";
    cin >> root;
    graph_combination parent_child_comb;
    int number;
    cout<<"Enter list of nodes"<<endl;
    for(int i = 0; i < no_of_nodes; i++)
    {
        cin>>number;
        list_of_nodes.push_back(number);
    }

    cout<<"enter parent child comb"<<endl;
    for(int i = 0; i < no_of_nodes-1; i++)
    {
        cin>>parent_child_comb.parent;
        cin>>parent_child_comb.child;
        list.push_back(parent_child_comb);
    }

    if(no_of_nodes >= 3)
    {
       for(int j = 2; j < no_of_nodes; j++)
       {
           int value = search_in_the_list(list,list_of_nodes[j]);
           int parent_value = list[value].parent;
           trace_path_till_root(root,list,value);

       }
    }

    return 0;
}