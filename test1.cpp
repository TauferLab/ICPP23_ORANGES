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

int trace_path_till_root(int root, vector<graph_combination> &list, int iterator)
{
    int current_node = list[iterator].parent;
    int iterator_parent = current_node;
    int minimum_node = current_node;
    while(list[current_node].parent ! = root)
    {
        iterator--;
        if(interator < 0)
        {
            return -999;
        }
        else
        {
            if(list[iterator].child == current_node)
            {
                if(list[iterator].parent < minimum_node)
                {
                    minimum_node = list[iterator].parent;
                }
            
            }
        }

    }
}

int main()
{
    int no_of_nodes;
    vector<graph_combination> list;
    cin >> no_of_nodes;
    int root;
    cin>> root;
    graph_combination parent_child_comb;
    for(int i = 0; i < no_of_nodes; i++)
    {
        cin>>parent_child_comb.parent;
        cin>>parent_child_comb.child;
        list.push_back(parent_child_comb);
    }
    if(no_of_nodes >= 3)
    {
        int iterator = 1;
        while(iterator < no_of_nodes)
        {
            int minimum_node = trace_path_till_root(root,list,iterator);
            cout<<list[iterator].child << minimum_node;
            iterator++;
        }    
    }

    return 0;
}