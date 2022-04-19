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
            //cout<<i;
            return i;
        }
    }

}

void trace_path_till_root(int root,vector<graph_combination> &list,vector<int> &list1,int value)
{
     list1.push_back(list[value].child);
     list1.push_back(list[value].parent);
     int current_node = list[value].parent;
     int current_value = value;
     current_value--;
     while(current_node != root)
     {
         while(list[current_value].child != current_node)
         {
             current_value--;
         }
         list1.push_back(list[current_value].parent);
         current_node = list[current_value].parent;
         current_value--;
    }
}

int main()
{
    int no_of_nodes;
    vector<graph_combination> list;
    graph_combination parent_child_comb;
    vector<int> list_of_nodes;
    cin >> no_of_nodes;
    int root;
    cout<<"Enter root";
    cin >> root;
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
     for(int i = 0; i < list.size(); i++)
    {
        cout<<list[i].parent;
        cout<<list[i].child;
        cout<<endl;
    }
    if(no_of_nodes >3)
    {
        int start_node,end_node;
        cout<<"Enter start node and end node"<<endl;
        cin>>start_node>>end_node;
        vector<int> list1;
        vector<int> list2;
        int start_value = search_in_the_list(list,start_node);
        int end_value = search_in_the_list(list,end_node);
        trace_path_till_root(root,list,list1,start_value);
        trace_path_till_root(root,list,list2,end_value);

        for(int i=0; i<list1.size(); i++)
        {
            cout<<list1[i]<<" ";
        }
        cout<<endl;
        for(int i=0; i<list2.size(); i++)
        {
            cout<<list2[i]<<" ";
        }
        int list1_start = list1.size()-1;
        int list2_start = list2.size()-1;
        while(list1[list1_start] == list2[list2_start])
        {
            list1_start--;
            list2_start--;
        }
        //int common_ancestor = list1[list1_start] + 1;
        int common_ancestor_list1 = list1_start +1;
        int common_ancestor_list2 = list2_start +1;
        //cout<<endl;
        //cout<<common_ancestor_list1<< " "<<common_ancestor_list2<<endl;
        int list1_min = list1[0];
        for(int i = 0; i <= common_ancestor_list1; i++)
        {
            if(list1[i]<list1_min)
            {
                list1_min = list1[i];
            }
        }
        cout<<list1_min;
        if(list1_min == list1[0])
        {
            int list2_min = list2[0];
            for(int i = 0; i <= common_ancestor_list2; i++)
            {
                if(list2[i]<list2_min)
                {
                    list2_min = list2[i];
                }
            }
            if(list2[0] == list2_min)
            {
                if(list1_min < list2_min)
                {
                    if(list2_min < list1[1])
                    {
                         cout<<"Cycle exist";
                    }
                }
                else
                {
                    if(list2[1] > list1_min)
                    {
                        cout<<"Cycle exist";
                    }
                }
            }
        }

    }
return 0;
}