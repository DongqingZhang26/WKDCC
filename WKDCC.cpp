/*
WKDCC algorithm's code
Language: C++
Author: Zhenyu Liu, Dongqing Zhang
Update Date: 2019-8-27
*/
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <string>
#include <bitset>
#include <fstream>
#include <cmath>
#include <list>
#include <cstring>
#include <queue>
#include <set>
using namespace std;
   double min_double=1e-10;
int lg2(int n)
{
    int c=0;
    while(n)
    {
        n>>=1;
        ++c;
    }
    return c;
}
struct Data
{
    int cols,rows,n_cols,c_cols,num_l;
    vector<string> id;
    vector<string> y;
    vector<vector<string> > c_x;
    vector<vector<double> > n_x;
    vector<double> min,max,mid,mean;
    vector<string> mode;
    Data(const char* fn,const char* subfn)
    {
        string path="E:\\data\\";
        path+=fn;
        path+="\\";
        string path_fn=path+"readme.txt";
        ifstream in(path_fn.c_str());
        in>>cols>>cols;
        char type[cols+1];
        int tmp[cols+1];
        int i ,k=0;
        n_cols=c_cols=0;
        for(i=0;i<cols+1;++i)
        {
            in>>type[i]>>tmp[i];
            if(type[i]=='L') num_l=tmp[i];
            if(type[i]=='N') ++n_cols;
            else ++c_cols;
        }
        --c_cols;
        if(c_cols)
        {
            c_x.reserve(c_cols);
            mode.reserve(c_cols);
        }
        if(n_cols)
        {
            min.reserve(n_cols);
            max.reserve(n_cols);
            mean.reserve(n_cols);
            mid.reserve(n_cols);
            n_x.reserve(n_cols);
        }
        for(i=0;i<n_cols;++i)
        {
            n_x.push_back(vector<double>());
            double d;
            in>>d;min.push_back(d);
            in>>d;max.push_back(d);
            in>>d;mean.push_back(d);
            in>>d;mid.push_back(d);
        }
        for(i=0;i<c_cols;++i)
        {
            c_x.push_back(vector<string>());
            string s;
            in>>s;
            mode.push_back(s);

        }
        in.close();
        if(subfn) path_fn=path+subfn;
        else path_fn=path+fn+".txt";
        in.open(path_fn.c_str());
        char str[4096];
        rows=0;
        while(in.getline(str,4096))++rows;
        y.reserve(rows);
        if(subfn) id.reserve(rows);
        for(i=0;i<c_cols;++i)
            c_x[i].reserve(rows);
        for(i=0;i<n_cols;++i)
            n_x[i].reserve(rows);
        in.close();
        in.open(path_fn.c_str());
        int ci=0,ni=0;
        while(k<rows*(cols+1))
        {
            string s;
            double d;
            if(subfn&&k%(cols+1)==0)
            {
                in>>s;
                id.push_back(s);
            }
            if(type[k%(cols+1)]=='L')
            {
                in>>s;
                y.push_back(s);
            }
            else if(type[k%(cols+1)]=='N')
            {
                in>>d;
                n_x[ni].push_back(d);
                ni=(ni+1)%n_cols;
            }
            else
            {
                in>>s;
                c_x[ci].push_back(s);
                ci=(ci+1)%c_cols;
            }
            ++k;
        }
        in.close();
    }

    void to1()
    {
        int i,j;
        for(i=0;i<n_cols;++i)
            for(j=0;j<rows;++j)
               if((max[i]-min[i]))
                      n_x[i][j]=(n_x[i][j]-min[i])/(max[i]-min[i]);
                else n_x[i][j]=1;
    }
    double get_nij(int i,int j)
    {
        return n_x[j][i];
    }
    const string& get_cij(int i,int j)
    {
        return c_x[j][i];
    }
};

struct Less
{
    Data & data;
    int index;
    Less(Data&d,int i):data(d),index(i){}
    double dis1(int j)
    {
        double sum=0;
        int i;
        for(i=0;i<data.n_cols;++i)
            sum+=fabs(data.n_x[i][j]-data.n_x[i][index]);
        for(i=0;i<data.c_cols;++i)
            sum+=data.c_x[i][j]!=data.c_x[i][index];
        return sum;
    }
    double dis2(int j)
    {
        double sum=0;
        int i;
        for(i=0;i<data.n_cols;++i)
            sum+=(data.n_x[i][j]-data.n_x[i][index])*(data.n_x[i][j]-data.n_x[i][index]);
        for(i=0;i<data.c_cols;++i)
            sum+=data.c_x[i][j]!=data.c_x[i][index];
        return sum;
    }
    bool operator()(int a,int b)
    {
        return dis1(a)<dis1(b);
    }
};
struct DSet
{
    vector<int> index;
    double gini;
    string label;
    map<string,vector<int> > distribution;
    vector<double> weight;

    vector<double> center_n;
    vector<string> center_c;
    Data* data;
    DSet(Data* _data=0):data(_data){}
    void clear()
    {
        data=0;
        weight.clear();

        center_c.clear();
        center_n.clear();
        index.clear();
        distribution.clear();
        gini=1;
    }
    void set_index()
    {
        if(!data) return;
        int i;
        index.clear();
        index.reserve(data->rows);
        for(int i=0;i<data->rows;++i)
            index.push_back(i);
    }
    void set_gini()
    {
        gini=1;
        int max=-1;
        for(map<string,vector<int> >::iterator p=distribution.begin();p!=distribution.end();++p)
        {
            gini-=1.0*p->second.size()*p->second.size()/index.size()/index.size();
            if(max==-1||p->second.size()>max)
            {
                label=p->first;
                max=p->second.size();
            }
        }

    }
    void distribute()
    {
        int i;
        distribution.clear();
        for(i=0;i<index.size();++i)
           distribution[data->y[index[i]]].push_back(index[i]);
        set_gini();

  }
    void sample(double per=1,DSet *other=0)
    {
        if(distribution.size()==0) return ;
        if(other)
        {
            other->clear();
            other->data=data;
        }

        index.clear();
        for(map<string,vector<int> >::iterator p=distribution.begin();p!=distribution.end();++p)
        {
            random_shuffle(p->second.begin(),p->second.end());
            int i;
            for(i=0;i<p->second.size()*per;++i)
                index.push_back(p->second[i]);
            if(other)
                for(;i<p->second.size();++i)
                    other->index.push_back(p->second[i]);
        }
        distribute();
        if(other) other->distribute();
    }

    void relief_f(double per,int k)
    {
        if(distribution.size()<2) return;
        weight.clear();
        weight.reserve(data->cols);
        int col;
        for(col=0;col<data->cols;++col)
            weight.push_back(0);

        map<string,vector<int> >sample;
        map<string,vector<int> >::iterator p,q;
        int sample_num=0;
        for( p=distribution.begin();p!=distribution.end();++p)
        {
        int num,num1=lg2(p->second.size()),num2=int(ceil(per*p->second.size()));
          //  int num,num1=11111,num2=int(ceil(per*p->second.size()));
            num=num1<num2?num1:num2;
            if(num>50) num=50;
            sample_num+=num;
            random_shuffle(p->second.begin(),p->second.end());
            sample[p->first].insert(sample[p->first].end(),p->second.begin(),p->second.begin()+num);
        }
        for( p=sample.begin();p!=sample.end();++p)
        {
            int i;
            for(i=0;i<p->second.size();++i)
            {
                for( q= distribution.begin();q!= distribution.end();++q)
                {
                    int K;
                    K=k+1<q->second.size()?k+1:q->second.size();
                    partial_sort(q->second.begin(),q->second.begin()+K,q->second.end(),Less(*data,p->second[i]));

                    int nn;
                    if(p->first==q->first&&K-1)
                        for(nn=0;nn<K;++nn)
                        {
                             for(col=0;col<data->cols;++col)
                                if(col<data->n_cols)
                                    weight[col]-=fabs(data->n_x[col][q->second[nn]]-data->n_x[col][p->second[i]])/(K-1);
                                else
                                    weight[col]-=1.0/(K-1)*(data->c_x[col-data->n_cols][q->second[nn]]!=data->c_x[col-data->n_cols][p->second[i]]);
                        }
                    else if(p->first!=q->first)
                    {

                       int KK=(K-1>0?K-1:1);
                        for(nn=0;nn<KK;++nn)
                        {
                            for(col=0;col<data->cols;++col)
                                if(col<data->n_cols)
                                    weight[col]+=fabs(data->n_x[col][q->second[nn]]-data->n_x[col][p->second[i]])*q->second.size()/(index.size()-distribution[p->first].size())/(KK);
                                else
                                    weight[col]+=double(data->c_x[col-data->n_cols][q->second[nn]]!=data->c_x[col-data->n_cols][p->second[i]])*q->second.size()/(index.size()-distribution[p->first].size())/(KK);

                        }
                    }
                }
            }
        }


        for(col=0;col<data->cols;++col)
        {
            weight[col]/=sample_num;
            if(weight[col]<0) weight[col]=0;
        }

       // 处理全负数权重
        for(col=0;col<data->cols;++col)
            if(weight[col]) break;
        if(col==data->cols)
            for(col=0;col<data->cols;++col)
                weight[col]=1;
    }


    double k_means(DSet** dset,int is_w,double min_max)
    {
        //聚成k簇，k为集合样本类别的数量

        int i,j,k=distribution.size();
        if(k<2) return 0;
        int sizes[k];
        list<int> clus[k];
        double center[2][k][data->n_cols],sum;
        map<string,vector<int> >::iterator p;

        //参数is_w&1时使用特征权重,否则不使用，
        //参数is_w&4时使用不平衡分簇,否则不使用，
        //参数is_w&2时使用kmeans多轮次迭代质心收敛后的划分,一次划分就结束,否则
        vector<double> w;
        if(is_w&1)
        {
            this->relief_f(0.4,1);
            w=weight;
            double max=0;
            for(i=0;i<w.size();++i) if(max<w[i]) max=w[i];
            min_max*=max;
        }
        else
        {
            w.reserve(data->n_cols);
            for(i=0;i<data->n_cols;++i)
                w.push_back(1.0);

        }


        memset(center[0],0,sizeof(double)*k*data->n_cols);
        for(p=distribution.begin(),i=0;i<k;++i,++p)
        {
            for(j=0;j<p->second.size();++j)
                for(int col=0;col<data->n_cols;++col)
                    center[0][i][col]+=data->n_x[col][p->second[j]]/p->second.size();
             sizes[i]=p->second.size();

        }
        //开始聚类迭代

        int turn=0;
        while(turn<6)
        {
            for(i=0;i<k;++i)
                clus[i].clear();
            memset(center[(turn+1)%2],0,sizeof(double)*k*data->n_cols);
            for(i=0;i<index.size();++i)
            {
                int min_k=0;
                if(is_w&4)
                    for(j=1;j<k;++j)
                    {
                        double dis=0,min_dis=0,dd=(double(sizes[j])-sizes[min_k])/(double(sizes[j])+sizes[min_k]);
                        int col;
                        for(col=0;col<data->n_cols;++col)
                        {
                            if(w[col]<min_max) continue;
                            double tmp=center[turn%2][j][col]+(center[turn%2][min_k][col]-center[turn%2][j][col])*dd;
                            dis+=w[col]*(tmp-data->n_x[col][index[i]])*(tmp-data->n_x[col][index[i]]);
                            min_dis+=w[col]*(center[turn%2][min_k][col]-data->n_x[col][index[i]])*(center[turn%2][min_k][col]-data->n_x[col][index[i]]);
                        }
                        if(dis<min_dis-min_double ) min_k=j;

                    }
                else
                {
                    min_k=-1;
                    double min_dis ;
                    for(j=0;j<k;++j)
                    {
                        double dis=0;
                        int col;
                        for(col=0;col<data->n_cols;++col)
                        {
                            if(w[col]<min_max) continue;
                            dis+=w[col]*(center[turn%2][j][col]-data->n_x[col][index[i]])*(center[turn%2][j][col]-data->n_x[col][index[i]]);

                        }

                        if(min_k==-1||dis<min_dis-min_double )
                        {
                            min_k=j;
                            min_dis=dis;
                        }
                    }
                }
                clus[min_k].push_back(index[i]);
                for(j=0;j<data->n_cols;++j)
                    center[(turn+1)%2][min_k][j]+=data->n_x[j][index[i]];
            }
            for(i=0;i<k;++i)
            {
                //sizes[i]=clus[i].size();
                if(clus[i].size())
                    for(j=0;j<data->n_cols;++j)
                        center[(turn+1)%2][i][j]/=clus[i].size();
                else
                    for(j=0;j<data->n_cols;++j)
                        center[(turn+1)%2][i][j]=-10000;
            }

            ++turn;
            if(!(is_w&2))
                break;
            sum=0;
            for(i=0;i<k;++i)
                for(j=0;j<data->n_cols;++j)
                    sum+=fabs(center[0][i][j]-center[1][i][j]);
            if(sum<0.000001) break;
        }

        sum=0;
        for(i=0;i<k;i++)
        {
            dset[i]->clear();
            dset[i]->data=data;
            dset[i]->index.insert(dset[i]->index.end(),clus[i].begin(),clus[i].end());
            dset[i]->distribute();
            dset[i]->center_n.reserve(data->n_cols);
            for(j=0;j<data->n_cols;++j)
                dset[i]->center_n.push_back(center[turn%2][i][j]);
            sum+=dset[i]->gini*dset[i]->index.size()/index.size();

        }

        return gini-sum;
    }
    double k_modes(DSet* *dset,int is_w,double min_max)
    {
        double sum=0;
        int i,j,k=distribution.size();
        list<int> clus[k];
        int sizes[k];
        vector<double> w;
        if(is_w&1)
        {
            this->relief_f(0.3,1);
            w=weight;
            double max=0;
            for(i=0;i<w.size();++i) if(max<w[i]) max=w[i];
            min_max*=max;
        }
        else
        {
            w.reserve(data->c_cols);
            for(i=0;i<data->c_cols;++i)
                w.push_back(1.0);
        }
        string center[2][k][data->c_cols];
        map<string,vector<int> >::iterator p=distribution.begin();
        for(i=0;i<k;++i,++p)
        {
           for(int col=0;col<data->c_cols;++col)
            {
                map<string,int> m;
                map<string,int>::iterator it;
                for(j=0;j<p->second.size();++j)
                    ++m[data->c_x[col][p->second[j]]];
                int max=0;
                string str;
                for(it=m.begin();it!=m.end();++it)
                    if(max<it->second)
                    {
                        max=it->second;
                        str=it->first;
                    }
                center[0][i][col]=str;
            }
            sizes[i]=p->second.size();
        }
        int turn=0;
        while(turn<6)
        {
            for(i=0;i<k;++i)
                clus[i].clear();
            for(i=0;i<index.size();++i)
            {
                int min_k=0;
                if(is_w&4)
                {

                    for(j=1;j<k;++j)
                    {
                        double j_dis=0,min_dis=0,dis=0;
                        int col;
                        for(col=0;col<data->c_cols;++col)
                        {   if(w[col]<min_max)continue;
                            j_dis+=w[col]*(center[turn%2][j][col]!=data->c_x[col][index[i]] );
                            min_dis+=w[col]*(center[turn%2][min_k][col]!=data->c_x[col][index[i]] );
                            dis+=w[col]*(center[turn%2][min_k][col]!=center[turn%2][j][col] );
                        }
                    //    cout<<dis<<endl;
                      //if(j_dis<dis&&min_dis>=dis || !(min_dis<dis&&j_dis>=dis )&& j_dis-min_dis<=-dis+2*dis*sizes[j]/(sizes[j]+sizes[min_k]))
                    if(j_dis==0 || j_dis<dis&&min_dis>=dis ) min_k=j;
                    else if(j_dis>=dis && min_dis>=dis && j_dis < min_dis ) min_k=j;
                   // else if (!(min_dis<dis&&j_dis>=dis )&& j_dis*sizes[min_k]<min_dis*sizes[j] ) min_k=j;



                    }
                }
                else
                {
                    min_k=-1;
                    double min_dis;
                    for(j=0;j<k;++j)
                    {
                        double dis=0;
                        int col;
                        for(col=0;col<data->c_cols;++col)
                        {
                            if(w[col]<min_max)continue;
                            dis+=w[col]*(center[turn%2][j][col]!=data->c_x[col][index[i]] );
                        }

                        if(min_k==-1||dis<min_dis)
                        {
                            min_k=j;
                            min_dis=dis;
                        }
                    }
                }
                clus[min_k].push_back(index[i]);
            }
            for(i=0;i<k;++i)
            {

                sizes[i]=clus[i].size();
                for(j=0;j<data->c_cols;++j)
                {
                    map<string,int> m;
                    for(list<int>::iterator modei=clus[i].begin();modei!=clus[i].end();++modei)
                       ++m[data->c_x[j][*modei]];
                    map<string,int>::iterator p=m.begin();
                    int max=0;
                    string str;
                    for(;p!=m.end();++p)
                        if(p->second>max)
                        {
                            max=p->second;
                            str=p->first;
                        }
                    center[(turn+1)%2][i][j]=str;
                }
            }
            ++turn;
            if(!(is_w&2))
                break;
            sum=0;
            for(i=0;i<k;++i)
                for(j=0;j<data->c_cols;++j)
                    sum+=fabs(center[0][i][j]!=center[1][i][j]);
            if(sum<0.0000001) break;
        }


        sum=0;
        for(i=0;i<k;i++)
        {   dset[i]->clear();
            dset[i]->data=data;
            dset[i]->index.insert(dset[i]->index.end(),clus[i].begin(),clus[i].end());
            dset[i]->distribute();
            dset[i]->center_c.reserve(data->c_cols);
            for(j=0;j<data->c_cols;++j)
                dset[i]->center_c.push_back(center[turn%2][i][j]);
            sum+=dset[i]->gini*dset[i]->index.size()/index.size();
        }
        //cout<<sum<<"\t";
        return gini-sum;
    }
    double k_prototype(DSet* *dset,int is_w,double a,double min_max)
    {

       if(data->c_cols==0) return k_means(dset,is_w,min_max);
       if(data->n_cols==0) return k_modes(dset,is_w,min_max);
        if(a==-1) a=double(data->c_cols)/data->cols;

        int i,j,k=distribution.size();
        list<int> clus[k];
        int sizes[k];
        vector<double> w;

        if(is_w&1) {
            this->relief_f(0.3,1);
            w=weight;
            double max=0;
            for(i=0;i<w.size();++i) if(max<w[i]) max=w[i];
             min_max*=max;
        }
        else
        {
            w.reserve(data->cols);
            for(i=0;i<data->cols;++i)
                w.push_back(1.0);
        }
        double center[2][k][data->n_cols];
        string c_center[2][k][data->c_cols];
        map<string,vector<int> >::iterator p=distribution.begin();
        memset(center[0],0,sizeof(double)*k*data->n_cols);
        for(i=0;i<k;++i,++p)
        {
            int col;
            for(j=0;j<p->second.size();++j)
                for(  col=0;col<data->n_cols;++col)
                    center[0][i][col]+=data->n_x[col][p->second[j]]/p->second.size();

            for(  col=0;col<data->c_cols;++col)
            {
                map<string,int> m;
                map<string,int>::iterator it;
                for(j=0;j<p->second.size();++j)
                    ++m[data->c_x[col][p->second[j]]];
                int max=0;
                string str;
                for(it=m.begin();it!=m.end();++it)
                    if(max<it->second)
                    {
                        max=it->second;
                        str=it->first;
                    }
                c_center[0][i][col]=str;
            }
            sizes[i]=p->second.size();
        }

        int turn=0;
        while(turn<6)
        {
            for(i=0;i<k;++i)
                clus[i].clear();
            memset(center[(turn+1)%2],0,sizeof(double)*k*data->n_cols);
            for(i=0;i<index.size();++i)
            {   int min_k;
                if(is_w&4)
                {
                    min_k=0;


                    for(j=1;j<k;++j)
                    {double j_disc=0,min_disc=0,disc=0;
                    double j_disn=0,min_disn=0,disn=0;
                    double j_dis,min_dis,dis;

                        int col;

                        for(col=0;col<data->n_cols;++col)
                        {
                            if(w[col]<min_max)continue;
                            j_disn+=w[col]*(center[turn%2][j][col]-data->n_x[col][index[i]])*(center[turn%2][j][col]-data->n_x[col][index[i]]);
                            min_disn+=w[col]*(center[turn%2][min_k][col]-data->n_x[col][index[i]])*(center[turn%2][min_k][col]-data->n_x[col][index[i]]);
                            disn+=w[col]*(center[turn%2][min_k][col]-center[turn%2][j][col])*(center[turn%2][min_k][col]-center[turn%2][j][col]);
                        }

                        for(col=0;col<data->c_cols;++col)
                        {
                            if(w[col+data->n_cols]<min_max)continue;
                            j_disc+=(w[col+data->n_cols])*(c_center[turn%2][j][col]!=data->c_x[col][index[i]] );
                            min_disc+=(w[col+data->n_cols])*(c_center[turn%2][min_k][col]!=data->c_x[col][index[i]] );
                            disc+=(w[col+data->n_cols])*(c_center[turn%2][min_k][col]!=c_center[turn%2][j][col] );
                        }

                        j_dis=a*j_disc+(1-a)*j_disn;
                        min_dis=a*min_disc+(1-a)*min_disn;
                        dis=a*disc+(1-a)*disn;

                        //if(j_dis<dis&&min_dis>=dis || !(min_dis<dis&&j_dis>=dis )&&j_dis-min_dis<=-dis+2*dis*sizes[j]/(sizes[j]+sizes[min_k] ))
                       if(j_dis<dis&&min_dis>=dis || !(min_dis<dis&&j_dis>=dis )&& (j_dis-min_dis)*(sizes[j]+sizes[min_k])<double(sizes[j]-sizes[min_k]) * (j_dis+min_dis) )
                            min_k=j;

                    }

                }
                else
                {
                    min_k=-1;
                    double min_dis=0;
                    for(j=0;j<k;++j)
                    {
                        double dis=0;
                        int col;
                        for(col=0;col<data->n_cols;++col)
                        {
                            if(w[col]<min_max) continue;
                            dis+=(1-a)*w[col]*(center[turn%2][j][col]-data->n_x[col][index[i]])*(center[turn%2][j][col]-data->n_x[col][index[i]]);
                        }
                         for(col=0;col<data->c_cols;++col)
                         {
                             if(w[col+data->n_cols]<min_max) continue;
                             dis+=a*(w[col+data->n_cols])*(c_center[turn%2][j][col]!=data->c_x[col][index[i]] );
                         }

                        if(min_k==-1||dis<min_dis)
                        {
                            min_k=j;
                            min_dis=dis;
                        }
                    }

                }

                clus[min_k].push_back(index[i]);
                for(j=0;j<data->n_cols;++j)
                    center[(turn+1)%2][min_k][j]+=data->n_x[j][index[i]];
            }
            for(i=0;i<k;++i)
            {
                sizes[i]=clus[i].size();
                if(clus[i].size())
                    for(j=0;j<data->n_cols;++j)
                        center[(turn+1)%2][i][j]/=clus[i].size();
                else
                    for(j=0;j<data->n_cols;++j)
                        center[(turn+1)%2][i][j]=-10000;
            }
            for(i=0;i<k;++i)
                for(j=0;j<data->c_cols;++j)
                {
                    map<string,int> m;
                    for(list<int>::iterator modei=clus[i].begin();modei!=clus[i].end();++modei)
                       ++m[data->c_x[j][*modei]];
                    map<string,int>::iterator p=m.begin();
                    int max=0;
                    string str;
                    for(;p!=m.end();++p)
                        if(p->second>max)
                        {
                            max=p->second;
                            str=p->first;
                        }
                    c_center[(turn+1)%2][i][j]=str;
                }
            ++turn;
            if(!(is_w&2))
                break;
            double sum=0;
            for(i=0;i<k;++i)
            {
                for(j=0;j<data->n_cols;++j)
                    sum+=fabs(center[0][i][j]-center[1][i][j]);
                for(j=0;j<data->c_cols;++j)
                    sum+=fabs(c_center[0][i][j]!=c_center[1][i][j]);
            }
            if(sum<0.000001) break;
        }
        //cout<<turn<<endl;
        double sum=0;
        for(i=0;i<k;i++)
        {
            dset[i]->clear();
            dset[i]->data=data;
            dset[i]->index.insert(dset[i]->index.end(),clus[i].begin(),clus[i].end());
            dset[i]->distribute();
            dset[i]->center_n.reserve(data->n_cols);
            for(j=0;j<data->n_cols;++j)
                dset[i]->center_n.push_back(center[turn%2][i][j]);
            dset[i]->center_c.reserve(data->c_cols);
            for(j=0;j<data->c_cols;++j)
                dset[i]->center_c.push_back(c_center[turn%2][i][j]);
            sum+=dset[i]->gini*dset[i]->index.size()/index.size();
        }
        return gini-sum;
    }
    void print_distribution(ostream& out)const
    {
        out<<"gini\t"<<gini<<"  "<<label<<" "<<index.size()<<endl;
        map<string,vector<int> > ::const_iterator p=distribution.begin();
        while(p!=distribution.end())
        {
            out<<p->first<<"\t"<<p->second.size()<<endl;
            for(int ii=0;ii<p->second.size();++ii)
              cout<<p->second[ii]<<" ";
            cout<<endl;
            ++p;
        }
    }
    void print_distribution1(ostream& out)const
    {
        out<<index.size()<<"\t"<<distribution.size()<<endl;
        map<string,vector<int> > ::const_iterator p=distribution.begin();
        while(p!=distribution.end())
        {
            out<<p->first<<"\t"<<p->second.size()<<endl;
            ++p;
        }
    }
    //
};
//int xx=0;
struct Node
{   //int x;
    DSet dset;
    int num_child;
    Node** children;
    Node():children(0),num_child(0){}
    Node(Data* data):dset(data),children(0),num_child(0){}
    ~Node()
    {
        int i;
        while(num_child--)
            delete children[num_child];
        delete[]children;
        //cout<<"xx"<<x<<endl;
    }
};
struct Tree
{
    Node root;
    int size,depth;
    int is_w;
    //vector<double> weights;
    Tree(Data* data=0,int w=0):root(data),size(1),is_w(w),depth(1)
    {
        if(data)
        {
            root.dset.set_index();
            root.dset.distribute();
        }
    }
    void grow(Node&cur,int _level,double min_gini,double aa,double min_max )
    {if(_level>depth) depth=_level;
      //if(cur.dset.index.size()<5||cur.dset.gini<min_gini) return ;

    if(cur.dset.index.size()==0)return;
       int k=cur.dset.distribution.size();
       if(k<2) return;
       cur.children=new Node*[k];
       DSet*dsp[k];
       for(int i=0;i<k;++i)
       {
          cur.children[i]=new Node;
          dsp[i]=&cur.children[i]->dset;
       }
       double delta=cur.dset.k_prototype(dsp,is_w,aa,min_max);

       int i,j=0;
       for(i=0;i<k;++i)
           if(cur.children[i]->dset.index.size())
           {
                if(i!=j)
                {
                    Node* t=cur.children[i];
                    cur.children[i]=cur.children[j];
                    cur.children[j]=t;
                }
                ++j;
            }


            cur.num_child=j;

            for( ;j<k;++j)
            {
                delete cur.children[j];
                cur.children[j]=0;
            }
            if(delta>min_double)
            {
               size+=cur.num_child;

               for(i=0;i<cur.num_child;++i)
               {
                   grow(*cur.children[i],_level+1,min_gini,aa,min_max);
               }

            }
            else
            {
                for(i=0;i<cur.num_child;++i)
                   delete cur.children[i];
                delete [] cur.children;
                cur.children=0;
                cur.num_child=0;

            }
    }
    void pr(ostream& out)
    {

        out<<depth<<"\t"<<size<<"\t"<<is_w<<endl;
        queue< Node* > q;
        q.push( &root );
        int i=0,j,k;
        while(!q.empty())
        {
            Node* tmp=q.front();q.pop();
            out<<tmp->num_child<<endl;

            if(tmp->num_child)
            {
                //tmp->dset.print_distribution1(cout);
                if(is_w&1)
                {
                    for(j=0;j<tmp->dset.data->cols;++j)
                        out<<"\t"<<tmp->dset.weight[j];
                    out<<endl;
                }
                map<string,vector<int> >::iterator pp=tmp->dset.distribution.begin();
                for(  j=0;j<tmp->num_child;++j)
                {
                    out<<++i;
                    if(tmp->dset.distribution.size()!=tmp->num_child)
                     out<<"\t"<<tmp->children[j]->dset.index.size();
                    else
                    {out<<"\t"<<pp->second.size();++pp;}
                    for(k=0;k<tmp->dset.data->n_cols;++k)
                        out<<"\t"<<tmp->children[j]->dset.center_n[k];
                    for(k=0;k<tmp->dset.data->c_cols;++k)
                        out<<"\t"<<tmp->children[j]->dset.center_c[k];
                    out<<endl;
                    q.push(tmp->children[j]);
                }
            }
            else
               tmp->dset.print_distribution1(out);
        }

    }
};

void train(const char* dataname,int is_w,int tree_num,double sa,double min_gini,double aa,double min_max)
{
    int x=1,y=1;
    char trainname[32];
    char modelname[128];
    for(x=1;x<=10;++x)
        for(y=1;y<=10;++y)
        {
            sprintf(trainname,"train_%d_%d.txt",x,y);
            sprintf(modelname,"e:\\data\\%s\\model_%d_%d.txt",dataname,x,y);
            Data data(dataname,trainname);
            ofstream out(modelname);
            data.to1();
            out<<tree_num ;
            int i,j,k;
            for(k=0;k<tree_num;++k)
            {
                Tree t(&data,is_w);
                if(sa<1) t.root.dset.sample(sa);
                if(!k)
                {
                    out<<"\t"<<t.root.dset.distribution.size()<<endl;
                    for(map<string,vector<int> >::iterator it=t.root.dset.distribution.begin();it!=t.root.dset.distribution.end();++it)
                       out<<it->first<<endl;
                }
                t.grow(t.root,1,min_gini,aa,min_max);
                t.pr(out);
             }
            out.close();
        }
}


struct Test_tree
{
   struct NLeaf
   {
      vector<double> weight;
      int *index;
      int *size ;
      vector<double>* center;
      vector<string>* center_c;
      int len;
      ~NLeaf()
      {
          delete []center;
          delete []center_c;
          delete []index;
          delete []size;
      }
   };
   struct Leaf
   {
      int size;
      map<string,int> distribution;

   };

   int is_w;

   int depth,size;
   int sample_count;
   vector<pair<bool,void*> > arr;
   void read(ifstream & in,int cols,int cols_c)
   {
       sample_count=0;
       in>>depth>>size>>is_w;
       arr.reserve(size);
       int i=0,j,k;
       for(;i<size;++i)
       {
           int num;
           in>>num;
           if(num)
           {
               NLeaf* p=new  NLeaf;
               p->len=num;
               if(is_w&1){
                   p->weight.reserve(cols+cols_c);
                   for(j=0;j<cols+cols_c;++j)
                   {
                       double dtmp;
                       in>>dtmp;
                       p->weight.push_back(dtmp);
                   }
               }
               p->index=new int[num];
               p->size=new int [num];
               p->center=new vector<double>[num];
               p->center_c=new vector<string>[num];
               for(j=0;j<num;j++)
               {

                   in>>p->index[j]>>p->size[j];
                   if(i==0)  sample_count+=p->size[j];
                   p->center[j].reserve(cols);
                   p->center_c[j].reserve(cols_c);
                   for(k=0;k<cols;++k)
                    {
                        double dtmp;
                        in>>dtmp;
                        p->center[j].push_back(dtmp);
                    }
                    for(k=0;k<cols_c;++k)
                    {
                        string stmp;
                        in>>stmp;
                        p->center_c[j].push_back(stmp);
                    }
                }
               arr.push_back(pair<bool,void*>(1,p));
           }
           else
           {
               Leaf* p=new Leaf;
               int tmp;

               in>>p->size>>tmp;
               if(i==0)  sample_count+=p->size;
               for(j=0;j<tmp;j++)
               {
                   string labeltmp;
                   int inttmp;
                   in>>labeltmp>>inttmp;
                   p->distribution[labeltmp]=inttmp;
                }
               arr.push_back(pair<bool,void*>(0,p));
           }
       }
   }
   ~Test_tree()
   {
       int i;
       for(i=0;i<size;++i)
          if(arr[i].first)
             delete ( NLeaf*)arr[i].second;
          else
              delete (Leaf*)arr[i].second;
    }
};


void test_help(Data& data,int row,Test_tree&tree,map<string,double>&v,int k ,double a,double min_max)
{
    int i=0,j;
    if(a==-1) a=double(data.c_cols)/data.cols;

    while(1)
    {
        if(tree.arr[i].first)
        {
            vector<double> weight;
            Test_tree::NLeaf& centers=*(Test_tree::NLeaf*)tree.arr[i].second;
            if(tree.is_w&1)
            {weight=centers.weight;
                            double max=0;
                            for(int ii=0;ii<weight.size();++ii) if(max<weight[ii]) max=weight[ii];
                            min_max*=max;

            }
            else
                for(j=0;j<data.cols;++j)
                   weight.push_back(1);

            int min=-1;
            double min_dis;
            if(tree.is_w&4)
            {
                if(a>min_double)
                {
                    min=0;
                    for(j=1;j<centers.len;++j)
                        {
                            double j_disc=0,j_disn=0,min_disc=0,min_disn=0,disc=0,disn=0,j_dis,dis;

                            int col;

                            for(col=0;col<data.n_cols;++col)
                            {
                                if(weight[col]<min_max) continue;
                                j_disn+=weight[col]*(data.n_x[col][row]-centers.center[j][col])*(data.n_x[col][row]-centers.center[j][col]);
                                min_disn+=weight[col]*(data.n_x[col][row]-centers.center[min][col])*(data.n_x[col][row]-centers.center[min][col]);
                                disn+=weight[col]*(centers.center[j][col]-centers.center[min][col])*(centers.center[j][col]-centers.center[min][col]);
                            }
                            for(col=0;col<data.c_cols;++col)
                            {    if(weight[col+data.n_cols]<min_max) continue;
                                j_disc+=weight[col+data.n_cols]*(data.c_x[col][row]!=centers.center_c[j][col]);
                                min_disc+=weight[col+data.n_cols]*(data.c_x[col][row]!=centers.center_c[min][col]);
                                disc+=weight[col+data.n_cols]*(centers.center_c[j][col]!=centers.center_c[min][col]);
                            }
                            j_dis=a*j_disc+(1-a)*j_disn;
                            min_dis=a*min_disc+(1-a)*min_disn;
                            dis=a*disc+(1-a)*disn;

                            if(j_dis==0 || j_dis<dis&&min_dis>=dis ) min=j;
                            else if(j_dis>=dis && min_dis>=dis && j_dis < min_dis ) min=j;
                            //else if (!(min_dis<dis&&j_dis>=dis )&& j_dis*centers.size[min]<min_dis*centers.size[j] )                               min=j;

                        //  if(j_dis<dis&&min_dis>=dis || !(min_dis<dis&&j_dis>=dis )&& j_dis-min_dis<=-dis+2*dis*centers.size[j]/(centers.size[j]+centers.size[min]))
                         }
                }
                else
                {
                     min=0;
                    for(j=1;j<centers.len;++j)
                    {
                        double dis=0;
                        min_dis=0;
                        int col;
                        double dd=double(centers.size[j]-centers.size[min])/(centers.size[j]+centers.size[min]);

                        for(col=0;col<data.n_cols;++col)
                        {
                            if(weight[col]<min_max) continue;
                            double tmp=centers.center[j][col]+(centers.center[min][col]-centers.center[j][col])*dd;
                            dis+=weight[col]*(data.n_x[col][row]-tmp)*(data.n_x[col][row]-tmp);
                            min_dis+=weight[col]*(data.n_x[col][row]-centers.center[min][col])*(data.n_x[col][row]-centers.center[min][col]);
                        }

                        if(dis<min_dis-min_double )
                            min=j;

                    }
                }
            }
            else
                for(j=0;j<centers.len;++j)
                {
                    double dis,disc=0,disn=0;
                    int col;
                    for(col=0;col<data.n_cols;++col)
                    {   if(weight[col]<min_max) continue;
                        disn+=weight[col]*(data.n_x[col][row]-centers.center[j][col])*(data.n_x[col][row]-centers.center[j][col]);
                    }

                    for(col=0;col<data.c_cols;++col)
                    {
                        if(weight[col+data.n_cols]<min_max) continue;
                        disc+=weight[col+data.n_cols]*(data.c_x[col][row]!=centers.center_c[j][col]);
                    }

                    dis=a*disc+(1-a)*disn;
                    if(min==-1||dis<min_dis)
                    {
                        min=j;
                        min_dis=dis;
                    }
                }
            i=centers.index[min];
        }
        else
        {
            map<string,double>::iterator p;
            Test_tree::Leaf& distribution=*(Test_tree::Leaf*)tree.arr[i].second;
            if(k<0)
            {
                double max=-1;
                string maxstring;
                map<string,int>::iterator q;
                for(q=distribution.distribution.begin();q!=distribution.distribution.end();++q)
                    if(q->second>max)
                    {
                        max=q->second;
                        maxstring=q->first;
                    }
                    v[maxstring]+=1;
            }
            else
                for(p=v.begin();p!=v.end();++p)
                {
                    if(k)
                        p->second+=(distribution.distribution[p->first]+1.0)/(distribution.size+k);
                    else p->second+=(distribution.distribution[p->first]+0.0)/distribution.size;
                }
            return;
        }
    }
}
void test(const char* dataname,int kk,double aa,double min_max)
{
    char testname[128];
    char modelname[128];
    char resultname[128];
    int x,y,i,j,k,is_w=-1;
    sprintf(resultname,"e:\\data\\%s\\result.txt",dataname);
    double a_sum1=0,a_sum2=0,dp_sum1=0,dp_sum2=0,node_sum1=0,node_sum2=0,sample_sum1=0,sample_sum2=0,p_sum1=0,p_sum2=0,r_sum1=0,r_sum2=0;
    map<string,double> ps_sum1,ps_sum2,rs_sum1,rs_sum2;
    ofstream out(resultname,ios_base::app);
    sprintf(resultname,"e:\\data\\%s\\result1.txt",dataname);
    ofstream out1(resultname,ios_base::app);

    for(x=1;x<=10;++x)
    {
        map<string,map<string,int> >confus;
        double dp_count=0,node_count=0,sample_count=0;
        int sum=0;
        int tree_num;
        for(y=1;y<=10;++y)
        {

            sprintf(modelname,"e:\\data\\%s\\model_%d_%d.txt",dataname,x,y);
            sprintf(testname,"test_%d_%d.txt",x,y);
            Data data(dataname,testname);
            ifstream in(modelname);
            //ofstream out(resultname);
            data.to1();

            in>> tree_num>>k;
            Test_tree trees[tree_num];
            string labels[k];
            for(i=0;i<k;++i)
            {
                in>>labels[i];

            }
            if(y==1)
            {
                for(i=0;i<k;++i)
                    for(j=0;j<k;++j)
                        confus[labels[i]][labels[j]]=0;

            }
            for(i=0;i<tree_num;++i)
            {
                trees[i].read(in,data.n_cols,data.c_cols);
                if(is_w==-1) is_w=trees[i].is_w;
                node_count+=trees[i].size;
                dp_count+=trees[i].depth;
                sample_count+=trees[i].sample_count;
                //cout<<i<<"\t"<<trees[i].sample_count<<endl;
            }
            cout<<modelname<<" "<<y<<endl;
            in.close();

            for(i=0;i<data.rows;++i)
            {
                map<string,double> re;
                for(j=0;j<k;++j) re[labels[j]]=0;
                for(j=0;j<tree_num;++j)
                    test_help(data,i,trees[j],re,kk<=0?kk:k,aa,min_max);

                double max=-1;
                string maxstring;
                for(map<string,double>::iterator p= re.begin();p!=re.end();++p)
                    if(p->second>max)
                    {
                        max=p->second;
                        maxstring=p->first;
                    }
                ++confus[data.y[i]][maxstring];
           }

       }cout<<"ok"<<endl;
       dp_count/=tree_num*10;
       node_count/=tree_num*10;
       sample_count/=tree_num*10;
       dp_sum1+=dp_count;
       dp_sum2+=dp_count*dp_count;
       node_sum1+=node_count;
       node_sum2+=node_count*node_count;
       sample_sum1+=sample_count;
       sample_sum2+=sample_count*sample_count;
       map<string,int> row_sum,col_sum;
       out <<is_w<<"\t"<<kk<<"\t"<<x<<endl;
       map<string,map<string,int> >::iterator p;
       for(p=confus.begin();p!=confus.end();++p)
       {
            out<<p->first<<"\t\t";
            for(map<string, int >::iterator q=p->second.begin();q!=p->second.end();++q)
            {
                out<<q->second<<"\t";
                row_sum[p->first]+=q->second;
                col_sum[q->first]+=q->second;
            }

            sum+=row_sum[p->first];
            out<<endl;
       }
       double ptmp_sum=0,atmp_sum=0,rtmp_sum=0;

       for(p=confus.begin();p!=confus.end();++p)
       {
             double tmp=double(confus[p->first][p->first]);
             if(col_sum[p->first]==0){col_sum[p->first]=1; }
             double ptmp=tmp/col_sum[p->first],rtmp=tmp/row_sum[p->first];
             ps_sum1[p->first]+=ptmp;
             ps_sum2[p->first]+=ptmp*ptmp;
             rs_sum1[p->first]+=rtmp;
             rs_sum2[p->first]+=rtmp*rtmp;
             ptmp_sum+=ptmp;
             rtmp_sum+=rtmp;
             atmp_sum+=tmp;
             out<<p->first<<"\t"<<ptmp<<"\t"<<rtmp<<endl;
       }
       out<<(atmp_sum/=sum)<<"\t"<<(ptmp_sum/=confus.size())<<"\t"<<(rtmp_sum/=confus.size())<<endl;
       a_sum1+=atmp_sum;
       a_sum2+=atmp_sum*atmp_sum;
       p_sum1+=ptmp_sum;
       p_sum2+=ptmp_sum*ptmp_sum;
       r_sum1+=rtmp_sum;
       r_sum2+=rtmp_sum*rtmp_sum;
    }

    out<<(dp_sum1/=10)<<"\t"<<(node_sum1/=10)<<"\t"<<(sample_sum1/=10)<<"\t"<<(a_sum1/=10)<<"\t"<<(p_sum1/=10)<<"\t"<<(r_sum1/=10)<<endl;
    cout<<is_w<<"\t"<<kk<<"\t"<<(dp_sum1)<<"\t"<<(node_sum1)<<"\t"<<(sample_sum1)<<"\t"<<(a_sum1)<<"\t"<<(p_sum1)<<"\t"<<(r_sum1)<<endl;
    out1<<is_w<<"\t"<<kk<<"\t"<<(dp_sum1)<<"\t"<<(node_sum1)<<"\t"<<(sample_sum1)<<"\t"<<(a_sum1)<<"\t"<<(p_sum1)<<"\t"<<(r_sum1)<<endl;
    out<<(sqrt(dp_sum2/10-dp_sum1*dp_sum1))
    <<"\t"<<(sqrt(node_sum2/10-node_sum1*node_sum1))
    <<"\t"<<(sqrt(fabs(sample_sum2/10-sample_sum1*sample_sum1)))
    <<"\t"<<(sqrt(a_sum2/10-a_sum1*a_sum1))
    <<"\t"<<(sqrt(p_sum2/10-p_sum1*p_sum1))
    <<"\t"<<(sqrt(r_sum2/10-r_sum1*r_sum1))
    <<endl;
    map<string,double>::iterator p;
    for(p=ps_sum1.begin();p!=ps_sum1.end();++p)
    {
        out<<p->first<<"\t"<<(ps_sum1[p->first]/=10);
        out<<"\t"<<sqrt(ps_sum2[p->first]/10-ps_sum1[p->first]*ps_sum1[p->first]);
        out<<"\t"<<(rs_sum1[p->first]/=10);
        out<<"\t"<<sqrt(rs_sum2[p->first]/10-rs_sum1[p->first]*rs_sum1[p->first])<<endl;
    }
}


int main()
{
    //train(const char* dataname,int is_w,int tree_num,double sa,double min_gini,double aa,double min_max)
//void test(const char* dataname,int kk,double aa,double min_max)

    char *fn[8]={"chess"};
    double min_gini=0.0,a=-1,min_max=0.2;

    for(int k=0;k<1;k++)
    {


   for(int i=10;i<8;i++)
     {
          train(fn[k],i,1,1,min_gini,a,min_max);
          test(fn[k],0,a,min_max);
     }

     for(int i=3;i<4;i++)
     {
          train(fn[k],i,10,0.7,min_gini,a,min_max );
           for(int j=-1;j<2;j++)
             test(fn[k],j,a,min_max);
     }
}


    return 0;
}



