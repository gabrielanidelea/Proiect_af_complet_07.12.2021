#include<iostream>
#include<fstream>
#include<queue>
#include<vector>
#include<stack>
#include <bits/stdc++.h>
using namespace std;
ifstream f("darb.in");
ofstream g("darb.out");
const int MAX=100001;
const int INF=0x3f3f3f3f;
class graf{
private:
    int n; //nr noduri //pt pb disjoint reprezinta nr de multimi
    int m; // nr muchii // pt pb disjoint repr nr de operatii
    vector <vector<pair<int, int>>> vector_costuri; //(pt apm)
    int vvect_disjoint[MAX]; // pt disjoint
    vector <vector<pair<int, int>>> vect_dijk;
    vector<vector<pair<int, int>>>vect_costuri_bellman;
    vector<int> graf_mf[1001]; // mf
    bool vizi_mf[1001]; // mf
    int capacitate[1001][1001], tati[1001]; // mf
    int mat[1001][1001]; // roy-floyd
    vector<int>v_darb[100001]; //darb
    int di[100001], nod_d; //vect de dst pt darb + nod

public:
    //---------Tema 1---------------
    //DFS - declarari
    int N, M, S, vizitat[100010], minime[100010];
    vector<int> stocare[100010];
    queue<pair<int, int>> coada;

    //DFS - declarari
    //int N, M;
    vector<vector<int>> D;
    vector<int>ok; //pt a marca nodul ca vizitat sau nu - 0 sau 1
    graf(){}

    //CTC - declarari
    //int N, M;
    int k=0, viz[100005], viz_t[100005];
    vector <int>  rez[100005], G[100005], G_t[100005];
    stack <int> stiv;

    //CRITICAL CONNECTIONS - declarari
    vector<vector<int>>D_C;
    vector<vector<int>>ans;
    vector<int>d;
    vector<int>l;
    vector<int>tata;
    int kk=0, ck=0;

    //BICONEX - declarari
    //int N, M;
    int q=0;
    vector<int>D_B[100005];
    vector<vector<int>>rezultat;
    stack<pair<int, int>> st;
    int d_B[100005]={0}, l_B[100005]={0};


    //HAVEL HAKIMI - declarari
    //int N; // nr de elem din vector
    vector<int> v; //vectorul de grade


    //SORTARE TOPOLOGICA - declarari
    vector<int>D_ST[100005];
    stack<int>s;
    bool vizitat_st[MAX]={0};

//-----------------------------------------------------------------------
    //BFS
   void bfs(int start)
   {
        pair<int, int> curent;
        coada.push({start, 0}); // dist de la start la el insusi e 0
        minime[start] = 0; // momentan, minime va avea val 0
        vizitat[start] = 1; // vizitez nodul curent
        while(!coada.empty())
        {
            curent = coada.front();
            minime[curent.first] = curent.second;
            coada.pop();
            for(auto i : stocare[curent.first]) // parcurgem primele elem din pereche
                if(!vizitat[i]){
                    coada.push({i, curent.second+1}); //cresc contorul
                    vizitat[i] = 1; // vizitez nodul i
                }
        }

    }
    void creare_graf(int N, int M, int S)
    {
        this-> N = N;
        this-> M = M;
        this-> S = S;
        for(int i = 0; i< 100010; i++)
                minime[i] = -1; //umplem tabl de dist minime cu -1
    }
    void creare_adiacente(int M)
    {
        int nod1, nod2;
        for(int i = 1; i<= M; i++)
        {
            f>>nod1>>nod2;
            stocare[nod1].push_back(nod2);
        }

    }
//-------------------------------------------------------------------------------------------
    //DFS
    void dfs(int start)
    {
        ok[start] = 1;
        for(auto i: D[start]) // parcurg D
            if(ok[i]==0)
            {
                ok[i] =1;
                dfs(i);
            }
    }
    void creare_graf_Dfs(int N, int M)
    {
        this-> N = N;
        this -> M =M;
        ok=vector<int>(N+1);
        D=vector<vector<int>>(N+1);
    }
    void creare_adiacente_Dfs(int M)
    {
         int nod1, nod2;
         for(int i = 1; i<= M; i++)
        {
            f>>nod1>>nod2;
            D[nod1].push_back(nod2);
            D[nod2].push_back(nod1);
        }
    }
    int comp_conexe(int N)
    {
        int k= 0;
        for(int i =1; i<= N; i++)
            if(ok[i]==0)
            {
                k++;
                dfs(i);
            }
        return k;
    }
//----------------------------------------------------------------------------------------
    //CTC
    void dfs_ctc1(int nod1)
    {
        viz[nod1]=1; //vizitam nodul1 de start
        for(auto nod2: G[nod1]) //iteram prin vectorul de noduri1 cu nod2
            if(!viz[nod2]) //daca nod2 nu a fost vizitat
                dfs_ctc1(nod2);//apelam recusiv fctia de parcurgere in adancime pt nod2
        stiv.push(nod1); //adaugam nodul1 in stiva
    }
    void dfs_ctc2(int nod1)
    {
        viz_t[nod1]=1; //marcam in tabl viz_t nodul1
        rez[k].push_back(nod1); //retinem in dreptul nr de comp tari conexe indicativul nodului marcat
        for(auto nod2: G_t[nod1])
            if(!viz_t[nod2])
                dfs_ctc2(nod2);
    }
    void creare_graf_ctc(int N, int M)
    {
        this-> N=N;
        this-> M=M;
    }
    void creare_adiacente_ctc(int M)
    {
        for(int i=1; i<=M; i++)
        {
            int nod1, nod2;
            f>>nod1>>nod2;
            G[nod1].push_back(nod2);
            G_t[nod2].push_back(nod1);
        }
    }
    void fct_ctc(int N)
    {
        for(int i=1; i<=N; i++)
            if(viz[i]==0)
                dfs_ctc1(i); //pt fiecare nod nevizitat aplicam dfs_ctc1
        while(!stiv.empty()) //cat timp mai avem elemente in stiva
        {
            if(viz_t[stiv.top()]==0) //daca in viz_t nu a fost vizitat elem din vf-ul stivei
            {
                dfs_ctc2(stiv.top()); //facem dfs_ctc2
                k=k+1; //increm nr de comp
            }
            stiv.pop();
        }
        g<<k<<"\n";
        for(int i=0; i<k; i++)
        {
            for(auto it: rez[i])
                g<<it<<" ";
            g<<"\n";
        }
    }
//------------------------------------------------------------------------------------

    //CRITICAL CONNECTIONS - leetcode
    void dfs_critical_connections(int nod)
    {
        d[nod]=kk++; //vizitam fiecare nod pt ambii vect
        l[nod]=kk++;
        for(auto it: D_C[nod]) //iteram prin nodurile grafului
        {
            if(d[it]==-1)//daca nodul curent nu este vizitat
            {
                tata[it]=nod; //ii atribuim nodului curent it ca parinte pe nod
                dfs_critical_connections(it); //facem dfs din nodul curent
                l[nod]=min(l[nod], l[it]); //mergem pe varianta minima
                if(d[nod]<l[it]) //verificam daca avem nod de intoarcere
                    ans.push_back({nod, it});
            }
            else if(it!=tata[nod]) //daca nodul adiacent nu este parinte
                l[nod]=min(l[nod], d[it]); //minimaliz valoarea l
        }
    }
    void criticalConnections(int N, int M) {
        tata.resize(N,-1);
        l.resize(N,-1);
        D_C.resize(N);
        d.resize(N,-1);
        for(int i=0; i<M; i++)
        {
            int nod1, nod2;
            f>>nod1>>nod2;
            D_C[nod1].push_back(nod2);
            D_C[nod2].push_back(nod1);
        }
        for(int i=0; i<N; i++)
            if(d[i]==-1)
                    dfs_critical_connections(i); //facem dfs pe nodurile nevizitate
        for(int i=0; i<ans.size(); i++)
        {
            g<<"[";
            for(auto it: ans[i])
                g<<it<<" ";
            g<<"]";
        }
    }
//-----------------------------------------------------------------------

    // BICONEX
    void bic(int nod1, int nod2)
    {
        vector<int>local;
        while(st.top().first!=nod1 && st.top().second!=nod2)
        {
            int y=st.top().second;
            local.push_back(y);
            st.pop();
        }
        local.push_back(nod2);
        local.push_back(nod1);
        st.pop();
        rezultat.push_back(local);
    }
    void dfs_bic(int nod_curent)
    {
        q++;
        d_B[nod_curent]=q; //marcam ambele tablouri
        l_B[nod_curent]=q;
        for(auto nod_adiac: D_B[nod_curent])//iteram prin graf
           {
               if(d_B[nod_adiac]) //daca nodul este vizitat
                    l_B[nod_curent]=min(l_B[nod_curent], d_B[nod_adiac]); //pastram min
                else
                {//daca nu este vizitat
                    st.push({nod_curent, nod_adiac}); //retinem in sitva nodul curent si perechea sa
                    dfs_bic(nod_adiac);//facem dfs din nod_adiac
                    l_B[nod_curent]=min(l_B[nod_adiac], l_B[nod_curent]); //pastram min
                    if(l_B[nod_adiac]>=d_B[nod_curent])
                        bic(nod_curent, nod_adiac);//aplicam fct bic pt nodul curent si perechea sa
                }
           }
    }
    void creare_graf_bic(int N, int M)
    {
        this->N=N;
        this->M=M;
    }
    void creare_adiacente_bic(int M)
    {
        int nod1, nod2;
         for(int i = 1; i<= M; i++)
        {
            f>>nod1>>nod2;
            D_B[nod1].push_back(nod2);
            D_B[nod2].push_back(nod1);
        }
    }
    void afisare_bic()
    {
        int t=0;
        for(int i=0; i<rezultat.size(); i++)
            t+=1;
        g<<t<<"\n";
        for(int i=0; i<rezultat.size(); i++) //pt fiecare componenta biconexa
        {
            for(int j=0; j<rezultat[i].size(); j++)
                g<<rezultat[i][j]<<' ';// afisam nodurile care ii apartin
            g<<endl;
        }
    }
//--------------------------------------------------------------------

    //HAVEL HAKIMI
    bool HavelOK(int N, vector<int>&v)
    {
        while(1)
        {
            sort(v.begin(), v.end(), greater<>()); // sortam in ord descresc
            if(!v[0])
                return 1; //prima conditie de stop=>daca toate elem sunt 0, doar atunci se verifica algo Havel Hakimi
            int x=v[0]; //retinem primul elem
            v.erase(v.begin());//in stergem pe x din vector
            if(x>v.size())
                return 0; //a doua conditie de stop=>daca nu mai am destule elem atunci
            for(int i=0; i<x;i++)
            {
                v[i]=v[i]-1;//decrementam urmatoarele x elem cu 1
                if(v[i]<0)
                    return 0;//a treia conditie de stop=>daca intampinam un element negativ
            }
        }
    }
//---------------------------------------------------------
    //SORTARE TOPOLOGICA
    void dfs_st(int nod)
    {
        vizitat_st[nod]=1;
        for(auto nod_adj: D_ST[nod])
            if(vizitat_st[nod_adj]==0) //daca nodul adiacent sau nu a fost vizitat, facem dfs din el
                        dfs_st(nod_adj);
        s.push(nod); // adaugam nodul curent in stiva
    }
    void sortare_topologica()
    {
        for(int i=1; i<=N; i++)
            if(vizitat_st[i]==0) // daca nodul nu a fost vizitat
                dfs_st(i);

    }
    void afisare_sortare_topo()
    {
        while(!s.empty())
        {
            int var = s.top();
            g<<var<<" ";
            s.pop();
        }
    }
    void creare_graf_st(int N, int M)
    {
        this->N=N;
        this->M=M;

    }
    void creare_adiacente_st(int M)
    {
        int nod1, nod2;
        for(int i=0; i<M; i++)
        {
            f>>nod1>>nod2;
            D_ST[nod1].push_back(nod2);
        }

    }

    //------------------------Tema 2 -----------------------------
    //APM
    void citire_APM(); // functie pt creare gf
    void algoritm_APM(); // functie pt implementarea algoritmului lui prim



    //Disjoint
    int gaseste_radacina(int nod);
    void uneste(int nod1, int nod2);
    void algo_disjuncte();

    //Dijkstra
    void algo_dijk();

    //Bellman_Ford
    void algo_bellman_ford();


//--------------------Tema 3---------------------------------
    //MaxFlow
    void citire_mf();
    bool bfs_mf();
    void maxfl();
    void fct_mf();


    //Roy-Floyd
    void citire_rf();
    void fct_royf();

    //Diametrul unui Arbore
    void citire_darb();
    void dfs_darb(int nod1, int nod2, int &dimax);
    void fct_darb();

};
//----------darb---------
void graf:: citire_darb()
{
    int nod1, nod2;
    f >> n;
    for(int i = 1; i < n; i++)
    {
        f >> nod1 >> nod2;
        v_darb[nod1].push_back(nod2);
        v_darb[nod2].push_back(nod1);
    }
}
void graf:: dfs_darb(int nod1, int nod2, int &dimax)
{
    for(auto i: v_darb[nod1])
    {
        if(i == nod2) // daca nodul i este egal cu nodul 2
            continue; //iesim din if
        else
        { //daca este diferit
            di[i] = di[nod1] + 1; // adaugam in dreptul distantei nodului i distanta nodului1 + 1
            if(di[i] > dimax) // daca distanta din e mai mare decat diametrul maxim
            {
                dimax = di[i];
                nod_d = i; // nodul declarat global va primi val din i
            }
            dfs_darb(i, nod1, dimax); // apelam recurs dfs
        }
    }

}
void graf::fct_darb()
{
    int dimax = 0; // diam max e initial 0
    dfs_darb(1, 1, dimax);
    di[nod_d] = 0;
    dimax = 0;
    dfs_darb(nod_d, nod_d, dimax);
    g << di[nod_d] + 1; // afisam dist afer nodului updatat global + 1
}
//---------darb----------
//---------roy-floyd----------
void graf:: citire_rf()
{
    f >> n;
    for(int i = 1; i<=n; i++)
        for(int j = 1; j<=n; j++)
            f >> mat[i][j];

}
void graf:: fct_royf()
{
    for(int p = 1; p <= n; p++)
    {
        for(int i = 1; i<= n; i++)
        {
            for(int j = 1; j <= n; j++)
            {
                if(mat[i][p] && mat[p][j] &&(!mat[i][j] || mat[i][j] > mat[i][p] + mat[p][j]) && i!=j)
                {
                    mat[i][j] = mat[i][p] + mat[p][j];
                }
            }
        }
    }
    for(int i = 1; i<= n; i++)
    {
        for(int j = 1; j<=n; j++)
        {
            g << mat[i][j] << " ";

        }
        g << "\n";
    }
}
//---------roy-floyd--------------------
//---------maxflow----------------------

void graf:: citire_mf()
{
    f >> n >> m;
    for(int i = 0; i < m; i++)
    {
        int nod1, nod2;
        f >> nod1 >> nod2;
        f >> capacitate[nod1][nod2];
        graf_mf[nod1].push_back(nod2);
        graf_mf[nod2].push_back(nod1);
    }
}
bool graf:: bfs_mf()
{
    queue<int> coada_mf;
    coada_mf.push(1); //adaug 1 in coada
    vizi_mf[1] = true; //il marchez
    while(!coada_mf.empty()) // cat timp coada este nevida
    {
        auto nod_curent = coada_mf.front(); // tinem in nod_curent prima var din coada
        coada_mf.pop();
        for(auto i: graf_mf[nod_curent])
        {
            if(vizi_mf[i]!=0) continue; // daca vizi_mf[i] este vizitat
            if(capacitate[nod_curent][i] == 0) continue; // daca in matrix nju sunt unite
            vizi_mf[i] = true; // marchez nodul vecin
            tati[i] = nod_curent; // pun in vect de tati in dreptul n odului vecin, nodul curent
            coada_mf.push(i); // adaug vecinul in coada

        }
    }
    bool ok = vizi_mf[n]; // ok va fi dat de ultimul nod viz
    return ok;
}
void graf:: maxfl()
{
    int mflx, nod_curent;
    for(auto j: graf_mf[n]) //parc graful
    {
        if(!vizi_mf[j]) continue; // daca nodul j nu a fost viz
        mflx = capacitate[j][n]; //tinem capacitatea aferenta lui j si n in mflx
        nod_curent = j; // tinem in var nod_curent variabila j
        while(tati[nod_curent]!=0) // cat timp tatal nodului curent are atribuita o valoare
        {
            mflx = min(mflx, capacitate[tati[nod_curent]][nod_curent]); //tinem in var mflx min dintre capacit si mflx actual
            nod_curent = tati[nod_curent]; // punem in nodul curent tatal nodului curent
        }
        capacitate[j][n] = capacitate[j][n] - mflx;
        capacitate[n][j] = capacitate[n][j] + mflx;
        nod_curent = j; // ii redam nodului curent valoarea j
        while(tati[nod_curent]!=0) //reluam while-ul
        {
            capacitate[tati[nod_curent]][nod_curent] = capacitate[tati[nod_curent]][nod_curent] - mflx;
            capacitate[nod_curent][tati[nod_curent]] = capacitate[nod_curent][tati[nod_curent]] + mflx;
            nod_curent = tati[nod_curent];

        }
    }
}
void graf:: fct_mf()
{
    int rezultat = 0;
    while(bfs_mf()!=0)//cat timp se poate apela bfs_mf()
    {
        maxfl();
        memset(vizi_mf, 0, (n+1)*sizeof(bool)); //setam vizi_mf si tati de fiecare data
        memset(tati, 0, (n+1)*sizeof(int));
    }
    for(auto i: graf_mf[1]) rezultat = rezultat + capacitate[i][1];
    g << rezultat << "\n";
}

//--------maxflow---------------
//--------bellmanford------------------
void graf:: algo_bellman_ford()
{
    f >> n >> m;
    int nod1, nod2, valoare;
    vect_costuri_bellman.resize(n + 1);
    for(int i = 0; i < m; i++)
    {
        f >> nod1 >> nod2 >> valoare;
        vect_costuri_bellman[nod1].push_back(make_pair(nod2, valoare));
    }
    vector<int> distante(n + 1, INF);
    vector<int> lazy(n + 1, 0);
    queue<int> coada_bell;
    vector<bool> ok_coada_bell(n + 1, false);
    int nod_curent, nod_vecin, cost_curent;

    ok_coada_bell[1] = true; // marcam nodul de plecare ca fiind continut de coada
    distante[1] = 0; // dist de la nodul 1 la el insusi e 0
    coada_bell.push(1);

    while(!coada_bell.empty()) // cat timp coada este nevida
    {
        nod_curent = coada_bell.front(); // luam primul nod din coada
        coada_bell.pop(); // ii dam pop
        ok_coada_bell[nod_curent] = false; // il de-marcam

        for(auto it: vect_costuri_bellman[nod_curent]) // parcurgem vectorul de costuri pt nodul curent
        {
            cost_curent = it.second;
            nod_vecin = it.first;
            if(distante[nod_curent] + cost_curent < distante[nod_vecin]) // daca suma dintre distanta aferenta nodului curent + costul curent este mai mica decat dist aferenta nodului vecin
            {
                distante[nod_vecin] = distante[nod_curent] + cost_curent; // adaugam la vect de distante pt nodul vecin aceasta suma
                lazy[nod_vecin]++; // incrementam vectorul folosit drept contor pt nodul vecin
                if(lazy[nod_vecin] == n) // daca ciclul are n
                {
                    g << "Ciclu negativ!";
                    return;
                }
                if (!ok_coada_bell[nod_vecin]) // daca ok_coada_bell[nod_vecin] nu a fost vizitat
                {
                    coada_bell.push(nod_vecin); // adaugam nodul vecin in coada
                    ok_coada_bell[nod_vecin] = true; // il marcam
                }


            }
        }

    }
    for(int i = 2; i<= n; i++)
    {
          if(distante[i] == INF) g << "0 ";
          else
            g << distante[i] << " ";
    }


}
// --------------bellman-ford----------------

//---------------dijkstra-----------
void graf:: algo_dijk()
{
    int nod1, nod2, cost;
    f >> n >> m;
    vect_dijk.resize(n + 1);
    for(int i = 1; i <= m; i++)
    {
        f >> nod1 >> nod2 >> cost;
        vect_dijk[nod1].push_back(make_pair(nod2, cost));
    }
    vector<int>costx(n + 1, 100001); //vect de costuri/distante
	priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> str_djk;
	vector<bool>viz_djk(n + 1, 0); //vector vizite
	int nod_curent;
	costx[1] = 0; //dist pt nodul 1 va fi 0 initial
	str_djk.push(make_pair(0, 1));
	while (!str_djk.empty()) //cat timp pq este nevida
	{
		nod_curent = str_djk.top().second; // luam nodul din top
		str_djk.pop(); // dam pop
		if(viz_djk[nod_curent] == 0) // daca nodul extras nu a fost vizitat
		{
			viz_djk[nod_curent] = 1; // il vizitam
			for(auto it: vect_dijk[nod_curent]) // parcurgem vect_dijk pt nodul curent
			{
				if(costx[nod_curent] + it.second < costx[it.first]) //verif daca suma dintre costul aferent nodului curent si costul curent este mai mica decat costul aferent nodului vecin(it)
				{ //daca da
					costx[it.first] = costx[nod_curent] + it.second; //adunam la costul aferent nodului vecin aceasta suma
					str_djk.push(make_pair(costx[it.first], it.first)); // adaugam in pq perechea formata din costul aferent nodului vecin s, respectiv, nodul vecin
				}
			}
		}

	}
	for (int i = 2; i <= n; i++)
	{
		if(costx[i] != 100001)g << costx[i] << " "; // afisam costurile
		else
			g<<"0 ";

	}
}

//--------------dijkstra--------------
//-------------disjoint-----------------
int graf:: gaseste_radacina(int nod) // cauta radacina nodului curent
{
    if (vvect_disjoint[nod] != nod) // cat timp parintele nodului e dif de nod aplicam recursiv fct
        vvect_disjoint[nod] = gaseste_radacina(vvect_disjoint[nod]);
    return vvect_disjoint[nod]; // returnam radacina noduluji curent
}
void graf:: uneste(int nod1, int nod2)
{
    int radacina1, radacina2; // luam cele doua noduri
    radacina1 = gaseste_radacina(nod1); // le gasim punctul de plecare
    radacina2 = gaseste_radacina(nod2);
    vvect_disjoint[radacina2] = radacina1; // pastram in vvect in dreptul radacinii2 radacina1
}
void graf:: algo_disjuncte()
{
    f >> n >> m;
    for(int i = 0; i < n; i++)
        vvect_disjoint[i] = i; // pt fiecare multime punem in dreptul ei in vvect indicele i
    for(int i = 0; i < m; i++) // parc nr de operatii
    {
        int op, nod1, nod2;
        f >> op >> nod1 >> nod2;
        if(op == 1) uneste(nod1, nod2); // daca suntem pe op 1 aplicam functia uneste
        else if (op == 2)
        {
            if(gaseste_radacina(nod2) == gaseste_radacina(nod1)) // daca suntem pe op 2 verificam daca cele doua noduri date fac parte din aceeasi multime
                g << "DA" << "\n";
            else
                g << "NU" << "\n";
        }
    }
}
//----------------disjoint----------------------
//-----------------apm-----------------------
void graf :: citire_APM()
{
    f>>n>>m;
    int nod1, nod2, cost;
    vector_costuri.resize(n + 5);
    for(int i = 0; i < m; i++)
    {
        f >> nod1 >> nod2 >> cost;
        vector_costuri[nod1].push_back(make_pair(nod2, cost));
        vector_costuri[nod2].push_back(make_pair(nod1, cost));
    }

}

void graf :: algoritm_APM()
{
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> coada_pri_APM;
    vector<int>tata(n+5, -1); // initial, vectorul de tati va contine doar -1
    vector<bool>ok_APM(n+5, 0); // ce va fi in apm este initializat cu 0
    int suma = 0;
    vector<int>costuri(n+5, 50005); // vector pt retinerea costurilor
    coada_pri_APM.push(make_pair(0, 1));
    costuri[1]=0;
    int nod_curent, cost_curent, nod_vecin;
    while(!coada_pri_APM.empty()) // cat timp avem elemente in priority queue
    {
        nod_curent = coada_pri_APM.top().second; // luma un nod din coada
        coada_pri_APM.pop(); // il extragem
        if(!ok_APM[nod_curent]) // daca modul curent nu este continut in apm
        {
            ok_APM[nod_curent] = 1; // il marcam
            for(auto it: vector_costuri[nod_curent]) //parcurgem vectorul de costuri
            {
                nod_vecin = it.first;
                cost_curent = it.second;
                if(cost_curent < costuri[nod_vecin] && !ok_APM[nod_vecin]) // daca costul curent este mai mic decat valoarea aferenta nodului vecin din vectorul de costuri si nodul vecin nu e continut in apm
                {
                    costuri[nod_vecin] = cost_curent; // valoarea corespuncatoare nodului vecin din vectorul de costuri ia valoarea costului curent
                    coada_pri_APM.push(make_pair(cost_curent, nod_vecin)); // noua pereche cost, nod
                    tata[nod_vecin] = nod_curent; // retinem in vectoul de tati, in dreptul nodului vecin, nodul curent
                }

            }
        }
    }
    for(int i = 1; i <= n; i++)
        suma = suma + costuri[i]; // parcurgem vectoruld e costuri pt a gasi suma
    g << suma << endl;
    g << n - 1 << endl;
    for(int i = 2; i <= n; i++)
        g << i << ' ' << tata[i] << endl;

}
//--------------- apm------------------
graf Gr;
int main()
{
    //BFS
/*
    int N, M, S;
    f>>N>>M>>S;
    Gr.creare_graf(N, M, S);
    Gr.creare_adiacente(M);
    Gr.bfs(S);
    for(int i = 1; i<=N; i++)
        g<< Gr.minime[i]<<" ";
*/
    //DFS
    /*
    int N, M;
    f>>N>>M;
    Gr.creare_graf_Dfs(N, M);
    Gr.creare_adiacente_Dfs(M);
    g<<Gr.comp_conexe(N);
    return 0;
*/

    //CTC
    /*
    int N, M;
    f>>N>>M;
    Gr.creare_graf_ctc(N, M);
    Gr.creare_adiacente_ctc(M);
    Gr.fct_ctc(N);
    */

    //critical connections
    /*
    int N, M;
    f>>N>>M;
    Gr.criticalConnections(N, M);
    */

    //Componente Biconexe
/*
    int N, M;
    f>>N>>M;
    Gr.creare_graf_bic(N, M);
    Gr.creare_adiacente_bic(M);
    Gr.dfs_bic(1);
    Gr.afisare_bic();
*/

    //Havel Hakimi
    /*
    vector<int>v;
    int N;
    int h;
    f>>N;
    for(int i=0; i<N; i++)
        {
            f>>h;
            v.push_back(h);
        }
    if(Gr.HavelOK(N, v)==1)
        g<<"Da";
    else
        g<<"Nu";
    */


    //Sortare Topologica
    /*
    int N, M;
    f>>N>>M;
    Gr.creare_graf_st(N, M);
    Gr.creare_adiacente_st(M);
    Gr.sortare_topologica();
    Gr.afisare_sortare_topo();
    */



    //APM
    /*
    Gr.citire_APM();
	Gr.algoritm_APM();
    */

	//Disjoint
	//Gr.algo_disjuncte();

	//Dijkstra
    //Gr.algo_dijk();

    //Bellman_ford
    //Gr.algo_bellman_ford();

    //MaxFlow
    //Gr.citire_mf();
    //Gr.fct_mf();

    //Roy-Floyd
    //Gr.citire_rf();
    //Gr.fct_royf();

    //Darb
    Gr.citire_darb();
    Gr.fct_darb();
    return 0;
}
