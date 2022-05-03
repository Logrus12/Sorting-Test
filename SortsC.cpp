// SortsC.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
#include<stdio.h>
#include <iostream>
#include <cmath>
#include <ctime>
#include <chrono>
//#include <profileapi.h>
#include <windows.h>
using namespace std;
#define razm 32768
#define srd 0.5;
const int karm = 8192;
int inter[karm + 1];
double step;
double m[razm];
double s[razm];

typedef double TBUF[razm / 16][32];

//class NormalRandom
//{
//    // сохранённое предыдущее значение
//    double prevSample = NAN;
//protected:
//    double Sample()
//    {
//        std::srand(std::time(nullptr));
//        // есть предыдущее значение? возвращаем его
//        if (prevSample != NAN)
//        {
//            double result = prevSample;
//            prevSample = NAN;
//            return result;
//        }
//
//        // нет? вычисляем следующие два
//        // Marsaglia polar method из википедии
//        double u, v, s;
//        do
//        {
//            u = 2 * (double)(std::rand() / RAND_MAX) - 1;
//            v = 2 * (double)(std::rand() / RAND_MAX) - 1; // [-1, 1)
//            s = u * u + v * v;
//        } while (u <= -1 || v <= -1 || s >= 1 || s == 0);
//        double r = std::sqrt(-2 * std::log(s) / s);
//
//        prevSample = r * v;
//        return r * u;
//    }
//};
double generateGaussian(double mean, double stdDev) {
    static double spare;
    static bool hasSpare = false;

    if (hasSpare) {
        hasSpare = false;
        return spare * stdDev + mean;
    }
    else {
        double u, v, s;
        do {
            u = (rand() / ((double)RAND_MAX)) * 2.0 - 1.0;
            v = (rand() / ((double)RAND_MAX)) * 2.0 - 1.0;
            s = u * u + v * v;
        } while (s >= 1.0 || s == 0.0);
        s = sqrt(-2.0 * log(s) / s);
        spare = v * s;
        hasSpare = true;
        return mean + stdDev * u * s;
    }
}
void soed(double s[], int n1, int n2, int dlina, int sn)
{//n1 -начало 1 подмассива, n2 - начало 2-го, dlina их длина,
 //sn - наало записи в массив s
    double temp = s[n1];
    for (int i = 0; i < dlina; i++)//dlina = 2*lim
    {
        s[sn + 2 * i] = temp;
        temp = s[sn + 2 * i + 1];
        s[sn + 2 * i + 1] = m[n2 + i];
    }
}
void soed2(double s[], int n1, int n2, int dlina, int sn)
{
    int i = 0, j = 0, k = 0;
    while ((i < dlina) && (j < dlina))
    {
        if (m[n1 + i] < m[n2 + j])
        {
            s[sn + k] = m[n1 + i]; i++; k++;
        }
        else
        {
            s[sn + k] = m[n2 + j]; j++; k++;
        }
    }
    if (i == dlina)
    {
        while (j < dlina)
        {
            s[sn + k] = m[n2 + j]; j++; k++;
        }
    }
    if (j == dlina)
    {
        while (i < dlina)
        {
            s[sn + k] = m[n1 + i]; i++; k++;
        }
    }

}
void insertsort(double m[], int n1, int n2) {
    int i, j; double temp;
    for (i = n1 + 1; i <= n2; i++)
        for (j = i; (j > n1) && (m[j - 1] > m[j]); j--) // пока j>n1 и элемент j-1 > j, m-массив double
        {
            temp = m[j - 1]; m[j - 1] = m[j]; m[j] = temp;
        }        // меняем местами элементы j и j-1
}
void insertsortb(double m[][32], int n, int n1, int n2) {
    int i, j; double temp;
    for (i = n1; i <= n2; i++)
        for (j = i; j > n1 && m[n][j - 1] > m[n][j]; j--) // пока j>n1 и элемент j-1 > j, m-массив double
        {
            temp = m[n][j - 1]; m[n][j - 1] = m[n][j]; m[n][j] = temp;
        }        // меняем местами элементы j и j-1
}
void Cpy(double s[], double m[]) {
    for (int i = 0; i < razm; i++) m[i] = s[i];
}

void sort2(double s[], int n1, int dlina)
{
    double temp;
    int n = 1, i = 0;
    while (n > 0)
    {
        n = 0;
        for (int j = n1; j < (n1 + dlina - 1 - i); j++)
        {
            if (s[j] > s[j + 1])
            {
                temp = s[j];
                s[j] = s[j + 1];
                s[j + 1] = temp;
                n++;
            }
        }
        i++;
    }
}

void quicksort(double s[], int n1, int n2) {
    int p, i, j; double temp, sp;
    i = n1; j = n2;
    p = ((n1 + n2) / 2); sp = s[p];

    do
    {
        while (s[i] < sp)  i++;
        while (s[j] > sp)  j--;
        if (i <= j) {
            temp = s[i]; s[i] = s[j]; s[j] = temp;
            i++; j--;
        }

    } while (i <= j);

    if (n1 < j)  quicksort(s, n1, j);
    if (i < n2)  quicksort(s, i, n2);
}
int comp(const void* a, const void* b)
{
    return (*(double*)a > *(double*)b) ? 1 : (*(double*)a < *(double*)b) ? -1 : 0;
}
void quicksort2(double s[], int n1, int n2, double sred, double diap) {
    int  i, j; double temp, sp;
    i = n1; j = n2; sp = sred;
    do
    {
        while (s[i] < sp)  i++;
        while (s[j] > sp)  j--;
        if (i <= j) {
            temp = s[i]; s[i] = s[j]; s[j] = temp;
            i++; j--;
        }

    } while (i <= j);
    if (n1 < j)  quicksort2(s, n1, j, sp - diap / 2.0, diap / 2.0);
    if (i < n2)  quicksort2(s, i, n2, sp + diap / 2.0, diap / 2.0);
}
void quicksort4(double s[], int n1, int n2, double sred, double diap) {
    int  i, i1, i2, j, j1, j2; double temp, sp1, sp2, sp3;
    i = n1; j = n2; sp2 = sred; sp1 = sred / 4; sp3 = 3 * sred / 4;
    do
    {
        while (s[i] < sp2)  i++;
        while (s[j] > sp2)  j--;
        if (i <= j) {
            temp = s[i]; s[i] = s[j]; s[j] = temp;
            i++; j--;
        }

    } while (i <= j); i1 = n1; j1 = j;
    do
    {
        while (s[i1] < sp1)  i1++;
        while (s[j1] > sp1)  j1--;
        if (i1 <= j1) {
            temp = s[i1]; s[i1] = s[j1]; s[j1] = temp;
            i++; j--;
        }

    } while (i1 <= j1); i2 = i; j2 = n2;
    do
    {
        while (s[i2] < sp2)  i2++;
        while (s[j2] > sp2)  j2--;
        if (i2 <= j2) {
            temp = s[i2]; s[i2] = s[j2]; s[j2] = temp;
            i++; j--;
        }

    } while (i2 <= j2);

    if (n1 < j1)  quicksort4(s, n1, j1, sp1 - diap / 4.0, diap / 4.0);
    if (i1 < j) quicksort4(s, i1, j, sp2 - diap / 4.0, diap / 4.0);
    if (i < j2) quicksort4(s, i, j2, sp2 + diap / 4.0, diap / 4.0);
    if (i2 < n2)  quicksort4(s, i2, n2, sp2 + diap / 4.0, diap / 4.0);
}
void pocketsort(double s[], int n1, int n2, double min, double max) {
    const int n = razm / 16; int k, sn = 0, inter[n]; double step = (max - min) / n; TBUF buff;
    for (int i = 0; i < n; i++) inter[i] = 0;
    for (int i = n1; i <= n2; i++) {
        k = (int)(m[i] / step); inter[k]++;
        buff[k][i] = m[i];
    }
    for (int i = 0; i < n; i++) {  /* (int j = 31; j > 0; j--) if (s[j] != 0.0)*/ insertsortb(buff, i, 0, inter[i] - 1); }
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < inter[i]; j++) {
            s[sn] = buff[i][j]; sn++;
        }
    }
}
void pocketsort2(double s[], int n1, int n2, double min, double max) {
    const int n = razm / 16; int k, sn = 0, inter[n];/*ind[n];*/ double step = (max - min) / n; double* buf[n];/*int[] ind = new int[n];*/
    for (int i = 0; i < n; i++) { inter[i] = 1; /*ind[i] = 0;*/ }
    for (int i = n1; i <= n2; i++) {
        k = (int)(m[i] / step); inter[k]++;
    }for (int i = 0; i < n; i++) buf[i] = new double[inter[i]]; //cout << "1" << endl;
    for (int i = 0; i < n; i++) inter[i] = 0;
    for (int i = n1; i <= n2; i++) {
        k = (int)(m[i] / step); buf[k][inter[k]] = m[i]; inter[k]++; //cout <<k<< " "<< inter[k] << endl;
    }//cout << "end cycle" << endl;

    for (int i = 0; i < n; i++) { insertsort(buf[i], 0, inter[i] - 1); }//cout << "end cycle 2" << endl;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < inter[i]; j++)
        {
            s[sn] = buf[i][j]; sn++;
        }//cout << "end cycle 3" << endl;
    }//for (int i = 0; i < n; i++) delete buf[i];
}
void flashsort(double m[], int n1, int n2, double min, double max)
{
    inter[karm] = razm;
    int k, k2 /*,sn = 0*/; double step = (max - min) / karm; double temp;
    for (int i = 0; i < karm; i++) inter[i] = 0;//Console.WriteLine(inter[1]);
    for (int i = n1; i <= n2; i++) { k = (int)(m[i] / step); inter[k]++; }
    inter[karm - 1] = n2 - inter[karm - 1] + 1;
    for (int i = karm - 2; i >= 0; i--) inter[i] = inter[i + 1] - inter[i];// Console.WriteLine(inter[0]);
    for (int i = n1; i <= n2; i++)
    {
    Restart: //cout<<"Hello";
        k = (int)(m[i] / step);//Интервал назначения I-го эл
        if ((i < inter[k] || i >= inter[k + 1]))  //не принадлежит своему интервалу
        {
            for (int j = inter[k]; j < inter[k + 1]; j++) {
                k2 = (int)(m[j] / step);
                if (k2 != k) {
                    temp = m[i]; m[i] = m[j]; m[j] = temp; i--; break;
                } //goto Restart;}// break; }//Если J элемент на месте двигаемся дальше,иначе меняем с i-м

            }//goto Restart; //возвращаемся и проверяем интервал нового I-го элемента
        }
    }
    for (int i = 0; i < (karm - 1); i++)
    {
        insertsort(m, inter[i], inter[i + 1] - 1); // Array.Sort(m, inter[i], inter[i + 1] - inter[i] + 1);
    }   insertsort(m, inter[karm - 1], n2);// Array.Sort(m, inter[karm - 1], n2 - inter[karm - 1] +1);
}
void vstavka(double m[], int i1, double* rezerv) {
    int k = (int)(m[i1] / step); int k2; double* rezervold = rezerv;
    for (int j = inter[k]; j < inter[k + 1]; j++) {
        k2 = (int)(m[j] / step);
        if (k2 != k) {
            rezerv = &m[j]; m[j] = m[i1];
            vstavka(m, j, rezerv);
        }
        if (rezervold) { m[i1] = *rezervold; }

        //Ватулин, Ратушня Методы сжатия данных
    }
}

void flashsort2(double m[], int n1, int n2, double min, double max)
{
    inter[karm] = razm;
    int k, k2 /*,sn = 0*/; step = (max - min) / karm; double temp;
    for (int i = 0; i < karm; i++) inter[i] = 0;//Console.WriteLine(inter[1]);
    for (int i = n1; i <= n2; i++) { k = (int)(m[i] / step); inter[k]++; }
    inter[karm - 1] = n2 - inter[karm - 1] + 1;
    for (int i = karm - 2; i >= 0; i--) inter[i] = inter[i + 1] - inter[i];// Console.WriteLine(inter[0]);
    for (int i = n1; i <= n2; i++)
    {
        //k = (int)(m[i] / step);//Интервал назначения I-го эл
        if ((i < inter[k] || i >= inter[k + 1]))  //не принадлежит своему интервалу
        {
            vstavka(m, i, NULL);//cout << i << endl;
        }
    }
    for (int i = 0; i < (karm - 1); i++)
    {
        insertsort(m, inter[i], inter[i + 1] - 1); // Array.Sort(m, inter[i], inter[i + 1] - inter[i] + 1);
    }   insertsort(m, inter[karm - 1], n2);// Array.Sort(m, inter[karm - 1], n2 - inter[karm - 1] +1);
}
void flashsort3(double m[], int n1, int n2, double min, double max)
{
    inter[karm] = razm;
    int j, l, k, k1, k2, sn = 0; step = (max - min) / karm; double temp = 0.0, rezerv = 0.0;
    for (int i = 0; i < karm; i++) inter[i] = 0;//Console.WriteLine(inter[1]);
    for (int i = n1; i <= n2; i++) { k = (int)(m[i] / step); inter[k]++; }
    inter[karm - 1] = n2 - inter[karm - 1] + 1;
    for (int i = karm - 2; i >= 0; i--) inter[i] = inter[i + 1] - inter[i];// Console.WriteLine(inter[0]);
    for (int i = n1; i <= n2; i++)
    {
        k = (int)(m[i] / step); //Интервал назначения I-го эл
        if ((i < inter[k] || i >= inter[k + 1]))  //не принадлежит своему интервалу
        {               //начало цепочки перестановок
            k1 = k;     //карман m[i] эл.
            temp = m[i];//значение, которое мы переносим в нужный карман
            do {
                for (int j = inter[k1]; j < inter[k1 + 1]; j++) {
                    k2 = (int)(m[j] / step);
                    if (k2 != k1) {//rezerv - замещаемое значение, которое мы переносим на след. этапе
                        rezerv = m[j];  m[j] = temp; temp = rezerv; k1 = k2/*(int)(temp / step)*/; break;
                    } //goto Restart;}// break; }//Если J элемент на месте двигаемся дальше
                }
            } while (k1 != k); //m[i] = rezerv;
            //while (true) {
            //    for (j = inter[k1]; j < inter[k1 + 1]; j++) {//цикл поиска позиции вставки и вставки
            //        k2 = (int)(m[j] / step);
            //        if (k2 != k1) {
            //            if (sn != 0) { temp = m[j]; m[j] = m[l]; k1 = (int)(temp / step); break; }
            //            else { temp = m[j]; m[j] = m[i]; k1 = (int)(temp / step); break; }//в первый раз вставляем i-й элемент
            //        }
            //    }
            //if ((int)(temp / step) == k) { m[i] = temp; break; }//если вытесненный элемент принадлежит к тому же карману что i-й, круг замкнулся
            //    for (l = inter[k1]; l < inter[k1 + 1]; l++) {// цикл поиска вставляемого элемента
            //        k2 = (int)(m[l] / step);
            //        if (k2 != k1) { k1 = k2; break; }
            //    }
            //sn++;
            //}sn = 0;
        }
    }
    for (int i = 0; i < (karm - 1); i++)
    {
        insertsort(m, inter[i], inter[i + 1] - 1); // Array.Sort(m, inter[i], inter[i + 1] - inter[i] + 1);
    }   insertsort(m, inter[karm - 1], n2);// Array.Sort(m, inter[karm - 1], n2 - inter[karm - 1] +1);
}
void merg(double m[], int a, int b, int c) {
    double temp; bool sign = m[a + 1] > m[a];
    if (sign) {
        for (int i = a; i < b;) {//левая половина
            for (int j = c; j >= b; j--) {//правая половина
                if ((m[i] < m[j]) || (j == b)) {
                    for (int k = i; k < j - 1; k++) {
                        temp = m[k + 1];
                        m[k + 1] = m[k];
                        m[k] = temp;
                    }
                    if ((j > b) && (b > a)) b--;
                }
                //else continue;
            }
        }
    }
    else {
        for (int i = a; i < b;) {//левая половина
            for (int j = c; j >= b; j--) {//правая половина
                if ((m[i] > m[j]) || (j == b)) {
                    for (int k = i; k < j; k++) {
                        temp = m[k + 1];
                        m[k + 1] = m[k];
                        m[k] = temp;
                    }
                    if ((j > b) && (b > a)) b--;
                }
                //else continue;
            }
        }
    }
}
void pmsort(double m[], int n1, int n2) {
    int i = 0, j = 0, vn = 2; unsigned short vert[razm / 2 + 1]; bool sign, sign2;
    sign = (m[1] > m[0]); vert[0] = n1; for (int i = 1; i < razm / 2 + 1; i++)vert[i] = 0;
    while (vn > 1) {
        vn = 1;
        for (i = n1 + 1; i < n2; i++) {
            sign2 = m[i + 1] > m[i];
            if (sign != sign2) { vert[vn] = i; vn++; sign = !sign; }
        }vert[vn] = n2 + 1;
        //        merg(m, 0, vert[0], vert[1]-1);
        while ((vert[j + 1] > 0) && (j < vn - 1)) {
            merg(m, vert[j], vert[j + 1], vert[j + 2] - 1); j += 2;
        }
        j = 0;
        //        merg(m, 0, vert[i], razm - 1);
        for (int k = 1; k < vn; k++) { vert[k] = 0; }
        //		sign = !sign;
    }
}
void shellsort(double m[], int n1, int n2)
{
    int i, j; int h[] = { 1750,701,301,132,57,23,10,4,1 }; double temp;
    for (int k = 0; k < 9; k++) {
        for (i = n1 + h[k]; i <= n2; i++) {
            for (j = i - h[k]; j >= n1 && m[j] > m[j + h[k]]; j -= h[k]) {
                temp = m[j]; m[j] = m[j + h[k]]; m[j + h[k]] = temp;
            }
        }
    }
}







int main()
{
    int i, j, k; int lim = 16;//текущая длина подмассивов
    int sn;//правая граница заполненности массива s
    LARGE_INTEGER start, start2, finish, finish2, CPS;//double r;
    setlocale(LC_ALL, "Russian");
    std::srand(std::time(NULL));
    QueryPerformanceFrequency(&CPS);
    for (i = 0; i < razm; i++) {
        if (i < razm / 2) { m[i] = (double)(razm - i - 1) / razm; }
        else m[i] = (double)(i - razm / 2) / razm;
        //m[i] = (std::rand() / (double)RAND_MAX);//m[i] = generateGaussian(0.0, 1.0);
    }
    for (i = 0; i < razm; i++) s[i] = m[i];

    //cout<<"Исходный массив : \n";
    //for (i = 0; i < razm; i++)
    //{
    //    cout<<fixed<<m[i]<<"  ";
    //}
    //auto start = std::chrono::high_resolution_clock::now();
    QueryPerformanceCounter(&start);
    //sort2(s,0, razm);
    //quicksort(s, 0, razm - 1);
    //quicksort4(s,0, razm-1,0.5,0.5);
    //pocketsort(s, 0, razm - 1, 0.0, 1.0);
    //flashsort3(m, 0, razm - 1, 0.0, 1.0);
    //insertsort(m, 0, razm - 1);
    //pmsort(m, 0, razm - 1);//merg(m, 0, 16, 31);
    shellsort(m, 0, razm - 1);
    //auto finish = std::chrono::high_resolution_clock::now();
    QueryPerformanceCounter(&finish);
    //std::chrono::duration<double> elapsed = finish - start;
    //cout.precision(4);
     //cout<<"Runtime1 "<<elapsed.count()<<endl;
     //cout<<"Runtime1 "<<(double)(finish2.QuadPart-start2.QuadPart)/CPS.QuadPart<<endl;
    printf("%-9s %e %c", "Runtime1 ", (double)(finish.QuadPart - start.QuadPart) / CPS.QuadPart, '\n');
    //cout.precision(3);
    //cout << "\n Отсортированный массив : \n";
    for (i = 0; i < razm; i++)
    {
        //      cout <<fixed<< m[i] << " ";
        //        printf("%6.3f %c",s[i],' ');
    }
    for (i = 0; i < razm; i++)
    {
        s[i] = m[i];
    }

    QueryPerformanceCounter(&start2);
    for (i = 0; i < razm / lim; i++)
    {
        std::qsort(&m[i * lim], lim, sizeof(double), comp);
    }
    //start = std::chrono::high_resolution_clock::now();
    //QueryPerformanceCounter(&start2);
    k = (int)std::log2(razm / lim);// Количество проходов соединений подмассивов
    for (i = 0; i < k; i++)
    {
        for (j = 0; j < (razm / lim); j++)// j - номер подмассива от 0 до (256/lim - 1)
        {
            if ((j % 2) == 0) { // соед. 1 со 2, 3 c 4 итд.
                sn = lim * j;//границы : 0..lim-1,lim..2lim-1,2lim..3lim-1
                soed2(s, j * lim, (j + 1) * lim, lim, sn);
                //quicksort(m, j * lim, (j + 2) * lim - 1);
                //sort2(s,j * lim, 2 * lim);
                //*if ((i % 2) == 1)*/quicksort(s,j * lim, (j+2) * lim-1);
                //Console.WriteLine("Hello World 2!");
            }

        }
        Cpy(s, m);
        //Console.WriteLine("Hello World!");
        lim *= 2; //начинаем работать с объединенными подмассивами
        //sn = 0;    //обнуляем границу
    }
    //sort2(s,0,lim);
    //finish = std::chrono::high_resolution_clock::now();
    //elapsed = finish - start;
    QueryPerformanceCounter(&finish2);
    //cout<<endl<< "Runtime2 " << elapsed.count()<<endl;
    //cout<<"Runtime2 "<<(double)(finish.QuadPart-start.QuadPart)/CPS.QuadPart<<endl;
    printf("%-9s %e %c", "Runtime2 ", (double)(finish2.QuadPart - start2.QuadPart) / CPS.QuadPart, '\n');
    //sort2(s,0,lim);


    //    cout<<"\n Отсортированный массив : \n";
    //for (i = 0; i < razm; i++)
    //{
    //    cout<<fixed<< s[i]<<" ";
    //}
   // cout<<"Hello World! \n";
    //cout<<"Длина подмассива : "<<lim<<" \n";
    printf("%s %d %c","Длина подмассива : ", lim, '\n');
    scanf_s("%d", k); int ini = int(55);
    //std::cout << "Hello World!\n";
}

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"
